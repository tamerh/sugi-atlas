#!/usr/bin/env python3
"""Over-representation analysis (ORA) — the deterministic statistics behind the
disease cohort enrichment sections (Reactome pathways first; reusable for GO /
protein-family backgrounds later).

Pure functions, no scipy: the hypergeometric tail is summed in log-space via
math.lgamma (stable at the genome universe N≈11k), and Benjamini-Hochberg FDR is
a sort. Same input → same output, so corpus builds stay byte-identical.

Model (one category, e.g. a pathway):
    N  universe size      — annotated genes (genes with >=1 category membership)
    K  category size       — genes in this category, genome-wide
    n  cohort size          — annotated cohort genes
    k  overlap              — cohort genes in this category
fold-enrichment = (k/n) / (K/N); p = P(X >= k), X ~ Hypergeometric(N, K, n).
"""
import functools
import json
import math
import os


@functools.lru_cache(maxsize=1)
def reactome_background():
    """(universe_n, {pathway_id: genome_wide_gene_count}) from the precomputed
    background (data/background/reactome.json, built by atlas.build_background).
    (0, {}) when absent — callers then fall back to raw-count ranking."""
    here = os.path.dirname(os.path.abspath(__file__))
    for path in ("data/background/reactome.json",
                 os.path.join(here, "..", "..", "data", "background", "reactome.json")):
        try:
            with open(path) as f:
                d = json.load(f)
            return int(d.get("universe_n") or 0), {k: int(v) for k, v in (d.get("sizes") or {}).items()}
        except (OSError, json.JSONDecodeError, ValueError):
            continue
    return 0, {}


def _log_choose(n, k):
    if k < 0 or k > n or n < 0:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k, N, K, n):
    """Upper tail P(X >= k) for X ~ Hypergeometric(population N, successes K,
    draws n). Summed in log-space for numerical stability. Returns 1.0 for k<=0
    (always at least 0 successes) and 0.0 for impossible k."""
    if k <= 0:
        return 1.0
    lo, hi = k, min(K, n)
    if lo > hi:
        return 0.0
    denom = _log_choose(N, n)
    terms = [_log_choose(K, i) + _log_choose(N - K, n - i) - denom
             for i in range(lo, hi + 1)]
    top = max(terms)                       # logsumexp
    p = math.exp(top) * sum(math.exp(t - top) for t in terms)
    return min(max(p, 0.0), 1.0)


def bh_fdr(pvals):
    """Benjamini-Hochberg adjusted p-values, returned in input order. Deterministic
    (ties broken by original index)."""
    m = len(pvals)
    if not m:
        return []
    order = sorted(range(m), key=lambda i: (pvals[i], i))   # ascending p
    adj = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):           # walk largest p -> smallest
        idx = order[rank - 1]
        prev = min(prev, pvals[idx] * m / rank)
        adj[idx] = prev
    return adj


def enrich(items, cohort_n, universe_n):
    """Annotate each category with fold-enrichment + hypergeometric p + BH-FDR.

    items: list of dicts each carrying at least an int `k` (cohort overlap) and an
    int `K` (genome-wide category size). cohort_n = annotated cohort size (n),
    universe_n = annotated universe (N). Returns a NEW list (input order) with
    `fold`, `p_value`, `fdr` added; rows with no/zero K are passed through with
    fold/p/fdr = None (can't be tested — keep the raw count). cohort_n or
    universe_n <= 0 → everything untested (degrades to raw counts upstream)."""
    out = [dict(it) for it in items]
    testable = [it for it in out if (it.get("K") or 0) > 0 and cohort_n and universe_n]
    pvals = [hypergeom_sf(it["k"], universe_n, it["K"], cohort_n) for it in testable]
    fdrs = bh_fdr(pvals)
    for it, p, q in zip(testable, pvals, fdrs):
        exp = cohort_n * it["K"] / universe_n          # expected overlap by chance
        it["fold"] = round(it["k"] / exp, 2) if exp else None
        it["p_value"] = p
        it["fdr"] = q
    for it in out:
        it.setdefault("fold", None)
        it.setdefault("p_value", None)
        it.setdefault("fdr", None)
    return out
