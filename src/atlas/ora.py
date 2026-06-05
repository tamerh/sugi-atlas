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
import gzip
import json
import math
import os
from collections import Counter

# Interaction-partner ORA tuning (see the gene-page interactome experiment): the
# standard gene-set-size band excludes umbrella categories (K>500) AND tiny-K
# leaf-pathway noise (K<15); min overlap + FDR filter; rank by fold. A floor on
# annotated partners degrades sparsely-connected genes gracefully.
_MIN_PARTNERS = 10
_K_LO, _K_HI, _MIN_K, _MAX_FDR = 15, 500, 5, 0.01


@functools.lru_cache(maxsize=8)
def background(name):
    """(universe_n, {category_id: genome_wide_gene_count}) from the precomputed
    background data/background/<name>.json (built by atlas.build_background).
    (0, {}) when absent — callers then fall back to raw-count ranking."""
    here = os.path.dirname(os.path.abspath(__file__))
    for path in (f"data/background/{name}.json",
                 os.path.join(here, "..", "..", "data", "background", f"{name}.json")):
        try:
            with open(path) as f:
                d = json.load(f)
            return int(d.get("universe_n") or 0), {k: int(v) for k, v in (d.get("sizes") or {}).items()}
        except (OSError, json.JSONDecodeError, ValueError):
            continue
    return 0, {}


def reactome_background():
    return background("reactome")


@functools.lru_cache(maxsize=1)
def load_membership():
    """(genes, names): genes = {symbol: {'reactome': [ids], 'go': [ids]}},
    names = {id: label}, from data/background/membership.json.gz (atlas.build_
    background). Powers interaction-partner ORA without a per-gene fan at corpus
    build. ({}, {}) when absent."""
    here = os.path.dirname(os.path.abspath(__file__))
    for path in ("data/background/membership.json.gz",
                 os.path.join(here, "..", "..", "data", "background", "membership.json.gz")):
        try:
            with gzip.open(path, "rt") as f:
                d = json.load(f)
            return d.get("genes") or {}, d.get("names") or {}
        except (OSError, json.JSONDecodeError, ValueError):
            continue
    return {}, {}


def interactome_enrichment(partners, kind, top_n=8):
    """ORA of a gene's interaction-partner set vs the genome background. partners:
    gene symbols (self excluded upstream). kind: 'reactome' | 'go'. Returns up to
    top_n enriched categories [{id,name,k,K,fold,fdr}], filtered to the standard
    gene-set size band + min overlap + FDR, ranked by fold. [] when too few
    annotated partners (caller shows a 'not enough partners' note)."""
    genes, names = load_membership()
    counts = Counter()
    annotated = 0
    for p in partners:
        ids = (genes.get(p) or {}).get(kind) or []
        if ids:
            annotated += 1
        for cid in ids:
            counts[cid] += 1
    if annotated < _MIN_PARTNERS:
        return []
    universe_n, sizes = background(kind)
    items = enrich([{"id": cid, "name": names.get(cid), "k": k, "K": sizes.get(cid, 0)}
                    for cid, k in counts.items()], cohort_n=annotated, universe_n=universe_n)
    band = [it for it in items
            if it["fdr"] is not None and it["fdr"] < _MAX_FDR
            and _K_LO <= (it["K"] or 0) <= _K_HI and it["k"] >= _MIN_K]
    band.sort(key=lambda it: (-(it["fold"] or 0.0), it["fdr"], it["id"]))
    return band[:top_n]


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
