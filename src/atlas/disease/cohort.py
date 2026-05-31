"""Cohort fan-out helper — runs a gene section's `collect_fn(GeneAnchors)`
over every gene in the disease cohort and returns the bundles per gene.

Disease sections §5–§12 and §14 don't define new biobtree chains; they
*aggregate* over the per-gene bundles produced by the existing gene
collectors. This module is the single place where the fan-out happens so the
disease section bodies stay tiny ("call cohort.fan(gene_section, anchors),
then aggregate") and we avoid re-implementing per-gene logic on the disease
side.

## Memoization

Multiple disease sections fan the same gene collector over the same cohort:
gene §3 (protein_ids) is fanned by disease §5, §6, §11 — three identical
fans, ~6 biobtree calls × 50 genes = ~900 redundant calls. Same for gene §10
(drug_targets) fanned by disease §10 + §11.

`fan()` memoizes per `(collect_fn identity, cohort identity)`. Within one
`collect_all(disease)` invocation, every gene collector runs at most once
over the cohort. `_reset_cache()` is called by `atlas.disease.collect.collect_all`
between diseases in batch runs so memory doesn't grow unbounded.

Pure perf — zero content change.

Cost: serial today (50 genes × N chains). Threaded fan-out would shrink wall
clock further at the cost of contention on biobtree's process-local cache —
introduce only if the cached path is still a bottleneck.
"""
from typing import Callable, List, Optional, Dict, Any
from atlas.gene.anchors import Anchors as GeneAnchors


# (id(collect_fn), id(cohort_tuple)) → list[bundle]. Tuples are immutable +
# per-anchor-instance unique, so id() is stable until _reset_cache().
_FAN_CACHE: Dict[tuple, List[dict]] = {}


def _reset_cache() -> None:
    """Called by disease.collect.collect_all between diseases in a batch run."""
    _FAN_CACHE.clear()


def fan(collect_fn: Callable[[GeneAnchors], dict],
        cohort: "tuple[GeneAnchors, ...]",
        on_error: str = "skip") -> List[dict]:
    """Run `collect_fn(gene_anchors)` for every gene in the cohort.

    Returns one bundle per gene, in the same order as `cohort`. The bundle
    includes the gene's symbol + hgnc_id at minimum (assumed by every gene
    collector), so disease sections can group/sort/filter without a second
    lookup.

    Results are memoized per (collect_fn, cohort) — see module docstring.

    on_error:
      'skip'  — log and drop the gene from the result (default; one bad gene
                shouldn't fail a 50-gene fan-out).
      'raise' — propagate; for debugging.
    """
    if not cohort:
        return []
    cache_key = (id(collect_fn), id(cohort))
    cached = _FAN_CACHE.get(cache_key)
    if cached is not None:
        return cached

    out: List[dict] = []
    for ga in cohort:
        try:
            b = collect_fn(ga)
            # Defensive — make sure symbol/hgnc_id always present in result.
            b.setdefault("symbol", ga.symbol)
            b.setdefault("hgnc_id", ga.hgnc_id)
            out.append(b)
        except Exception as e:
            if on_error == "raise":
                raise
            # Keep going — partial coverage beats a hard failure on rare HGNC.
            out.append({"symbol": ga.symbol, "hgnc_id": ga.hgnc_id,
                        "_error": str(e)})
    _FAN_CACHE[cache_key] = out
    return out


def by_symbol(bundles: List[dict]) -> Dict[str, dict]:
    """Index a fan-out result by gene symbol for fast lookup in render code."""
    return {b["symbol"]: b for b in bundles if "symbol" in b}


def filter_evidence(cohort: "tuple[GeneAnchors, ...]",
                    cohort_evidence: Dict[str, Dict[str, bool]],
                    route: str) -> "tuple[GeneAnchors, ...]":
    """Subset the cohort to genes hit by a specific evidence route
    ('gwas' / 'gencc' / 'clinvar' / 'civic_evidence'). Used by §4 Mendelian
    overlap and §16 druggability slicing."""
    return tuple(g for g in cohort
                 if (cohort_evidence.get(g.hgnc_id) or {}).get(route))
