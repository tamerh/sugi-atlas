"""Cohort fan-out helper — runs a gene section's `collect_fn(GeneAnchors)`
over every gene in the disease cohort and returns the bundles per gene.

Disease sections §5–§12 and §14 don't define new biobtree chains; they
*aggregate* over the per-gene bundles produced by the existing gene
collectors. This module is the single place where the fan-out happens so the
disease section bodies stay tiny ("call cohort.fan(gene_section, anchors),
then aggregate") and we avoid re-implementing per-gene logic on the disease
side.

Cost: serial today (50 genes × N chains). Threaded fan-out would shrink wall
clock at the cost of contention on biobtree's process-local cache —
introduce only if a section becomes the bottleneck.
"""
from typing import Callable, List, Optional, Dict, Any
from atlas.gene.anchors import Anchors as GeneAnchors


def fan(collect_fn: Callable[[GeneAnchors], dict],
        cohort: "tuple[GeneAnchors, ...]",
        on_error: str = "skip") -> List[dict]:
    """Run `collect_fn(gene_anchors)` for every gene in the cohort.

    Returns one bundle per gene, in the same order as `cohort`. The bundle
    includes the gene's symbol + hgnc_id at minimum (assumed by every gene
    collector), so disease sections can group/sort/filter without a second
    lookup.

    on_error:
      'skip'  — log and drop the gene from the result (default; one bad gene
                shouldn't fail a 50-gene fan-out).
      'raise' — propagate; for debugging.
    """
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
