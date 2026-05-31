"""Target fan-out helper — the drug analog of disease's cohort.fan().

Drug §8 (target_pathways) reuses the gene §7 pathways collector over the
drug's *target* genes, rather than re-implementing pathway logic. The cohort
here is the drug's GtoPdb-curated primary targets (~1-9 genes), so the fan is
cheap — no cap needed (disease fans up to 50 genes; drugs far fewer).

§7 (related_molecules) works directly off the target UniProt accessions (no
GeneAnchors needed), so it doesn't use this helper.
"""
from typing import Callable, List
from atlas.gene.anchors import Anchors as GeneAnchors, resolve as resolve_gene_anchors


def target_gene_anchors(anchors) -> "List[GeneAnchors]":
    """Resolve one GeneAnchors per unique primary-target gene symbol. Targets
    without a resolved gene symbol (e.g. protein complexes) are skipped."""
    seen, out = set(), []
    for t in anchors.targets:
        sym = t.gene_symbol
        if not sym or sym in seen:
            continue
        seen.add(sym)
        try:
            out.append(resolve_gene_anchors(sym))
        except Exception:
            continue
    return out


def over_targets(collect_fn: Callable[[GeneAnchors], dict],
                 gene_anchors: "List[GeneAnchors]",
                 on_error: str = "skip") -> List[dict]:
    """Run a gene section's `collect_fn(GeneAnchors)` over each target gene.
    Returns one bundle per gene (symbol/hgnc_id always present). One bad gene
    is skipped rather than failing the whole fan."""
    out: List[dict] = []
    for ga in gene_anchors:
        try:
            b = collect_fn(ga)
            b.setdefault("symbol", ga.symbol)
            b.setdefault("hgnc_id", ga.hgnc_id)
            out.append(b)
        except Exception as e:
            if on_error == "raise":
                raise
            out.append({"symbol": ga.symbol, "hgnc_id": ga.hgnc_id, "_error": str(e)})
    return out
