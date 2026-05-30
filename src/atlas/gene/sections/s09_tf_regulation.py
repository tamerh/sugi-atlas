"""§9 — tf_regulation: CollecTRI downstream targets (direction-filtered) +
upstream regulators, JASPAR motifs. is_tf inferred from having targets or motifs."""
from atlas.biobtree import map_all
from atlas.gene.sections.base import Section

CHAINS = (
    '>>hgnc>>collectri[tf_gene=="<symbol>"]',
    '>>hgnc>>collectri[target_gene=="<symbol>"]',
    ">>uniprot>>jaspar",
)
DATASETS = ("collectri", "jaspar", "hgnc", "uniprot")

def collect(a):
    bundle = {"section": "09_tf_regulation", "symbol": a.symbol}
    # direction-filtered + fully paginated -> high-degree TFs (TP53: 1207 CollecTRI
    # rows) get complete downstream targets, not a 100-row truncation.
    down = map_all(a.hgnc_id, f'>>hgnc>>collectri[tf_gene=="{a.symbol}"]')
    up   = map_all(a.hgnc_id, f'>>hgnc>>collectri[target_gene=="{a.symbol}"]')
    bundle["downstream_targets"] = [{"target": r.get("target_gene"),
                                     "regulation": r.get("regulation")} for r in down]
    bundle["downstream_count"] = len(down)
    bundle["upstream_regulators"] = [{"regulator": r.get("tf_gene"),
                                      "regulation": r.get("regulation")} for r in up]
    bundle["jaspar_motifs"] = [{"id": t["id"], "name": t.get("name"),
                                "class": t.get("class"), "family": t.get("family")}
                               for t in (map_all(a.canonical_uniprot, ">>uniprot>>jaspar")
                                         if a.canonical_uniprot else [])]
    bundle["is_transcription_factor"] = bool(down or bundle["jaspar_motifs"])
    return bundle

SECTION = Section(
    id="9", name="tf_regulation",
    description="Transcription regulation — CollecTRI downstream targets + upstream regulators (direction-filtered) and JASPAR motifs",
    needs=("hgnc_id", "canonical_uniprot", "symbol"),
    produces=("downstream_targets", "upstream_regulators", "jaspar_motifs", "is_transcription_factor"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
