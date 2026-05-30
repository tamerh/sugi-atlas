"""Per-section registry. Each `s<NN>_*` module exposes a SECTION (atlas.gene.sections.base.Section).

The orchestrator in atlas.gene.collect imports REGISTRY and runs each section
against a single shared atlas.gene.anchors.Anchors record.

Adding a section = drop a new s<NN>_*.py here that exports SECTION."""
from atlas.gene.sections.base import Section
from atlas.gene.sections import (
    s01_gene_ids, s02_transcripts, s03_protein_ids, s04_structure,
    s05_orthologs, s06_variants, s07_pathways, s08_interactions,
    s09_tf_regulation, s10_drugs, s11_expression, s12_diseases,
)

_MODULES = (
    s01_gene_ids, s02_transcripts, s03_protein_ids, s04_structure,
    s05_orthologs, s06_variants, s07_pathways, s08_interactions,
    s09_tf_regulation, s10_drugs, s11_expression, s12_diseases,
)

REGISTRY: "dict[str, Section]" = {m.SECTION.id: m.SECTION for m in _MODULES}

__all__ = ["Section", "REGISTRY"]
