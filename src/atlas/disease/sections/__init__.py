"""Per-section registry for disease pages — mirrors atlas.gene.sections.

Each s<NN>_*.py module exposes a SECTION (atlas.section.Section). The
orchestrator in atlas.disease.collect imports REGISTRY and runs each section
against a single shared atlas.disease.anchors.DiseaseAnchors record.

§15 (drug_repurposing), §16 (druggability_pyramid), §17 (undrugged_target_profiles)
are render-only — they compute from §1–§14 bundles and don't have their own
collect_fn, so they live in atlas.disease.render, not here."""
from atlas.section import Section
from atlas.disease.sections import (
    s01_disease_ids, s02_gwas_landscape, s03_variant_details,
    s04_mendelian_overlap, s05_genes_proteins, s06_protein_families,
    s07_expression_context, s08_protein_interactions, s09_structural_data,
    s10_drug_targets, s11_bioactivity_enzyme, s12_pharmacogenomics,
    s13_clinical_trials, s14_pathways,
)

_MODULES = (
    s01_disease_ids, s02_gwas_landscape, s03_variant_details,
    s04_mendelian_overlap, s05_genes_proteins, s06_protein_families,
    s07_expression_context, s08_protein_interactions, s09_structural_data,
    s10_drug_targets, s11_bioactivity_enzyme, s12_pharmacogenomics,
    s13_clinical_trials, s14_pathways,
)

REGISTRY: "dict[str, Section]" = {m.SECTION.id: m.SECTION for m in _MODULES}

__all__ = ["Section", "REGISTRY"]
