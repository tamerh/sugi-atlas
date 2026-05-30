"""§9 — structural data: per-cohort-gene PDB + AlphaFold availability.
REUSE wrapper over gene §4 + aggregate.

Structure availability drives druggability (a kinase with PDB + co-crystal
is more druggable than one with AlphaFold only)."""
from atlas.section import Section

CHAINS   = (">>uniprot>>pdb", ">>uniprot>>alphafold")  # via gene §4
DATASETS = ("uniprot", "pdb", "alphafold")

def collect(a):
    bundle = {"section": "09_structural_data", "mondo_id": a.mondo_id,
              "_todo": "Implement: fan gene §4 over cohort + aggregate counts"}
    return bundle

SECTION = Section(
    id="9", name="structural_data",
    description=("Per-cohort-gene PDB + AlphaFold structure availability. "
                 "PDB-present / AlphaFold-only / no-structure split feeds §16."),
    needs=("cohort",),
    produces=("per_gene_structure", "pdb_count", "alphafold_only_count",
              "no_structure_count"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
