"""§11 — bioactivity & enzyme data: per-cohort-gene ChEMBL assay depth + BRENDA
EC + PubChem activity. REUSE wrapper over gene §3 (brenda) + gene §10
(chembl_assay).

Bioactivity counts are a 'studied-ness' signal that helps rank undrugged
targets for §17 (high assay count = easier starting point)."""
from atlas.section import Section

CHAINS   = (">>uniprot>>chembl_target>>chembl_assay", ">>uniprot>>brenda",
            ">>uniprot>>pubchem_activity")
DATASETS = ("uniprot", "chembl_target", "chembl_assay", "brenda",
            "pubchem_activity")

def collect(a):
    bundle = {"section": "11_bioactivity_enzyme", "mondo_id": a.mondo_id,
              "_todo": "Implement: fan gene §3+§10 over cohort; bioactivity + EC"}
    return bundle

SECTION = Section(
    id="11", name="bioactivity_enzyme",
    description=("Per-cohort-gene ChEMBL assay depth + BRENDA enzyme "
                 "classification + PubChem activity. Bioactivity signal for "
                 "undrugged-target prioritisation."),
    needs=("cohort",),
    produces=("per_gene_bioactivity", "top_studied_genes", "enzyme_genes",
              "undrugged_starting_points"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
