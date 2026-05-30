"""§10 — drug targets: per-cohort-gene ChEMBL targets + approved/phased
molecules. REUSE wrapper over gene §10 subset (chembl_targets + molecules) +
aggregate.

Key derived stat for §16: 'genes with ≥1 approved drug' vs 'undrugged'."""
from atlas.section import Section

CHAINS   = (">>uniprot>>chembl_target",
            ">>chembl_target>>chembl_molecule[highestDevelopmentPhase>=1]")
DATASETS = ("uniprot", "chembl_target", "chembl_molecule")

def collect(a):
    bundle = {"section": "10_drug_targets", "mondo_id": a.mondo_id,
              "_todo": "Implement: fan gene §10 over cohort + drugged/undrugged split"}
    return bundle

SECTION = Section(
    id="10", name="drug_targets",
    description=("Per-cohort-gene ChEMBL targets + phased molecules. "
                 "Approved / Phase ≥3 / Phase ≥1 buckets per gene. "
                 "Drugged-vs-undrugged split."),
    needs=("cohort",),
    produces=("per_gene_drugs", "approved_count", "phased_count",
              "undrugged_count", "drugs"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
