"""§6 — protein-family classification: cohort proteins classified by druggable
families (Kinase / GPCR / Ion channel / Nuclear receptor / Protease / TF / ...)
via InterPro. REUSE wrapper over gene §3 + custom classifier.

Druggable vs difficult split drives §16 (druggability pyramid)."""
from atlas.section import Section

CHAINS   = (">>uniprot>>interpro",)  # reused via gene §3
DATASETS = ("uniprot", "interpro", "pfam")

def collect(a):
    bundle = {"section": "06_protein_families", "mondo_id": a.mondo_id,
              "_todo": "Implement: classify by InterPro families; druggable/difficult"}
    return bundle

SECTION = Section(
    id="6", name="protein_families",
    description=("Cohort proteins classified by druggable family (Kinase, GPCR, "
                 "Ion channel, Nuclear receptor, Protease, Phosphatase, Enzyme, "
                 "TF, Scaffold) via InterPro. Druggable vs difficult split."),
    needs=("cohort",),
    produces=("family_assignments", "family_counts", "druggable_count",
              "difficult_count"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
