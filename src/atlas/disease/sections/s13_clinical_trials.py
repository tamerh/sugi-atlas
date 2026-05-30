"""§13 — clinical trials: mondo → clinical_trials (disease-level, not gene
fanout). Surfaces total trial count, top trials by phase/status, and the
drugs tested across them via clinical_trials → chembl_molecule.

NEW collector — disease anchors directly to trials, no gene cohort fanout."""
from atlas.section import Section

CHAINS   = (">>mondo>>clinical_trials",
            ">>mondo>>clinical_trials>>chembl_molecule",
            ">>mesh>>clinical_trials>>chembl_molecule")
DATASETS = ("mondo", "mesh", "clinical_trials", "chembl_molecule")

def collect(a):
    bundle = {"section": "13_clinical_trials", "mondo_id": a.mondo_id,
              "_todo": "Implement: mondo->clinical_trials + per-trial attrs"}
    return bundle

SECTION = Section(
    id="13", name="clinical_trials",
    description=("Disease-level clinical trials from mondo→clinical_trials + "
                 "drugs tested (clinical_trials → chembl_molecule). Phase / "
                 "status / lead-drug attrs."),
    needs=("mondo_id", "mesh_ids", "xref_counts"),
    produces=("trial_count", "top_trials", "trial_drugs", "phase_counts",
              "status_counts"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
