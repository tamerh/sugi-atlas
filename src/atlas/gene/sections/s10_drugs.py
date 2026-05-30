"""§10 — drugs: ChEMBL targets + phased molecules, PharmGKB, BindingDB sample,
clinical trials via the DISEASE route (gene → MONDO → trials)."""
from atlas.biobtree import map_all
from atlas.gene.sections.base import Section

CHAINS = (
    ">>uniprot>>chembl_target",
    '>>chembl_target>>chembl_molecule[highestDevelopmentPhase>=1]',
    ">>hgnc>>pharmgkb_gene",
    ">>uniprot>>bindingdb",
    ">>hgnc>>gencc>>mondo>>clinical_trials",
    ">>hgnc>>clinvar>>mondo>>clinical_trials",
)
DATASETS = ("chembl_target", "chembl_molecule", "pharmgkb_gene", "bindingdb",
            "clinical_trials", "mondo", "gencc", "clinvar", "uniprot", "hgnc")

def _phase(d):
    return int(d) if (d or "").isdigit() else 0

def collect(a):
    bundle = {"section": "10_drugs", "symbol": a.symbol}
    uni = a.canonical_uniprot

    targets = map_all(uni, ">>uniprot>>chembl_target") if uni else []
    bundle["chembl_targets"] = [{"id": t["id"], "title": t.get("title"), "type": t.get("type")}
                                for t in targets]

    # Phased drugs only (highestDevelopmentPhase>=1) — the raw chembl_molecule edge
    # returns thousands of ID-ordered screening compounds; filter + rank for the real drugs.
    drugs = {}
    for t in targets:
        for m in map_all(t["id"], '>>chembl_target>>chembl_molecule[highestDevelopmentPhase>=1]', cap=10):
            drugs[m["id"]] = {"id": m["id"], "name": m.get("name"), "type": m.get("type"),
                              "phase": m.get("highestDevelopmentPhase")}
    bundle["molecules"] = sorted(drugs.values(), key=lambda d: _phase(d["phase"]), reverse=True)
    bundle["molecule_count"] = len(drugs)

    bundle["pharmgkb"] = [{"id": t["id"], "vip": t.get("is_vip"),
                           "cpic_guideline": t.get("has_cpic_guideline")}
                          for t in map_all(a.hgnc_id, ">>hgnc>>pharmgkb_gene")]

    bd = map_all(uni, ">>uniprot>>bindingdb", cap=2) if uni else []
    bundle["bindingdb_sample"] = [{"ligand": t.get("ligand_name"), "ki": t.get("ki"),
                                   "ic50": t.get("ic50")} for t in bd[:30]]
    bundle["bindingdb_sampled"] = len(bd)

    # Clinical trials via the DISEASE route — chembl_molecule>>clinical_trials would
    # pollute with off-target drugs (ChEMBL target↔molecule is bioactivity-based).
    # Disease-level is biobtree's intended pattern.
    trials = {}
    for chain in (">>hgnc>>gencc>>mondo>>clinical_trials",
                  ">>hgnc>>clinvar>>mondo>>clinical_trials"):
        for t in map_all(a.hgnc_id, chain, cap=2):
            trials[t["id"]] = {"id": t["id"], "title": t.get("brief_title"),
                               "phase": t.get("phase"), "status": t.get("overall_status")}
    def _ph(t):
        p = (t.get("phase") or "").upper().replace("PHASE", "").strip()
        return int(p) if p.isdigit() else 0
    bundle["disease_trials"] = sorted(trials.values(), key=_ph, reverse=True)
    bundle["disease_trial_count"] = len(trials)
    bundle["is_drug_target"] = bool(targets)
    return bundle

SECTION = Section(
    id="10", name="drugs",
    description="ChEMBL targets + phased targeting molecules, PharmGKB, BindingDB affinities, clinical trials via disease route",
    needs=("hgnc_id", "canonical_uniprot"),
    produces=("chembl_targets", "molecules", "pharmgkb", "bindingdb_sample",
              "disease_trials", "is_drug_target"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
