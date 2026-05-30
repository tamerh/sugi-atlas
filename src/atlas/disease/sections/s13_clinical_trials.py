"""§13 — clinical trials: mondo → clinical_trials (disease-level, not gene
fanout). Surfaces total trial count, top trials by phase/status, and the
drugs tested across them via clinical_trials → chembl_molecule.

NEW collector — disease anchors directly to trials, no gene cohort fanout."""
from collections import Counter
from atlas.biobtree import map_all, bbmap, map_targets
from atlas.section import Section

CHAINS   = (">>mondo>>clinical_trials",
            ">>mondo>>clinical_trials>>chembl_molecule",
            ">>mesh>>clinical_trials>>chembl_molecule",
            ">>chembl_molecule>>clinical_trials")
DATASETS = ("mondo", "mesh", "clinical_trials", "chembl_molecule")

# Ordering helpers ----------------------------------------------------------
# Trials: PHASE4 > PHASE3 > PHASE2 > PHASE1 > EARLY_PHASE1 > NA/blank.
_PHASE_RANK = {"PHASE4": 4, "PHASE3": 3, "PHASE2": 2,
               "PHASE1/PHASE2": 2, "PHASE2/PHASE3": 3,
               "PHASE1": 1, "EARLY_PHASE1": 0}
_ACTIVE_STATUS = {"RECRUITING", "ACTIVE_NOT_RECRUITING", "ENROLLING_BY_INVITATION",
                  "NOT_YET_RECRUITING", "AVAILABLE"}

def _phase_key(p):
    return _PHASE_RANK.get((p or "").upper(), -1)

def _trial_sort_key(t):
    # phase desc, then active status, then completed
    return (-_phase_key(t.get("phase")),
            0 if (t.get("overall_status") or "").upper() in _ACTIVE_STATUS else 1,
            t.get("id") or "")

def _drug_max_phase(d):
    try:
        return int(d.get("highestDevelopmentPhase") or 0)
    except (ValueError, TypeError):
        return 0

def collect(a):
    # ---- 1. all trials (cheap target dicts; carry phase/status/title) -----
    trials = map_all(a.mondo_id, ">>mondo>>clinical_trials", cap=10)
    trial_count = (a.xref_counts or {}).get("clinical_trials") or len(trials)

    # ---- 2. top 20 trials by phase/status (no entry call needed) ----------
    trials_sorted = sorted(trials, key=_trial_sort_key)
    top_trials = [
        {"id": t.get("id"),
         "title": t.get("brief_title"),
         "phase": t.get("phase"),
         "status": t.get("overall_status"),
         "sponsor": None}                     # not exposed by biobtree
        for t in trials_sorted[:20]
    ]

    phase_counts  = dict(Counter(t["phase"]  or "NA" for t in top_trials))
    status_counts = dict(Counter(t["status"] or "NA" for t in top_trials))

    # ---- 3. trial drugs: mondo + (optional) mesh union -------------------
    drugs = {}                                # molecule_id -> dict
    for d in map_all(a.mondo_id, ">>mondo>>clinical_trials>>chembl_molecule", cap=10):
        mid = d.get("id")
        if mid and mid not in drugs:
            drugs[mid] = {"molecule_id": mid,
                          "name": d.get("name"),
                          "type": d.get("type"),
                          "max_phase": _drug_max_phase(d),  # int 0..4
                          "trial_count": None}
    for mesh in (a.mesh_ids or ()):
        for d in map_all(mesh, ">>mesh>>clinical_trials>>chembl_molecule", cap=10):
            mid = d.get("id")
            if not mid:
                continue
            if mid not in drugs:
                drugs[mid] = {"molecule_id": mid,
                              "name": d.get("name"),
                              "type": d.get("type"),
                              "max_phase": _drug_max_phase(d),
                              "trial_count": None}

    # ---- 4. per-drug trial_count: how many of THIS disease's trials
    # reference each drug. Walk ">>chembl_molecule>>clinical_trials" for every
    # candidate drug (batched + paginated) and intersect with `trial_ids`.
    trial_ids = {t.get("id") for t in trials if t.get("id")}
    if drugs and trial_ids:
        ids = ",".join(drugs.keys())
        per_drug_trials: dict = {mid: set() for mid in drugs}
        page = None
        for _ in range(30):                          # hard cap on pagination
            resp = bbmap(ids, ">>chembl_molecule>>clinical_trials", page)
            for m in (resp.get("mappings") or []):
                mid = (m.get("input") or "").strip()
                if mid not in per_drug_trials:
                    continue
                for tgt in (m.get("targets") or []):
                    nct = tgt.split("|", 1)[0]      # raw "NCT...|title|..."
                    if nct in trial_ids:
                        per_drug_trials[mid].add(nct)
            pg = resp.get("pagination", {}) or {}
            if not pg.get("has_next") or not pg.get("next_token"):
                break
            page = pg.get("next_token")
        for mid, d in drugs.items():
            d["trial_count"] = len(per_drug_trials[mid])

    # Top-30 by max_phase desc, then trial_count desc, then name.
    top_drugs = sorted(drugs.values(),
                       key=lambda x: (-(x["max_phase"] or 0),
                                      -(x["trial_count"] or 0),
                                      x["name"] or ""))[:30]

    return {"section": "13_clinical_trials",
            "mondo_id": a.mondo_id,
            "trial_count": trial_count,
            "top_trials": top_trials,
            "trial_drugs": top_drugs,
            "phase_counts": phase_counts,
            "status_counts": status_counts}

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
