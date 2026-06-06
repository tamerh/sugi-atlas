"""§5 — clinical trials for the drug (>>chembl_molecule>>clinical_trials,
direct). Total count + true phase/status distribution over all retrieved
trials + top 100 by phase."""
from collections import Counter
from atlas.biobtree import map_all, entry, xref_counts
from atlas.render_common import phase_label
from atlas.section import Section

_PHASE_RANK = {"PHASE4": 4, "PHASE3": 3, "PHASE2/PHASE3": 3, "PHASE2": 2,
               "PHASE1/PHASE2": 2, "PHASE1": 1, "EARLY_PHASE1": 0}
_ACTIVE = {"RECRUITING", "ACTIVE_NOT_RECRUITING", "ENROLLING_BY_INVITATION",
           "NOT_YET_RECRUITING", "AVAILABLE"}


def _phase_key(t):
    return _PHASE_RANK.get((t.get("phase") or "").upper(), -1)


def collect(a):
    trials = map_all(a.chembl_id, ">>chembl_molecule>>clinical_trials", cap=25)
    # True total from the molecule's xref counts — the fetch is capped at 15
    # pages (~1,600), so len(trials) under-reports high-trial drugs (metformin,
    # gemcitabine both hit exactly 1,600). The phase/status distribution + top
    # list below are over the sampled trials; trial_count is the real total.
    total = len(trials)
    try:
        total = (xref_counts(entry(a.chembl_id, "chembl_molecule")) or {}).get(
            "clinical_trials") or total
    except Exception:
        pass
    phase_counts = dict(Counter(phase_label(t.get("phase")) for t in trials))
    status_counts = dict(Counter((t.get("overall_status") or "NA").upper() for t in trials))
    top = sorted(trials,
                 key=lambda t: (-_phase_key(t),
                                0 if (t.get("overall_status") or "").upper() in _ACTIVE else 1,
                                t.get("id") or ""))[:100]
    return {
        "section": "05_clinical_trials",
        "trial_count": total,
        "sampled_trials": len(trials),
        "phase_counts": phase_counts,
        "status_counts": status_counts,
        "top_trials": [{"id": t.get("id"), "title": t.get("brief_title"),
                        "phase": t.get("phase"), "status": t.get("overall_status")}
                       for t in top],
    }


SECTION = Section(
    id="5", name="clinical_trials",
    description=("Clinical trials for the drug (chembl_molecule→clinical_trials); "
                 "count, phase/status distribution, top 100 by phase"),
    needs=("chembl_id",),
    produces=("trial_count", "sampled_trials", "phase_counts", "status_counts", "top_trials"),
    datasets=("chembl_molecule", "clinical_trials"),
    chains=(">>chembl_molecule>>clinical_trials",),
    collect_fn=collect,
)
