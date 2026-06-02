"""§13 — clinical trials: mondo → clinical_trials (disease-level, not gene
fanout). Surfaces total trial count, top trials by phase/status, and the
drugs tested across them via clinical_trials → chembl_molecule.

NEW collector — disease anchors directly to trials, no gene cohort fanout."""
from collections import Counter
from atlas.biobtree import map_all, bbmap, map_targets, entry
from atlas.civic import aggregate_predictive
from atlas.render_common import phase_label
from atlas.section import Section

CHAINS   = (">>mondo>>clinical_trials",
            ">>mondo>>clinical_trials>>chembl_molecule",
            ">>mesh>>clinical_trials>>chembl_molecule",
            ">>chembl_molecule>>clinical_trials",
            ">>mondo>>civic_evidence")
DATASETS = ("mondo", "mesh", "clinical_trials", "chembl_molecule", "civic_evidence")

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


def _name_tokens(a):
    """Disease name tokens for trial-title validation: canonical name + Mondo
    synonyms, ≥4 chars (drops ambiguous 2–3 char acronyms that match anything)."""
    names = [a.canonical_name] + list(a.synonyms or ())
    return [n.lower() for n in names if n and len(n.strip()) >= 4]


def _title_matches(title, tokens):
    t = (title or "").lower()
    return bool(t) and any(tok in t for tok in tokens)


def collect(a):
    # ---- 1. trials — VALIDATED by title. biobtree's mondo→clinical_trials edge
    # is contaminated: it links trials whose actual conditions don't match the
    # disease (Vici syndrome → 1,156 glaucoma/cataract studies; see
    # BIOBTREE_ISSUES). The raw map + xref count are therefore untrustworthy, so
    # we keep only trials whose brief_title names the disease (or a synonym).
    # brief_title is in the map projection → no per-trial fetch. Drugs/counts
    # below all derive from this validated set (a contaminated drug ends up in 0
    # validated trials → trial_count 0 → dropped).
    raw_trials = map_all(a.mondo_id, ">>mondo>>clinical_trials", cap=10)
    tokens = _name_tokens(a)
    trials = [t for t in raw_trials if _title_matches(t.get("brief_title"), tokens)] if tokens else []
    trial_count = len(trials)
    trial_count_raw = (a.xref_counts or {}).get("clinical_trials") or len(raw_trials)

    # ---- 2. true phase/status distribution over ALL trials (not the
    # top-20 sample, which is biased toward PHASE4 by sort).
    phase_counts  = dict(Counter(phase_label(t.get("phase")) for t in trials))
    status_counts = dict(Counter((t.get("overall_status") or "NA").upper() for t in trials))

    # ---- 3. top 20 trials by phase/status (no entry call needed) ----------
    trials_sorted = sorted(trials, key=_trial_sort_key)
    top_trials = [
        {"id": t.get("id"),
         "title": t.get("brief_title"),
         "phase": t.get("phase"),
         "status": t.get("overall_status"),
         "sponsor": None}                     # not exposed by biobtree
        for t in trials_sorted[:40]
    ]

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

    # ---- 5. parent/child salt-form dedupe via biobtree's `childs` field.
    # ChEMBL treats salt/anhydrous/etc. forms as distinct molecules (CHEMBL92
    # = "DOCETAXEL ANHYDROUS" parent, CHEMBL3545252 = "DOCETAXEL" child).
    # Both can appear in trial_drugs because both are referenced by trials.
    # We fold every child into its parent, summing trial_count and keeping
    # the parent's name. Build child→parent map only for the top candidates
    # (bounds entry calls to ~30 per disease).
    candidate_ids = sorted(drugs.keys(),
                           key=lambda mid: -(drugs[mid].get("trial_count") or 0))[:60]
    child_to_parent = {}
    for mid in candidate_ids:
        try:
            en = entry(mid, "chembl_molecule")
            attrs = (((en.get("Attributes") or {}).get("Chembl") or {}))
            mol = attrs.get("molecule") if isinstance(attrs.get("molecule"), dict) else attrs
            for child_id in (mol.get("childs") or []):
                child_to_parent[child_id] = mid
        except Exception:
            continue
    # Apply: replace child entries with parent (sum trial_count).
    folded = {}
    for mid, d in drugs.items():
        parent_id = child_to_parent.get(mid, mid)
        if parent_id in folded:
            folded[parent_id]["trial_count"] = (folded[parent_id]["trial_count"] or 0) \
                                               + (d.get("trial_count") or 0)
            folded[parent_id]["max_phase"] = max(folded[parent_id]["max_phase"] or 0,
                                                  d.get("max_phase") or 0)
        else:
            # Use the parent's metadata if we have it loaded, else the child's.
            if parent_id != mid and parent_id in drugs:
                folded[parent_id] = dict(drugs[parent_id])
                folded[parent_id]["trial_count"] = (drugs[parent_id].get("trial_count") or 0) \
                                                   + (d.get("trial_count") or 0)
                folded[parent_id]["max_phase"] = max(drugs[parent_id].get("max_phase") or 0,
                                                     d.get("max_phase") or 0)
            else:
                folded[parent_id] = dict(d)
                folded[parent_id]["molecule_id"] = parent_id

    # Top-30 by max_phase desc, then trial_count desc, then name. Keep ONLY
    # drugs referenced by ≥1 VALIDATED trial (trial_count>0) — a drug seen only
    # in title-mismatched (contaminated) trials would otherwise pollute the
    # lead's "top interventions" (e.g. cataract drugs for Vici syndrome).
    top_drugs = sorted((d for d in folded.values() if (d.get("trial_count") or 0) > 0),
                       key=lambda x: (-(x["max_phase"] or 0),
                                      -(x["trial_count"] or 0),
                                      x["name"] or ""))[:30]

    # ---- 6. CIViC precision-subtype map — the drug × variant × indication
    # triple at disease level. mondo→civic_evidence yields predictive
    # associations stratified by molecular subtype (e.g. NSCLC → EGFR T790M →
    # Osimertinib, ALK fusion → Alectinib, KRAS G12C → Sotorasib). Same
    # aggregation as gene §10 (dedup + rank by CIViC evidence level); here the
    # `profile` column reads as the molecular subtype of the disease. Empty
    # for non-cancer diseases so the block elides cleanly.
    civic = map_all(a.mondo_id, ">>mondo>>civic_evidence", cap=15)
    civic_ranked, civic_stats = aggregate_predictive(civic, top=30)

    return {"section": "13_clinical_trials",
            "mondo_id": a.mondo_id,
            "trial_count": trial_count,
            "trial_count_raw": trial_count_raw,
            "top_trials": top_trials,
            "trial_drugs": top_drugs,
            "phase_counts": phase_counts,
            "status_counts": status_counts,
            "civic_evidence": civic_ranked,
            "civic_evidence_total": civic_stats["evidence_total"],
            "civic_predictive_total": civic_stats["predictive_total"],
            "civic_association_total": civic_stats["association_total"],
            "civic_evidence_type_counts": civic_stats["evidence_type_counts"]}

SECTION = Section(
    id="13", name="clinical_trials",
    description=("Disease-level clinical trials from mondo→clinical_trials + "
                 "drugs tested (clinical_trials → chembl_molecule). Phase / "
                 "status / lead-drug attrs. Plus CIViC precision-subtype map "
                 "(drug × molecular subtype × indication via mondo→civic_evidence)."),
    needs=("mondo_id", "mesh_ids", "xref_counts"),
    produces=("trial_count", "trial_count_raw", "top_trials", "trial_drugs", "phase_counts",
              "status_counts", "civic_evidence", "civic_evidence_total",
              "civic_predictive_total", "civic_association_total",
              "civic_evidence_type_counts"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
