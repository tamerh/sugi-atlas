#!/usr/bin/env python3
"""Decide whether a drug→disease indication is an APPROVED indication vs an
investigational trial.

The problem: ChEMBL's per-indication `highestDevelopmentPhase` is the highest
*clinical-trial* phase recorded, and it caps at 3 for many approved drugs —
imatinib is FDA-approved for chronic myeloid leukaemia, yet ChEMBL logs that
indication at phase 3. So "indication phase == 4" under-states real approvals
for ~20% of approved drugs (measured on the v1.1.3 corpus).

We can't recover "approved for THIS disease" in general from this data, and ATC
matching across all 14 therapeutic classes is fragile (a third of the affected
drugs carry no ATC code at all). But it IS reliable for ONCOLOGY: cancer disease
names are unambiguous and ATC L01/L02 = anticancer. v1.12 upgraded such phase-3
rows to "approved"; that over-claimed (see approved_indication), so v1.13 keeps
three tiers instead: approved (phase 4), approved drug with phase-3 evidence for
the use (late), and investigational trials.
"""

# ATC 1st/2nd level for anticancer therapy: L01 (antineoplastic agents),
# L02 (endocrine therapy). L03 (immunostimulants) / L04 (immunosuppressants)
# are excluded — they treat autoimmune/transplant, not cancer.
_ONCOLOGY_ATC = ("L01", "L02")

# Cancer disease-name keywords. Distinctive enough to avoid benign "-oma"
# collisions (glaucoma, trachoma, granuloma, adenoma, lipoma are NOT matched).
_CANCER_KW = (
    "cancer", "carcinoma", "adenocarcinoma", "neoplasm", "neoplasia",
    "leukemia", "leukaemia", "lymphoma", "sarcoma", "melanoma", "glioma",
    "glioblastoma", "blastoma", "myeloma", "mesothelioma", "malignan",
    "tumor", "tumour", "hodgkin", "myelodysplastic", "myeloproliferative",
    "carcinoid", "metastatic", "oncolog",
)


def is_oncology_drug(atc_codes):
    return any((c or "").upper().startswith(_ONCOLOGY_ATC) for c in (atc_codes or []))


def _int_phase(mp):
    try:
        return int(float(mp)) if mp is not None else None
    except (TypeError, ValueError):
        return None


def has_phase4_trial(b5):
    """True if the molecule has a registered Phase-4 (post-marketing) trial.

    A Phase-4 trial only exists for an approved drug, so this recovers approvals
    that ChEMBL's `max_phase` under-states — notably non-oncology oligonucleotides
    like inclisiran/Leqvio (FDA-approved Dec 2021, yet ChEMBL caps it at phase 3),
    which carry no PubChem is_fda_approved flag. Reads §5's phase distribution
    (`phase_counts`, keyed by the raw uppercased trial phase, e.g. "PHASE4")."""
    pc = (b5 or {}).get("phase_counts") or {}
    return any(v and str(k).upper() == "PHASE4" for k, v in pc.items())


def molecule_approved(b1, b5=None):
    """Single source of truth for "is this an approved drug?" at the molecule
    level (status bullet, lead class clause, "not yet approved" callout).

    Approved when PubChem says so, or ChEMBL max_phase is 4, or the drug is
    late-stage (max_phase >= 3) AND has a registered Phase-4 trial — the last
    branch catches drugs ChEMBL under-phases. Conservative: it only ever upgrades
    to approved. NOTE: this is molecule-level; per-indication tiering is
    approved_indication (phase 4) / late_stage_use (phase 3 of an approved drug)."""
    if b1.get("is_fda_approved") or _int_phase(b1.get("max_phase")) == 4:
        return True
    mp = _int_phase(b1.get("max_phase"))
    return mp is not None and mp >= 3 and has_phase4_trial(b5)


def is_cancer_disease(name):
    n = (name or "").lower()
    return any(kw in n for kw in _CANCER_KW)


def approved_indication(atc_codes, phase, disease_name, molecule_approved):
    """True if this drug→disease indication is recorded at approval phase (4).

    v1.13: the oncology upgrade (approved L01/L02 drug × phase-3 cancer row →
    "approved") is gone. ChEMBL carries no label source per indication, and the
    upgrade turned every phase-3 cancer trial of an approved cancer drug into an
    approval — palbociclib "approved for 9" (label: HR+/HER2− breast cancer),
    cisplatin 100; 28% of all approved rows. Phase-3 uses of an approved drug are
    the separate `late` tier (late_stage_use) — reported as phase-3 evidence,
    never as approval. Signature kept for callers."""
    return (phase or 0) >= 4


def late_stage_use(phase, molecule_approved):
    """True for a phase-3 indication of an APPROVED molecule — ChEMBL often caps
    real labelled uses at phase 3 (imatinib→CML, abacavir→HIV), but the same tier
    holds phase-3 trials of other uses, so it's reported as "approved drug,
    phase-3 evidence for this use", not as an approval."""
    return (phase or 0) == 3 and bool(molecule_approved)


def named_indications(b4):
    """Indications with a real disease name (raw EFO/MeSH/MP ids dropped), deduped
    by (name, MONDO) — the source repeats a disease across xrefs."""
    from atlas.render_common import is_ontology_id
    named, seen = [], set()
    for i in (b4 or {}).get("indications") or []:
        nm = (i.get("name") or "").strip()
        if not nm or is_ontology_id(nm):
            continue
        key = (nm.lower(), i.get("mondo_id") or "")
        if key in seen:
            continue
        seen.add(key)
        named.append(i)
    return named


def indication_tiers(b4):
    """(approved, late, trials) — the drug's named indications deduped by (name,
    MONDO): approved = phase 4; late = phase 3 of an approved molecule (the
    `late` flag set in s04); trials = the remaining phase 2–3 (the disease-side
    index's floor). The single source for the indications table, the lead's
    counts, the FAQ and JSON-LD `treats`, so they can't disagree."""
    named = named_indications(b4)
    approved = [i for i in named if i.get("approved")]
    late = [i for i in named if not i.get("approved") and i.get("late")]
    trials = [i for i in named if not i.get("approved") and not i.get("late")
              and 2 <= (i.get("max_phase") or 0) <= 3]
    return approved, late, trials
