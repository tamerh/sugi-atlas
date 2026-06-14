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
names are unambiguous and ATC L01/L02 = anticancer. So we upgrade a phase-3
indication to "approved" only when the drug is an antineoplastic AND the disease
is a cancer. Everything else is reported honestly by phase (see the renderers'
caveat) — never falsely "approved", never falsely "in trials only".
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
    to approved. NOTE: this is molecule-level; per-indication "approved for THIS
    disease" tiering stays oncology-only (see approved_indication)."""
    if b1.get("is_fda_approved") or _int_phase(b1.get("max_phase")) == 4:
        return True
    mp = _int_phase(b1.get("max_phase"))
    return mp is not None and mp >= 3 and has_phase4_trial(b5)


def is_cancer_disease(name):
    n = (name or "").lower()
    return any(kw in n for kw in _CANCER_KW)


def approved_indication(atc_codes, phase, disease_name, molecule_approved):
    """True if this drug→disease indication counts as an APPROVED indication.

    A phase-4 indication is always approved. A phase-3 indication is upgraded to
    approved only when ALL of: the molecule is FDA-approved (PubChem
    is_fda_approved, or ChEMBL molecule max_phase 4), it is an anticancer drug
    (ATC L01/L02), and the disease is a cancer. The molecule-approval gate is
    load-bearing: without it an *experimental* anticancer drug in phase-3 cancer
    trials would be falsely marked approved. Phase <3 is never approved.
    """
    p = phase or 0
    if p >= 4:
        return True
    if (p >= 3 and molecule_approved and is_oncology_drug(atc_codes)
            and is_cancer_disease(disease_name)):
        return True
    return False
