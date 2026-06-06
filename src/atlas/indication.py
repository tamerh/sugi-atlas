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


def is_cancer_disease(name):
    n = (name or "").lower()
    return any(kw in n for kw in _CANCER_KW)


def approved_indication(atc_codes, phase, disease_name):
    """True if this drug→disease indication counts as an APPROVED indication.

    phase 4 is always approved. A phase-3 indication is upgraded to approved
    only for an anticancer drug (ATC L01/L02) against a cancer — the one place
    ChEMBL's phase-3 reliably means an approved labelled use rather than an
    ongoing trial. Phase <3 is never approved.
    """
    p = phase or 0
    if p >= 4:
        return True
    if p >= 3 and is_oncology_drug(atc_codes) and is_cancer_disease(disease_name):
        return True
    return False
