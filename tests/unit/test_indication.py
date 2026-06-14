"""Regression guard for the drug→disease "approved vs in-trials" tiering.

The bug this catches: ChEMBL logs many approved indications at phase 3 (imatinib
is FDA-approved for CML/GIST, yet ChEMBL records those at phase 3). A naive
"indication phase == 4 ⇒ approved" tiering therefore reports imatinib with ZERO
approved indications and files its real approvals under "in clinical trials".
A corpus sanity check caught it; these tests would have caught it earlier.
"""
from atlas.indication import (approved_indication, is_oncology_drug,
                              is_cancer_disease, has_phase4_trial,
                              molecule_approved)


def test_phase4_is_always_approved():
    # phase-4 indication is approved regardless of the molecule-approval gate
    assert approved_indication([], 4, "anything", False)
    assert approved_indication(["B01AC06"], 4, "stroke disorder", True)


def test_phase_le2_is_never_approved():
    assert not approved_indication(["L01EA01"], 2, "chronic myeloid leukemia", True)
    assert not approved_indication([], 2, "asthma", True)


def test_imatinib_phase3_cancer_is_approved():
    # THE regression: APPROVED anticancer drug (ATC L01) + cancer + phase 3 ⇒ approved.
    assert approved_indication(["L01EA01"], 3, "chronic myeloid leukemia", True)
    assert approved_indication(["L01EA01"], 3, "gastrointestinal stromal tumor", True)


def test_experimental_anticancer_phase3_cancer_is_NOT_approved():
    # molecule-approval gate (PubChem is_fda_approved): an UNAPPROVED anticancer
    # drug in phase-3 cancer trials must NOT be marked approved.
    assert not approved_indication(["L01EA01"], 3, "chronic myeloid leukemia", False)


def test_aspirin_phase3_cancer_stays_investigational():
    # antithrombotic/analgesic (B01/N02), NOT anticancer — its phase-3 cancer
    # chemoprevention trials must NOT be promoted to "approved" (even though the
    # molecule itself is FDA-approved).
    assert not approved_indication(["B01AC06", "N02BA01"], 3, "squamous cell carcinoma", True)
    assert not approved_indication(["B01AC06"], 3, "colorectal neoplasm", True)


def test_oncology_drug_phase3_noncancer_is_not_approved():
    # imatinib's phase-3 stroke / pneumonia trials are not approved indications.
    assert not approved_indication(["L01EA01"], 3, "stroke disorder", True)
    assert not approved_indication(["L01EA01"], 3, "pneumonia", True)


def test_has_phase4_trial_reads_phase_counts():
    assert has_phase4_trial({"phase_counts": {"PHASE4": 8, "PHASE3": 18}})
    assert not has_phase4_trial({"phase_counts": {"PHASE3": 18, "PHASE2": 4}})
    assert not has_phase4_trial({})
    assert not has_phase4_trial(None)


def test_molecule_approved_inclisiran_regression():
    # THE regression: inclisiran (Leqvio) is FDA-approved (Dec 2021) but ChEMBL
    # caps it at max_phase 3 and PubChem carries no is_fda_approved flag. A
    # registered Phase-4 trial (post-marketing) recovers the approval.
    b1 = {"max_phase": 3, "is_fda_approved": None}
    b5 = {"phase_counts": {"PHASE4": 8, "PHASE3": 18}}
    assert molecule_approved(b1, b5)
    # without the Phase-4 signal it correctly stays not-approved
    assert not molecule_approved(b1, {"phase_counts": {"PHASE3": 18}})


def test_molecule_approved_direct_signals():
    assert molecule_approved({"max_phase": 4})
    assert molecule_approved({"is_fda_approved": True, "max_phase": 1})
    assert not molecule_approved({"max_phase": 2}, {"phase_counts": {"PHASE2": 3}})
    # a stray Phase-4 trial on an early molecule is NOT enough (guard: max_phase>=3)
    assert not molecule_approved({"max_phase": 2}, {"phase_counts": {"PHASE4": 1}})


def test_helpers_are_specific():
    assert is_oncology_drug(["L01EA01"]) and is_oncology_drug(["L02BX03"])
    assert not is_oncology_drug(["B01AC06"]) and not is_oncology_drug([])
    assert is_cancer_disease("chronic myeloid leukemia")
    assert is_cancer_disease("breast carcinoma") and is_cancer_disease("glioblastoma")
    # benign "-oma" words and non-cancers must not match
    assert not is_cancer_disease("glaucoma")
    assert not is_cancer_disease("asthma") and not is_cancer_disease("trachoma")
