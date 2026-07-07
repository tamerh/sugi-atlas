"""Display normalizers in render_common — trial-phase + drug-name casing.
Quality fixes (2026-07-07): PHASE enum → 'Phase N', and multi-word shouting drug
names with digits get title-cased while single-token codes are preserved."""
from atlas.render_common import phase_label, display_name


# ── phase_label: raw CT.gov enum → human, idempotent ─────────────────────────
def test_phase_label_enum_to_human():
    assert phase_label("PHASE1") == "Phase 1"
    assert phase_label("PHASE2") == "Phase 2"
    assert phase_label("PHASE3") == "Phase 3"
    assert phase_label("PHASE4") == "Phase 4"
    assert phase_label("EARLY_PHASE1") == "Early Phase 1"
    assert phase_label("PHASE1/PHASE2") == "Phase 1/2"     # combined


def test_phase_label_idempotent_and_blank():
    assert phase_label("Phase 3") == "Phase 3"             # already-formatted passes through
    assert phase_label("") == "Not specified"
    assert phase_label("NaN") == "Not specified"
    assert phase_label("NA") == "Not specified"


# ── display_name: de-SHOUT names, preserve codes/symbols ─────────────────────
def test_display_name_deshouts_shouting_names():
    assert display_name("FEDRATINIB") == "Fedratinib"
    assert display_name("SUNITINIB") == "Sunitinib"
    # multi-word name WITH a digit is still a name → title-cased (the bug fix)
    assert display_name("INTERFERON ALFA-2B") == "Interferon Alfa-2B"
    assert display_name("PEGINTERFERON ALFA-2B") == "Peginterferon Alfa-2B"


def test_display_name_preserves_codes_symbols_and_mixed():
    assert display_name("ABT199") == "ABT199"     # single-token code (digit, no space)
    assert display_name("K-877") == "K-877"       # single-token code
    assert display_name("TP53") == "TP53"         # gene symbol (upper, no digit-guard needed... )
    assert display_name("Imatinib") == "Imatinib" # already mixed-case, untouched
    assert display_name(None) is None


# ── disease lead: trials clause grammar (no "is a disease AND N trials") ──────
def test_trials_clause_grammar():
    from atlas.page.disease_declarative import _trials_clause
    b13 = {"trial_count": 193}
    # standalone (no preceding evidence) → opens a "with …" (not "and …")
    assert _trials_clause({}, b13, joined=False) == " with 193 registered clinical trials"
    # continuing a preceding "with N cohort genes" list → "and …"
    assert _trials_clause({}, b13, joined=True) == " and 193 registered clinical trials"
    assert _trials_clause({}, {"trial_count": 1}, joined=False) == " with 1 registered clinical trial"
    assert _trials_clause({}, {"trial_count": 0}) == ""
    assert _trials_clause({}, {}) == ""
