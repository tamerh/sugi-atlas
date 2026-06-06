"""Unit tests for helpers added in the 2026-06 review pass — pure logic, no
network. Each locks in a behavior that was subtle to get right."""
from atlas.disease.cohort import causal_genes
from atlas.render_common import table
from atlas.drug.render import _split_chebi_roles


# ── causal_genes: exact-id GenCC (Definitive/Strong) + OMIM overlap ──────────
def _bundle(disease_gencc=(), omim_genes=()):
    return {"1": {"canonical_name": "x"},
            "4": {"disease_gencc": list(disease_gencc),
                  "omim_genes": list(omim_genes)}}


def test_causal_genes_keeps_definitive_and_strong():
    cg = causal_genes(_bundle(disease_gencc=[
        {"symbol": "CFTR", "classification": "Definitive"},
        {"symbol": "GENE2", "classification": "Strong"},
    ]))
    syms = [s for s, _ in cg]
    assert syms == ["CFTR", "GENE2"]            # Definitive first (rank desc)
    assert cg[0] == ("CFTR", "GenCC Definitive")


def test_causal_genes_drops_moderate_and_limited():
    """The T2D fix — weak per-disease GenCC must NOT be claimed as causal."""
    cg = causal_genes(_bundle(disease_gencc=[
        {"symbol": "WEAK1", "classification": "Moderate"},
        {"symbol": "WEAK2", "classification": "Limited"},
        {"symbol": "WEAK3", "classification": "Supportive"},
    ]))
    assert cg == []


def test_causal_genes_dedups_to_strongest_classification():
    cg = causal_genes(_bundle(disease_gencc=[
        {"symbol": "G", "classification": "Strong"},
        {"symbol": "G", "classification": "Definitive"},
    ]))
    assert cg == [("G", "GenCC Definitive")]


def test_causal_genes_includes_omim_overlap():
    cg = causal_genes(_bundle(omim_genes=[{"symbol": "HBB"}]))
    assert cg == [("HBB", "OMIM Mendelian")]


def test_causal_genes_gencc_wins_over_omim_for_same_gene():
    cg = causal_genes(_bundle(
        disease_gencc=[{"symbol": "G", "classification": "Definitive"}],
        omim_genes=[{"symbol": "G"}]))
    assert cg == [("G", "GenCC Definitive")]    # not duplicated as OMIM


def test_causal_genes_empty_for_polygenic():
    assert causal_genes(_bundle()) == []


# ── table(): empty / all-blank guard ─────────────────────────────────────────
def test_table_empty_rows_renders_nothing():
    assert table(["A", "B"], []) == ""


def test_table_all_blank_rows_renders_nothing():
    assert table(["A", "B"], [(None, ""), ("", None)]) == ""


def test_table_with_data_renders_header_and_row():
    out = table(["A", "B"], [("x", "y")])
    assert "| A | B |" in out and "| x | y |" in out


def test_table_drops_blank_rows_but_keeps_real_ones():
    out = table(["A", "B"], [("x", "y"), (None, None)])
    assert "| x | y |" in out
    assert out.count("\n") == 2                 # header + separator + 1 data row


# ── _split_chebi_roles: pharmacological vs other ─────────────────────────────
def test_split_chebi_roles_separates_non_pharma():
    pharma, other = _split_chebi_roles(
        ["EC 1.1.1.34 (HMG-CoA reductase) inhibitor", "environmental contaminant",
         "xenobiotic", "anticholesteremic drug"])
    assert "environmental contaminant" in other and "xenobiotic" in other
    assert "anticholesteremic drug" in pharma
    assert "EC 1.1.1.34 (HMG-CoA reductase) inhibitor" in pharma


def test_split_chebi_roles_all_other_leaves_pharma_empty():
    pharma, other = _split_chebi_roles(["environmental contaminant", "xenobiotic"])
    assert pharma == [] and len(other) == 2


def test_split_chebi_roles_empty():
    assert _split_chebi_roles([]) == ([], [])


# ── drug alt-name splitting / fragment rejection (gemcitabine IUPAC bug) ──────
from atlas.drug.anchors import _split_synonym, _alt_name_ok


def test_split_synonym_keeps_chemistry_whole():
    iupac = "(2R,3R,4R,5R)-3,3-difluoro-2'-deoxycytidine"
    assert _split_synonym(iupac) == [iupac]                 # parens → never split
    assert _split_synonym("3,3-difluoro-2-deoxycytidine") == \
        ["3,3-difluoro-2-deoxycytidine"]                    # comma-then-digit → kept


def test_split_synonym_splits_brand_code_pairs():
    assert _split_synonym("GLEEVEC,STI-571") == ["GLEEVEC", "STI-571"]
    # digit before comma but a LETTER after → separator, still splits
    assert _split_synonym("STI571,GLEEVEC") == ["STI571", "GLEEVEC"]


def test_alt_name_ok_rejects_chemistry_fragments():
    for frag in ["2", "2'", "2′", "4R", "5R)-3", "(unbalanced", "3,3-"]:
        assert not _alt_name_ok(frag), frag


def test_alt_name_ok_keeps_real_synonyms():
    for ok in ["Gleevec", "STI-571", "LY-188011", "Gemzar", "Imatinib"]:
        assert _alt_name_ok(ok), ok


# ── admission gate 2: non-human / veterinary disease exclusion ───────────────
from atlas.disease.corpus import non_human_ids, NON_HUMAN_ROOT


def test_non_human_ids_subtree_and_flags():
    terms = [
        {"id": NON_HUMAN_ROOT, "name": "non-human animal disease", "parents": [], "non_human": False},
        {"id": "MONDO:9000001", "name": "achondroplasia, cattle",
         "parents": [NON_HUMAN_ROOT], "non_human": True},          # under subtree
        {"id": "MONDO:9000002", "name": "x, cattle (deeper)",
         "parents": ["MONDO:9000001"], "non_human": False},        # descendant of subtree
        {"id": "MONDO:9000003", "name": "y, dog", "parents": ["MONDO:8888"],
         "non_human": True},                                       # flagged, outside subtree
        {"id": "MONDO:9000004", "name": "cystic fibrosis", "parents": ["MONDO:8888"],
         "non_human": False},                                      # human — keep
    ]
    nh = non_human_ids(terms)
    assert NON_HUMAN_ROOT in nh
    assert "MONDO:9000001" in nh and "MONDO:9000002" in nh         # subtree + descendant
    assert "MONDO:9000003" in nh                                   # flagged outside subtree
    assert "MONDO:9000004" not in nh                               # human kept


# ── drug development-status label (approved vs phase-3 seed expansion) ────────
def test_drug_status_label():
    from atlas.page.drug_at_a_glance import _status
    assert _status({"is_fda_approved": True}) == "Approved (max clinical phase 4)"
    assert _status({"max_phase": 4}).startswith("Approved")
    assert _status({"max_phase": 3}) == "Max clinical phase 3 (not approved)"
    assert _status({"max_phase": 2}) == "Max clinical phase 2 (not approved)"
    assert _status({"max_phase": None}) == ""


# ── disease→drug indication section (reverse-indication, tiered phase 3+ / 2) ─
from atlas.disease.render import r_drugs_indicated


def test_r_drugs_indicated_empty_elides():
    assert r_drugs_indicated({}) == ""
    assert r_drugs_indicated({"_indicated_drugs": []}) == ""


def test_r_drugs_indicated_only_phase4_is_indicated_phase3_is_in_trials():
    md = r_drugs_indicated({"_indicated_drugs": [
        {"name": "Rituximab", "url": "/atlas/drug/rituximab/", "max_phase": 4},
        {"name": "Belimumab", "url": "/atlas/drug/belimumab/", "max_phase": 3},
    ]})
    assert "Drugs indicated or in trials for this disease" in md
    # ONLY phase-4 is an approved indication; phase-3 is investigational (in trials)
    assert "1 approved drug" in md and "Approved (phase 4)" in md
    assert "in clinical trials for this disease" in md and "Phase 3" in md
    assert "not* an indication" in md
    assert "/atlas/drug/rituximab/" in md and "/atlas/drug/belimumab/" in md
    approved_part, trials_part = md.split("in clinical trials for this disease")
    assert "Rituximab" in approved_part and "Belimumab" not in approved_part


def test_r_drugs_indicated_tiers_investigational_separately():
    md = r_drugs_indicated({"_indicated_drugs": [
        {"name": "Belimumab", "url": "/atlas/drug/belimumab/", "max_phase": 4},
        {"name": "Atorvastatin", "url": "/atlas/drug/atorvastatin/", "max_phase": 2},
    ]})
    # phase-4 approved table; phase-2 in the investigational block, never "indicated"
    assert "Approved (phase 4)" in md
    approved_part, trials_part = md.split("in clinical trials for this disease")
    assert "Belimumab" in approved_part and "Atorvastatin" not in approved_part
    assert "Atorvastatin" in trials_part


def test_r_drugs_indicated_phase2_only_does_not_claim_indicated():
    md = r_drugs_indicated({"_indicated_drugs": [
        {"name": "Atorvastatin", "url": "/atlas/drug/atorvastatin/", "max_phase": 2},
    ]})
    # no approved indication → explicit disclaimer; shown only as an investigational trial
    assert "No drug is approved" in md
    assert "in clinical trials for this disease" in md and "not* an indication" in md
    assert "Atorvastatin" in md
    assert "Approved (phase 4)" not in md
