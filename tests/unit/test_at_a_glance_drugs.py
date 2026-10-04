"""At-a-glance drug headline bullets (v1.11.8): approved MOA drugs + PharmGKB
guidelines on gene pages; approved indicated drugs on disease pages; the gene
"no CIViC evidence" Notable line no longer fires on pharmacogenes / Mendelian genes."""
from atlas.page import at_a_glance as G
from atlas.page import disease_at_a_glance as DG


def _gene(**b10):
    return {"6": {"clinvar_total": 500, "clinvar_breakdown": {}}, "10": b10, "12": {}}


def test_gene_approved_moa_bullet():
    md = G.at_a_glance(_gene(moa_drugs=[
        {"name": "GEFITINIB", "phase": "4"}, {"name": "ERLOTINIB", "phase": "4.0"},
        {"name": "OSIMERTINIB", "phase": 4}, {"name": "AFATINIB", "phase": "4"},
        {"name": "TRIALDRUG", "phase": "2"}, {"name": None, "phase": "4"}]))
    assert "**Approved drugs (ChEMBL mechanism of action):** 4 — Gefitinib, Erlotinib, Osimertinib, +1 more" in md
    assert "TRIALDRUG" not in md and "Trialdrug" not in md


def test_gene_pgx_bullet_and_notable_gating():
    pgx = G.at_a_glance(_gene(pharmgkb_guideline=[
        {"chemical_names": "clopidogrel"}, {"chemical_names": "citalopram, escitalopram"},
        {"chemical_names": "Clopidogrel"}]))
    assert "**Pharmacogenomics (PharmGKB):** 3 CPIC/DPWG dosing guidelines (e.g. clopidogrel, citalopram, escitalopram)" in pgx
    assert "Notable" not in pgx                                     # pharmacogene
    plain = G.at_a_glance(_gene())
    assert "no curated precision-oncology (CIViC) evidence yet" in plain
    mendel = G.at_a_glance({**_gene(), "12": {"gencc": [{"disease": "cystic fibrosis",
                                                         "classification": "Definitive"}]}})
    assert "Notable" not in mendel                                  # Mendelian gene


def test_disease_approved_indicated_bullet():
    b = {"1": {}, "13": {"trial_count": 40}, "10": {"phased_count": 0},
         "_indicated_drugs": [{"name": "IVACAFTOR", "approved": True},
                              {"name": "trial-x", "approved": False}]}
    md = DG.at_a_glance(b)
    assert "**Approved drugs (ChEMBL indications):** 1 — Ivacaftor" in md
    assert "no approved drug yet" not in md                         # no contradiction
    b["_indicated_drugs"] = []
    assert "no approved drug yet" in DG.at_a_glance(b)
