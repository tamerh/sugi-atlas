"""Mendelian diseases lead #variant-tiers with ClinVar, GWAS tiers after."""
import atlas.disease.cohort as C
from atlas.disease import render as DR

B3 = {"top_variants": [{"rsid": "rs1", "tier": "coding"}], "tier_counts": {"coding": 1},
      "clinvar_total": 1, "clinvar_class_counts": {"Pathogenic": 1},
      "clinvar_variants": [{"id": "7105", "hgvs": "NM_000492.3(CFTR):c.1521_1523del (p.Phe508del)",
                            "gene": "CFTR", "classification": "Pathogenic",
                            "review_status": "practice guideline"}]}


def test_default_gwas_first():
    md = DR.r_variant_details(B3)
    assert md.index("{#top-variants}") < md.index("{#clinvar-variants}")


def test_clinvar_first_for_mendelian():
    md = DR.r_variant_details(B3, clinvar_first=True)
    assert md.index("{#clinvar-variants}") < md.index("{#top-variants}")
    assert md.count("{#clinvar-variants}") == 1


def test_mendelian_disease_flag(monkeypatch):
    monkeypatch.setattr(C, "causal_genes", lambda b: [("CFTR", "GenCC Definitive")])
    assert DR._mendelian_disease({})
    monkeypatch.setattr(C, "causal_genes", lambda b: [("SHROOM4", "GenCC Strong")])
    assert not DR._mendelian_disease({})
