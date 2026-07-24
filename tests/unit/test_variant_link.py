"""Sugi Variant cross-link (gene §6) — the URL builder + the P/LP coverage gate
in r_variants. Mirrors the Sugi Predict pattern: link only genes in Sugi Variant's
pathogenic-gated corpus, so a VUS-only gene never links to a non-existent index."""
from atlas.variant import gene_variants_url
from atlas.gene.render import r_variants


def test_gene_variants_url():
    assert gene_variants_url("PTEN") == "https://sugi.bio/variant/gene/PTEN"
    assert gene_variants_url("  TP53 ") == "https://sugi.bio/variant/gene/TP53"
    assert gene_variants_url("") is None
    assert gene_variants_url(None) is None


def _bundle(patho=0, lp=0):
    return {"symbol": "PTEN", "clinvar_total": 10,
            "clinvar_breakdown": {"Pathogenic": patho, "Likely pathogenic": lp,
                                  "Uncertain significance": 5, "Likely benign": 0, "Benign": 0},
            "top_pathogenic": [], "top_spliceai": [], "top_alphamissense": [], "dbsnp_sample": []}


def test_crosslink_emitted_for_pathogenic_gene():
    md = r_variants(_bundle(patho=3))
    assert "### Per-variant reference — Sugi Variant {#sugi-variant}" in md
    assert "https://sugi.bio/variant/gene/PTEN" in md
    assert "predictor-disagreement QC" in md


def test_crosslink_emitted_for_lp_only_gene():
    md = r_variants(_bundle(patho=0, lp=2))     # LP but no Pathogenic → still covered
    assert "sugi.bio/variant/gene/PTEN" in md


def test_crosslink_elided_for_vus_only_gene():
    md = r_variants(_bundle(patho=0, lp=0))     # not in Sugi Variant's corpus → no link
    assert "Sugi Variant" not in md
    assert "sugi-variant" not in md
