"""Sugi Variant cross-link (gene §6) — the URL builder + the P/LP coverage gate
in r_variants. Mirrors the Sugi Predict pattern: link only genes in Sugi Variant's
pathogenic-gated corpus, so a VUS-only gene never links to a non-existent index."""
from atlas.variant import gene_variants_url, variant_page_url, in_variant_corpus
from atlas.gene.render import r_variants


# Pinned canaries for the SHARED sugislug contract — byte-identical to Sugi
# Variant's own slugs. If these break, the shared slug logic changed and Atlas→
# Variant deep links would 404: coordinate a sugislug version bump, don't just
# re-baseline. (Values verified against Sugi Variant's live index.)
def test_variant_page_url_pinned_slugs():
    cases = [
        ("NM_001100.4(ACTA1):c.925C>G (p.Pro309Ala)", "https://sugi.bio/variant/acta1-p-pro309ala"),
        ("NM_002485.5(NBN):c.1889C>A (p.Ser630Ter)",  "https://sugi.bio/variant/nbn-p-ser630ter"),
        ("NM_000051.4(ATM):c.1009C>T (p.Arg337Cys)",  "https://sugi.bio/variant/atm-p-arg337cys"),
    ]
    for name, expected in cases:
        assert variant_page_url(name) == expected


def test_variant_page_url_uses_transcript_gene():
    # name_gene wins over the passed fallback (the transcript gene is authoritative)
    assert variant_page_url("NM_002485.5(NBN):c.1889C>A (p.Ser630Ter)", gene="WRONG") \
        == "https://sugi.bio/variant/nbn-p-ser630ter"


def test_variant_page_url_none_when_unparseable():
    assert variant_page_url("") is None
    assert variant_page_url("no hgvs here", gene=None) is None


def test_in_variant_corpus_gate():
    for c in ("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic",
              "Conflicting classifications of pathogenicity", "likely pathogenic"):
        assert in_variant_corpus(c)
    for c in ("Uncertain significance", "Benign", "Likely benign", "", None, "risk factor"):
        assert not in_variant_corpus(c)


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
