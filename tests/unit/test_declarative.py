"""Unit tests for atlas.page.declarative.

Pure function — exercise the gene-class mapping, the optional-field handling,
and the protein-coding vs ncRNA branch."""
import pytest

from atlas.page.declarative import declarative_sentence, _gene_class, BIOTYPE_HUMAN


# ───── fixtures (synthetic bundles — no biobtree) ─────────────────────

def _b1(symbol, **hgnc_extra):
    """Build a minimal §1 (gene_ids) bundle slice."""
    b = {
        "section": "01_gene_ids",
        "symbol": symbol,
        "hgnc_id": hgnc_extra.pop("hgnc_id", None),
        "hgnc": {
            "symbol": symbol,
            "name": hgnc_extra.pop("name", None),
            "location": hgnc_extra.pop("location", None),
            "locus_type": hgnc_extra.pop("locus_type", None),
        },
        "ensembl_id": hgnc_extra.pop("ensembl_id", None),
        "ensembl": {"biotype": hgnc_extra.pop("biotype", None)},
    }
    return b

def _b3(canonical_uniprot=None, reviewed=()):
    return {
        "section": "03_protein_ids",
        "reviewed_uniprot": list(reviewed),
        "canonical_uniprot": canonical_uniprot,
    }

TP53_BUNDLE = {
    "1": _b1("TP53", name="tumor protein p53", hgnc_id="HGNC:11998",
             location="17p13.1", locus_type="gene with protein product",
             biotype="protein_coding"),
    "3": _b3(canonical_uniprot="P04637", reviewed=("P04637",)),
}

MALAT1_BUNDLE = {
    "1": _b1("MALAT1", name="metastasis associated lung adenocarcinoma transcript 1",
             hgnc_id="HGNC:29665", location="11q13.1",
             locus_type="RNA, long non-coding", biotype="lncRNA"),
    "3": _b3(reviewed=()),
}

CDKN2A_BUNDLE = {  # dual-product gene
    "1": _b1("CDKN2A", name="cyclin dependent kinase inhibitor 2A",
             hgnc_id="HGNC:1787", location="9p21.3",
             locus_type="gene with protein product", biotype="protein_coding"),
    "3": _b3(canonical_uniprot="P42771", reviewed=("P42771", "Q8N726")),
}

MIR21_BUNDLE = {
    "1": _b1("MIR21", name="microRNA 21", hgnc_id="HGNC:31586",
             location="17q23.1", locus_type="RNA, micro", biotype="miRNA"),
    "3": _b3(),
}


# ───── _gene_class ────────────────────────────────────────────────────

@pytest.mark.parametrize("biotype,expected", [
    ("protein_coding", "protein-coding gene"),
    ("lncRNA", "long non-coding RNA gene"),
    ("miRNA", "microRNA gene"),
    ("snoRNA", "small nucleolar RNA gene"),
    ("pseudogene", "pseudogene"),
    ("Mt_tRNA", "mitochondrial tRNA gene"),
    ("IG_V_gene", "immunoglobulin V gene"),
    ("TR_V_gene", "T-cell receptor V gene"),
])
def test_gene_class_via_biotype(biotype, expected):
    b1 = _b1("X", biotype=biotype)
    assert _gene_class(b1) == expected

def test_gene_class_locus_fallback_protein_coding():
    # biotype missing, locus_type names it
    b1 = _b1("X", locus_type="gene with protein product")
    assert _gene_class(b1) == "protein-coding gene"

def test_gene_class_locus_fallback_lncrna():
    b1 = _b1("X", locus_type="RNA, long non-coding")
    assert _gene_class(b1) == "long non-coding RNA gene"

def test_gene_class_locus_fallback_pseudogene():
    b1 = _b1("X", locus_type="pseudogene")
    assert _gene_class(b1) == "pseudogene"

def test_gene_class_default_to_gene():
    # neither biotype nor locus_type give us anything actionable
    b1 = _b1("X")
    assert _gene_class(b1) == "gene"

def test_biotype_human_table_covers_common_classes():
    # smoke that the table at least has the headline cases
    for k in ("protein_coding", "lncRNA", "miRNA", "pseudogene"):
        assert k in BIOTYPE_HUMAN


# ───── declarative_sentence — full sentence shape ─────────────────────

def test_protein_coding_full_sentence():
    s = declarative_sentence(TP53_BUNDLE)
    assert s == ("**TP53** (tumor protein p53, HGNC:11998) is a protein-coding "
                 "gene on chromosome 17p13.1, encoding the reviewed UniProt "
                 "protein P04637.")

def test_lncrna_sentence_omits_protein_clause():
    s = declarative_sentence(MALAT1_BUNDLE)
    assert s.startswith("**MALAT1** (metastasis associated lung adenocarcinoma transcript 1, HGNC:29665)")
    assert "long non-coding RNA gene" in s
    assert "chromosome 11q13.1" in s
    assert "UniProt" not in s   # no protein clause for ncRNA

def test_microrna_sentence():
    s = declarative_sentence(MIR21_BUNDLE)
    assert "microRNA gene" in s
    assert "UniProt" not in s

def test_dual_product_uses_canonical_only():
    """For CDKN2A (P42771 + Q8N726), the lead cites the canonical only —
    the dual-product detail belongs in §3 tables, not the headline."""
    s = declarative_sentence(CDKN2A_BUNDLE)
    assert "P42771" in s
    assert "Q8N726" not in s

def test_sentence_ends_in_single_period():
    for b in (TP53_BUNDLE, MALAT1_BUNDLE, CDKN2A_BUNDLE, MIR21_BUNDLE):
        s = declarative_sentence(b)
        assert s.endswith(".")
        assert not s.endswith("..")

def test_sentence_is_single_line():
    for b in (TP53_BUNDLE, MALAT1_BUNDLE, CDKN2A_BUNDLE, MIR21_BUNDLE):
        assert "\n" not in declarative_sentence(b)


# ───── optional-field tolerance ───────────────────────────────────────

def test_missing_name_drops_from_parens():
    b = {"1": _b1("KRAS", hgnc_id="HGNC:6407", location="12p12.1",
                  biotype="protein_coding"),
         "3": _b3(canonical_uniprot="P01116")}
    s = declarative_sentence(b)
    # no name -> parens contain only HGNC id
    assert "**KRAS** (HGNC:6407)" in s
    assert "P01116" in s

def test_missing_hgnc_id_drops_from_parens():
    b = {"1": _b1("X", name="protein x", biotype="protein_coding"),
         "3": _b3()}
    s = declarative_sentence(b)
    assert "**X** (protein x)" in s
    assert "HGNC" not in s

def test_missing_everything_optional_still_well_formed():
    b = {"1": _b1("UNKNOWNGENE"), "3": _b3()}
    s = declarative_sentence(b)
    # no parens, no chromosome, default class — but well-formed sentence
    assert s == "**UNKNOWNGENE** is a gene."

def test_missing_location_omits_chromosome_clause():
    b = {"1": _b1("TP53", name="tumor protein p53", hgnc_id="HGNC:11998",
                  biotype="protein_coding"),
         "3": _b3(canonical_uniprot="P04637")}
    s = declarative_sentence(b)
    assert "chromosome" not in s
    assert "encoding the reviewed UniProt protein P04637" in s

def test_missing_canonical_uniprot_omits_protein_clause():
    b = {"1": _b1("X", name="protein x", hgnc_id="HGNC:99",
                  location="1p1", biotype="protein_coding"),
         "3": _b3()}  # no canonical
    s = declarative_sentence(b)
    assert "UniProt" not in s
    assert s.endswith("chromosome 1p1.")

def test_name_equal_to_symbol_drops_from_parens():
    # if HGNC `name` happens to equal `symbol`, don't redundantly include it
    b = {"1": _b1("TP53", name="TP53", hgnc_id="HGNC:11998",
                  location="17p13.1", biotype="protein_coding"),
         "3": _b3(canonical_uniprot="P04637")}
    s = declarative_sentence(b)
    assert "**TP53** (HGNC:11998)" in s
    # name "TP53" should NOT appear in parens since it equals symbol
    assert "(TP53," not in s


def test_missing_bundle_sections_default_to_unknown_gene():
    """Even a near-empty bundle should not crash."""
    s = declarative_sentence({})
    assert s == "**?** is a gene."
