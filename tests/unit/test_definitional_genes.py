"""Definitional gene fallback — when no structured gene route resolves, a rare
disease's causal gene is recovered from the Orphanet clinical definition, which
italicizes gene symbols (e.g. "<i>PIK3CA</i>"). Pure helpers; the live cohort
seeding (anchors.resolve) hits biobtree and is exercised by the corpus build."""
from atlas.disease.anchors import _definitional_genes
from atlas.disease.sections.s05_genes_proteins import _classify


# ── italic-gene extraction ───────────────────────────────────────────────────
def test_definitional_genes_extracts_italic_symbols():
    d = ("A rare PIK3CA-related overgrowth syndrome due to post-zygotic "
         "activating mutations in the <i>PIK3CA</i> gene.")
    assert _definitional_genes(d) == ["PIK3CA"]


def test_definitional_genes_dedups_in_order():
    d = "Caused by <i>BRCA1</i> or <i>BRCA2</i>; <i>BRCA1</i> is most common."
    assert _definitional_genes(d) == ["BRCA1", "BRCA2"]


def test_definitional_genes_ignores_non_italic_and_lowercase():
    # The bare "PIK3CA" (not italicized) and lowercase prose are not candidates;
    # only Orphanet's italic gene-symbol markup is trusted.
    assert _definitional_genes("PIK3CA without markup, post-zygotic <i>mosaic</i>") == []


def test_definitional_genes_empty_and_none():
    assert _definitional_genes("") == []
    assert _definitional_genes(None) == []


# ── evidence partition bucket ────────────────────────────────────────────────
def test_classify_definitional_only():
    ev = {"gwas": False, "gencc": False, "clinvar": False,
          "civic_evidence": False, "definitional": True}
    assert _classify(ev) == "definitional"


def test_classify_structured_evidence_beats_definitional():
    # definitional only fires when nothing structured exists; if some structured
    # flag is present it must NOT be bucketed as definitional.
    ev = {"gwas": True, "gencc": False, "clinvar": False,
          "civic_evidence": False, "definitional": True}
    assert _classify(ev) == "gwas_only"
