"""P2/P3 frontmatter key-facts (identifier, alt_names, tldr, section_defaults)."""
from atlas.page.meta_facts import entity_facts, _clean_aliases


def test_clean_aliases_splits_commas_and_dedups_case():
    out = _clean_aliases(["GLEEVEC,STI-571", "Imatinib", "IMATINIB", "imatinib", ""])
    assert out == ["GLEEVEC", "STI-571", "Imatinib"]   # split, first-form kept, dedup


def test_gene_facts():
    b = {"1": {"symbol": "EGFR", "hgnc_id": "HGNC:3236",
               "hgnc": {"aliases": ["ERBB1", "ERRP"]}},
         "3": {"protein_name": "Epidermal growth factor receptor",
               "canonical_uniprot": "P00533"},
         "10": {"is_drug_target": True}}
    f = entity_facts("gene", b)
    assert f["identifier"] == "EGFR"                   # typed id (not the slug)
    assert f["alt_names"] == ["ERBB1", "ERRP"]
    assert f["tldr"][0] == "Encodes Epidermal growth factor receptor (UniProt P00533)"
    assert f["section_defaults"] == {"summary": "open", "drugs": "open"}  # druggable


def test_disease_identifier_is_mondo():
    f = entity_facts("disease", {"1": {"mondo_id": "MONDO:0007959",
                                       "synonyms": ["brain medulloblastoma"]}})
    assert f["identifier"] == "MONDO:0007959"
    assert f["alt_names"] == ["brain medulloblastoma"]
    assert f["section_defaults"]["summary"] == "open"


def test_drug_identifier_is_chembl_and_section_default_indications():
    f = entity_facts("drug", {"1": {"chembl_id": "CHEMBL941",
                                    "alt_names": ["Gleevec,STI-571"]}})
    assert f["identifier"] == "CHEMBL941"
    assert f["alt_names"] == ["Gleevec", "STI-571"]
    assert f["section_defaults"] == {"summary": "open", "indications": "open"}
