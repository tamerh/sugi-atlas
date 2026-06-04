"""P2/P3 frontmatter key-facts (identifier, alt_names, tldr, section_defaults)."""
from atlas.page.meta_facts import entity_facts, _clean_aliases


def test_clean_aliases_dedups_case_without_splitting():
    # _clean_aliases no longer splits on comma (that shattered IUPAC names);
    # splitting moved to the drug anchor. Here it only strips + case-dedups,
    # keeping the first-seen form. A comma-joined value is left intact.
    out = _clean_aliases(["GLEEVEC,STI-571", "Imatinib", "IMATINIB", "imatinib", ""])
    assert out == ["GLEEVEC,STI-571", "Imatinib"]      # NOT split; first-form, deduped


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
    # alt_names arrive already split from the drug anchor (which can tell a
    # separator comma from a chemistry comma); entity_facts just normalizes +
    # dedups and must NOT re-split (that shattered IUPAC names — gemcitabine).
    f = entity_facts("drug", {"1": {"chembl_id": "CHEMBL941",
                                    "alt_names": ["Gleevec", "STI-571", "Gleevec"]}})
    assert f["identifier"] == "CHEMBL941"
    assert f["alt_names"] == ["Gleevec", "STI-571"]            # deduped, not re-split
    assert f["section_defaults"] == {"summary": "open", "indications": "open"}


def test_clean_aliases_does_not_split_chemical_commas():
    """Regression: a comma-joined IUPAC fragment must NOT be split by the
    frontmatter normalizer (that's what produced gemcitabine's "2'"/"5R)-3")."""
    from atlas.page.meta_facts import _clean_aliases
    out = _clean_aliases(["(2R,3R,4R,5R)-3,3-difluoro-2'-deoxycytidine"])
    assert out == ["(2R,3R,4R,5R)-3,3-difluoro-2'-deoxycytidine"]
