"""Unit tests for atlas.page.jsonld — schema.org Gene JSON-LD emitter.

Exercises: bundle → JSON-LD dict composition, sameAs URL generation,
encodes-protein branching (single / dual / none), missing-field tolerance,
and the two output forms (script tag + sidecar string)."""
import json
import pytest

from atlas.page.jsonld import (
    build_jsonld, same_as_urls, as_script_tag, as_jsonld_string, BASE_URL,
)


# ───── fixtures ───────────────────────────────────────────────────────

def _bundle(**spec):
    """Build a minimal multi-section bundle. Pass any of:
      symbol, name, hgnc_id, location, locus_type, biotype, ensembl_id,
      entrez (list), mim (list), reviewed_uniprot (list), canonical_uniprot,
      aliases (list)"""
    sym = spec.get("symbol", "X")
    return {
        "1": {
            "section": "01_gene_ids", "symbol": sym,
            "hgnc_id": spec.get("hgnc_id"),
            "hgnc": {
                "symbol": sym, "name": spec.get("name"),
                "location": spec.get("location"),
                "locus_type": spec.get("locus_type"),
                "aliases": spec.get("aliases", []),
            },
            "ensembl_id": spec.get("ensembl_id"),
            "ensembl": {"biotype": spec.get("biotype")},
            "entrez": spec.get("entrez", []),
            "mim": spec.get("mim", []),
        },
        "3": {
            "section": "03_protein_ids",
            "reviewed_uniprot": spec.get("reviewed_uniprot", []),
            "canonical_uniprot": spec.get("canonical_uniprot"),
        },
    }

TP53 = _bundle(
    symbol="TP53", name="tumor protein p53", hgnc_id="HGNC:11998",
    location="17p13.1", biotype="protein_coding",
    ensembl_id="ENSG00000141510", entrez=["7157"], mim=["191170"],
    aliases=["p53", "LFS1"],
    reviewed_uniprot=["P04637"], canonical_uniprot="P04637",
)
CDKN2A = _bundle(
    symbol="CDKN2A", name="cyclin dependent kinase inhibitor 2A",
    hgnc_id="HGNC:1787", location="9p21.3", biotype="protein_coding",
    ensembl_id="ENSG00000147889", entrez=["1029"], mim=["600160"],
    reviewed_uniprot=["P42771", "Q8N726"], canonical_uniprot="P42771",
)
MALAT1 = _bundle(
    symbol="MALAT1", name="metastasis associated lung adenocarcinoma transcript 1",
    hgnc_id="HGNC:29665", location="11q13.1", biotype="lncRNA",
    ensembl_id="ENSG00000251562", entrez=["378938"],
    reviewed_uniprot=[],
)


# ───── same_as_urls ───────────────────────────────────────────────────

def test_same_as_full_set_for_protein_coding():
    urls = same_as_urls(TP53)
    # all 5 expected federated identities
    assert any("genenames.org" in u and "HGNC:11998" in u for u in urls)
    assert any("ncbi.nlm.nih.gov/gene/7157" in u for u in urls)
    assert any("uniprot.org/uniprotkb/P04637" in u for u in urls)
    assert any("ensembl.org" in u and "ENSG00000141510" in u for u in urls)
    assert any("omim.org/entry/191170" in u for u in urls)

def test_same_as_omits_missing_sources():
    urls = same_as_urls(MALAT1)
    # MALAT1 has no MIM and no UniProt -> those URLs absent
    assert all("uniprot.org" not in u for u in urls)
    assert all("omim.org" not in u for u in urls)
    # but HGNC, NCBI, Ensembl are present
    assert any("genenames.org" in u for u in urls)
    assert any("ncbi.nlm.nih.gov" in u for u in urls)
    assert any("ensembl.org" in u for u in urls)

def test_same_as_orders_hgnc_first():
    urls = same_as_urls(TP53)
    assert "genenames.org" in urls[0]


# ───── build_jsonld — shape + key fields ─────────────────────────────

def test_jsonld_context_and_type():
    j = build_jsonld(TP53)
    assert j["@context"] == "https://schema.org"
    assert j["@type"] == "Gene"

def test_jsonld_id_url_uses_canonical_atlas_path():
    j = build_jsonld(TP53)
    assert j["@id"] == "https://sugi.bio/atlas/gene/TP53/"
    assert j["url"] == j["@id"]

def test_jsonld_id_respects_custom_base_url():
    j = build_jsonld(TP53, base_url="https://example.org/atlas")
    assert j["@id"] == "https://example.org/atlas/gene/TP53/"

def test_jsonld_name_and_identifier():
    j = build_jsonld(TP53)
    assert j["name"] == "TP53"
    assert j["identifier"] == "HGNC:11998"

def test_jsonld_alternate_name_includes_name_plus_aliases():
    j = build_jsonld(TP53)
    assert "tumor protein p53" in j["alternateName"]
    assert "p53" in j["alternateName"]
    assert "LFS1" in j["alternateName"]

def test_jsonld_description_is_plain_text_no_markdown():
    j = build_jsonld(TP53)
    assert "**" not in j["description"]
    assert j["description"].startswith("TP53")
    assert j["description"].endswith(".")

def test_jsonld_sameAs_present_and_complete():
    j = build_jsonld(TP53)
    assert isinstance(j["sameAs"], list)
    assert len(j["sameAs"]) >= 5


# ───── encodesBioChemEntity branches ──────────────────────────────────

def test_encodes_single_protein_is_dict():
    j = build_jsonld(TP53)
    assert isinstance(j["encodesBioChemEntity"], dict)
    enc = j["encodesBioChemEntity"]
    assert enc["@type"] == "Protein"
    assert enc["name"] == "P04637"
    assert enc["identifier"] == "UniProtKB:P04637"
    assert "uniprot.org/uniprotkb/P04637" in enc["url"]

def test_encodes_dual_product_is_list():
    j = build_jsonld(CDKN2A)
    enc = j["encodesBioChemEntity"]
    assert isinstance(enc, list)
    assert len(enc) == 2
    assert {p["name"] for p in enc} == {"P42771", "Q8N726"}

def test_encodes_absent_for_ncrna():
    j = build_jsonld(MALAT1)
    assert "encodesBioChemEntity" not in j


# ───── layer A: typed, addressable Protein nodes (MOLECULAR_ENRICHMENT) ─────

def _bundle_with_molecular():
    """TP53-shaped bundle carrying §3 ufeatures + §4 structure so the Protein
    node populates its Bioschemas seed."""
    b = _bundle(symbol="TP53", name="tumor protein p53", hgnc_id="HGNC:11998",
                location="17p13.1", biotype="protein_coding",
                reviewed_uniprot=["P04637"], canonical_uniprot="P04637")
    b["3"]["protein_name"] = "Cellular tumor antigen p53"
    b["3"]["ufeatures"] = [
        {"uniprot": "P04637", "type": "active site", "description": "x", "begin": "120", "end": "120"},
        {"uniprot": "P04637", "type": "binding site", "begin": "176", "end": "179"},
        # an ortholog feature on a different accession must NOT leak in
        {"uniprot": "Q00000", "type": "active site", "begin": "9", "end": "9"},
        {"uniprot": "P04637", "type": "sequence conflict", "begin": "1", "end": "1"},  # not a residue type
    ]
    b["4"] = {"section": "04_structure", "reviewed_uniprot": ["P04637"],
              "pdb": [{"id": "1TUP", "method": "X-RAY DIFFRACTION", "resolution": "2.2"}],
              "alphafold": [{"id": "AF-P04637-F1", "uniprot": "P04637", "present": True}]}
    return b


def test_protein_node_has_addressable_id_and_reciprocal_edge():
    j = build_jsonld(_bundle_with_molecular())
    p = j["encodesBioChemEntity"]
    assert p["@id"] == f"{BASE_URL}/gene/TP53/#protein-P04637"
    assert p["isEncodedByBioChemEntity"] == {"@id": f"{BASE_URL}/gene/TP53/"}
    assert p["name"] == "Cellular tumor antigen p53"
    assert p["identifier"] == "UniProtKB:P04637"


def test_protein_sequence_annotations_filtered_and_typed():
    p = build_jsonld(_bundle_with_molecular())["encodesBioChemEntity"]
    anns = p["hasSequenceAnnotation"]
    names = [a["name"] for a in anns]
    assert names == ["active site", "binding site"]          # priority order
    assert "sequence conflict" not in names                  # non-residue dropped
    # the Q00000 ortholog feature did not leak onto P04637
    assert all(a.get("sequenceLocation", {}).get("rangeStart") != 9 for a in anns)
    assert anns[0]["sequenceLocation"] == {"@type": "SequenceRange",
                                           "rangeStart": 120, "rangeEnd": 120}


def test_protein_representations_pdb_and_alphafold():
    p = build_jsonld(_bundle_with_molecular())["encodesBioChemEntity"]
    reps = {r["name"] for r in p["hasRepresentation"]}
    assert "AlphaFold AF-P04637-F1" in reps
    assert "PDB 1TUP" in reps


def test_dual_product_each_gets_own_protein_id():
    enc = build_jsonld(CDKN2A)["encodesBioChemEntity"]
    ids = {p["@id"] for p in enc}
    assert ids == {f"{BASE_URL}/gene/CDKN2A/#protein-P42771",
                   f"{BASE_URL}/gene/CDKN2A/#protein-Q8N726"}
    assert all(p["isEncodedByBioChemEntity"] == {"@id": f"{BASE_URL}/gene/CDKN2A/"}
               for p in enc)


# ───── chromosome / partOf ────────────────────────────────────────────

def test_is_part_of_chromosome():
    j = build_jsonld(TP53)
    assert j["isPartOfBioChemEntity"] == {"@type": "Chromosome", "name": "17p13.1"}

def test_is_part_of_omitted_when_location_missing():
    b = _bundle(symbol="X", biotype="protein_coding", hgnc_id="HGNC:99",
                reviewed_uniprot=["P00001"], canonical_uniprot="P00001")
    j = build_jsonld(b)
    assert "isPartOfBioChemEntity" not in j


# ───── taxonomicRange (human) ─────────────────────────────────────────

def test_taxonomic_range_is_human():
    j = build_jsonld(TP53)
    assert j["taxonomicRange"] == "https://www.ncbi.nlm.nih.gov/taxonomy/9606"


# ───── missing-field tolerance + cleanliness ──────────────────────────

def test_empty_bundle_yields_minimal_but_valid_jsonld():
    j = build_jsonld({})
    assert j["@type"] == "Gene"
    # at minimum we should still have the context + type + canonical URL
    assert j["@id"].endswith("/gene/?/")
    # no None values leak through
    assert None not in j.values()

def test_no_empty_lists_or_none_in_output():
    j = build_jsonld(MALAT1)
    for k, v in j.items():
        assert v not in (None, [], "")


# ───── serialization forms ────────────────────────────────────────────

def test_script_tag_wraps_valid_json():
    j = build_jsonld(TP53)
    tag = as_script_tag(j)
    assert tag.startswith('<script type="application/ld+json">')
    assert tag.endswith("</script>")
    # body parses as the inline-compacted graph (audit #6: over-long arrays are
    # capped for the inline copy; the full graph is the entity.jsonld sidecar).
    from atlas.page.jsonld_inline import compact_for_inline
    body = tag[len('<script type="application/ld+json">'):-len("</script>")].strip()
    assert json.loads(body) == compact_for_inline(j)

def test_sidecar_string_is_parseable_json_with_trailing_newline():
    j = build_jsonld(TP53)
    s = as_jsonld_string(j)
    assert s.endswith("\n")
    assert json.loads(s) == j


# ───── sanity on dual-product encoding ────────────────────────────────

def test_dual_product_sameAs_links_canonical_uniprot_only():
    """sameAs lists only ONE UniProt URL (the canonical) — the dual-product
    detail is encoded via the list-shaped `encodesBioChemEntity`. Avoids
    inflating sameAs into an ambiguous identity claim."""
    urls = same_as_urls(CDKN2A)
    uniprot_urls = [u for u in urls if "uniprot.org" in u]
    assert len(uniprot_urls) == 1
    assert "P42771" in uniprot_urls[0]
