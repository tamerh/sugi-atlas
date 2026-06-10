"""Unit tests for atlas.page.provenance — per-page upstream-source trail."""
import json

import pytest

from atlas.page.provenance import (
    UPSTREAM, build_provenance, as_provenance_string, _section_provenance,
)
from atlas.gene.sections import REGISTRY


# ───── fixtures ───────────────────────────────────────────────────────

def _bundle(symbol="TP53", hgnc_id="HGNC:11998", ensembl_id="ENSG00000141510",
            canonical_uniprot="P04637", canonical_transcript="ENST00000269305"):
    return {
        "1": {"section": "01_gene_ids", "symbol": symbol, "hgnc_id": hgnc_id,
              "ensembl_id": ensembl_id, "hgnc": {"name": "tumor protein p53",
              "location": "17p13.1"}, "entrez": ["7157"], "mim": ["191170"]},
        "3": {"section": "03_protein_ids", "reviewed_uniprot": [canonical_uniprot],
              "canonical_uniprot": canonical_uniprot},
        "6": {"section": "06_variants", "canonical_transcript": canonical_transcript},
    }

META = {
    "generated_at": "2026-05-30T12:00:00+00:00",
    "atlas_version": "0.1.0",
    "biobtree_version": "unknown",
}


# ───── UPSTREAM table coverage ────────────────────────────────────────

def test_upstream_has_entries_for_every_section_dataset():
    """Every dataset referenced by any section's metadata should have an
    upstream-source entry — otherwise an agent looking at the page can't
    follow the trail. Catches drift when a new dataset is added to a
    section but its upstream isn't registered."""
    missing = set()
    for sec in REGISTRY.values():
        for ds in sec.datasets:
            if ds not in UPSTREAM:
                missing.add(ds)
    assert not missing, (
        f"UPSTREAM map missing entries for: {sorted(missing)}. "
        "Add each to atlas.page.provenance.UPSTREAM with (display name, source URL).")

def test_upstream_urls_are_http_https():
    for name, url in UPSTREAM.values():
        assert url.startswith(("http://", "https://")), f"{name}: {url}"


# ───── _section_provenance ────────────────────────────────────────────

def test_section_provenance_shape():
    sec = REGISTRY["1"]   # gene_ids
    p = _section_provenance(sec)
    assert p["id"] == "1"
    assert p["name"] == "gene_ids"
    assert "datasets" in p and "chains" in p
    assert isinstance(p["upstream_sources"], list)
    # at least HGNC + Ensembl + OMIM + NCBI Gene present for §1
    names = {u["name"] for u in p["upstream_sources"]}
    assert "HGNC" in names
    assert "Ensembl" in names
    assert "OMIM" in names

def test_section_provenance_upstream_dedupes_repeated_datasets():
    # §7 uses both uniprot and ensembl for go AND reactome; both map to the
    # same upstream-source brand entries — those should appear once each.
    sec = REGISTRY["7"]
    p = _section_provenance(sec)
    names = [u["name"] for u in p["upstream_sources"]]
    assert len(names) == len(set(names))


# ───── build_provenance — anchors + sections + meta ──────────────────

def test_build_provenance_top_level_keys():
    p = build_provenance(_bundle(), meta=META)
    for k in ("@context", "@type", "name", "isPartOf", "url",
              "anchors", "sections", "data_access", "generated_at",
              "atlas_version", "biobtree_version"):
        assert k in p

def test_build_provenance_is_schema_org_Dataset():
    p = build_provenance(_bundle(), meta=META)
    assert p["@context"] == "https://schema.org"
    assert p["@type"] == "Dataset"

def test_build_provenance_anchors_drawn_from_bundle():
    p = build_provenance(_bundle(), meta=META)
    a = p["anchors"]
    assert a["symbol"] == "TP53"
    assert a["hgnc_id"] == "HGNC:11998"
    assert a["ensembl_id"] == "ENSG00000141510"
    assert a["canonical_uniprot"] == "P04637"
    assert a["canonical_transcript"] == "ENST00000269305"

def test_build_provenance_anchors_tolerate_missing_sections():
    # ncRNA-style: no §3 / §6 entries
    b = {"1": {"symbol": "MALAT1", "hgnc_id": "HGNC:29665"}}
    p = build_provenance(b)
    assert p["anchors"]["symbol"] == "MALAT1"
    assert p["anchors"]["hgnc_id"] == "HGNC:29665"
    assert p["anchors"]["canonical_uniprot"] is None
    assert p["anchors"]["canonical_transcript"] is None

def test_build_provenance_meta_threads_through():
    p = build_provenance(_bundle(), meta=META)
    assert p["generated_at"] == META["generated_at"]
    assert p["atlas_version"] == META["atlas_version"]
    assert p["biobtree_version"] == META["biobtree_version"]

def test_build_provenance_no_meta_yields_null_meta_fields():
    p = build_provenance(_bundle())
    assert p["generated_at"] is None
    assert p["atlas_version"] is None
    assert p["biobtree_version"] is None

def test_build_provenance_sections_list_has_all_sections():
    p = build_provenance(_bundle(), meta=META)
    assert len(p["sections"]) == 13          # §1–§12 + §13 Human Protein Atlas
    ids = [s["id"] for s in p["sections"]]
    assert ids == [str(n) for n in range(1, 14)]

def test_build_provenance_each_section_has_required_fields():
    p = build_provenance(_bundle(), meta=META)
    for s in p["sections"]:
        for k in ("id", "name", "description", "datasets", "chains",
                  "upstream_sources", "produces"):
            assert k in s, f"section {s.get('id')} missing {k}"

def test_build_provenance_chains_are_strings():
    p = build_provenance(_bundle(), meta=META)
    for s in p["sections"]:
        for ch in s["chains"]:
            assert isinstance(ch, str)
            # canonical biobtree chain syntax begins with ">>"
            # (or contains it for compound; the filter-bracket case still has ">>")
            assert ">>" in ch

def test_build_provenance_urls_use_atlas_base():
    p = build_provenance(_bundle(), meta=META)
    assert p["isPartOf"] == "https://sugi.bio/atlas/gene/TP53/"
    assert p["url"] == "https://sugi.bio/atlas/gene/TP53/provenance.json"

def test_build_provenance_respects_custom_base_url():
    p = build_provenance(_bundle(), meta=META, base_url="https://x.example/atlas")
    assert p["isPartOf"] == "https://x.example/atlas/gene/TP53/"


# ───── particular upstream-source signal we promise ───────────────────

def test_clinvar_section_routes_to_ncbi_clinvar():
    """§6 variants should expose ClinVar as one of its upstreams."""
    p = build_provenance(_bundle(), meta=META)
    s6 = [s for s in p["sections"] if s["id"] == "6"][0]
    names = {u["name"] for u in s6["upstream_sources"]}
    assert "ClinVar" in names
    cv = [u for u in s6["upstream_sources"] if u["name"] == "ClinVar"][0]
    assert "ncbi.nlm.nih.gov/clinvar" in cv["url"]

def test_reactome_section_routes_to_reactome_org():
    p = build_provenance(_bundle(), meta=META)
    s7 = [s for s in p["sections"] if s["id"] == "7"][0]
    names = {u["name"] for u in s7["upstream_sources"]}
    assert "Reactome" in names

def test_disease_section_routes_to_gencc_mondo_hpo_gwas():
    p = build_provenance(_bundle(), meta=META)
    s12 = [s for s in p["sections"] if s["id"] == "12"][0]
    names = {u["name"] for u in s12["upstream_sources"]}
    for expected in ("GenCC", "Mondo Disease Ontology",
                     "Human Phenotype Ontology", "GWAS Catalog", "OMIM"):
        assert expected in names, f"missing upstream {expected!r} in §12"


# ───── serialization ──────────────────────────────────────────────────

def test_as_provenance_string_is_valid_json_with_trailing_newline():
    p = build_provenance(_bundle(), meta=META)
    s = as_provenance_string(p)
    assert s.endswith("\n")
    assert json.loads(s) == p
