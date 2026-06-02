"""Phase B audit-polish fixes (#6–#13) — lock in the new helpers' behavior."""
from atlas.render_common import is_ontology_id, display_name, gencc_rank
from atlas.drug.roles import pharma_class
from atlas.page.jsonld_inline import compact_for_inline, INLINE_CAP
from atlas.page import links


# ── #6 inline JSON-LD cap ──────────────────────────────────────────────────
def test_compact_for_inline_caps_long_arrays():
    j = {"name": "X", "associatedGene": [{"name": f"G{i}"} for i in range(50)],
         "sameAs": ["a", "b"]}
    c = compact_for_inline(j)
    assert len(c["associatedGene"]) == INLINE_CAP
    assert c["sameAs"] == ["a", "b"]          # short arrays untouched
    assert "comment" in c                      # truncation noted
    assert j["associatedGene"][0]["name"] == "G0"  # original not mutated


def test_compact_for_inline_noop_when_short():
    j = {"name": "X", "associatedGene": [{"name": "G0"}]}
    c = compact_for_inline(j)
    assert "comment" not in c
    assert c == j


def test_compact_for_inline_caps_reverse_edges():
    j = {"name": "X", "@reverse": {"associatedGene": list(range(40))}}
    c = compact_for_inline(j)
    assert len(c["@reverse"]["associatedGene"]) == INLINE_CAP


# ── #9 pharmacological-role selection ──────────────────────────────────────
def test_pharma_class_skips_environmental_and_toxicology():
    assert pharma_class(["environmental contaminant", "xenobiotic"],
                        fallback="small molecule") == "small molecule"
    # caffeine: leads with 'mutagen', real class deeper
    assert pharma_class(["mutagen", "central nervous system stimulant",
                         "adenosine receptor antagonist"]) == "central nervous system stimulant"


def test_pharma_class_keeps_real_class():
    assert pharma_class(["tyrosine kinase inhibitor"]) == "tyrosine kinase inhibitor"
    assert pharma_class([], fallback="drug") == "drug"
    assert pharma_class(["mutagen"]) is None       # only junk → fallback (None)


# ── #11 raw-ontology-id detection ──────────────────────────────────────────
def test_is_ontology_id():
    assert is_ontology_id("MONDO:0004992")
    assert is_ontology_id("EFO:0010282")
    assert is_ontology_id("MP:0001914")
    assert not is_ontology_id("breast carcinoma")
    assert not is_ontology_id("")
    assert not is_ontology_id(None)


# ── #12 de-SHOUT display name ──────────────────────────────────────────────
def test_display_name_de_shouts_but_spares_symbols():
    assert display_name("IMATINIB") == "Imatinib"
    assert display_name("WATER") == "Water"
    assert display_name("acitretin") == "acitretin"      # mixed/lower untouched
    # gene symbols are upper but must NOT be routed through here — verified by
    # build_meta exempting genes; the helper itself would title-case them.
    assert display_name("Imatinib") == "Imatinib"


# ── #13 GenCC rank + destination-canonical relabel ─────────────────────────
def test_gencc_rank_order():
    assert gencc_rank("Definitive") > gencc_rank("Strong") > gencc_rank("Limited")
    assert gencc_rank("Refuted") == 0
    assert gencc_rank(None) == 0


def test_drug_mesh_label_de_shouts_and_salt_strips():
    from atlas.page.links import _drug_display
    assert _drug_display("CISPLATIN") == "Cisplatin"
    assert _drug_display("IMATINIB MESYLATE") == "Imatinib"       # de-SHOUT + salt-strip
    assert _drug_display("CANDESARTAN CILEXETIL") == "Candesartan Cilexetil"


def test_gencc_dedup_prefers_on_disease_record():
    from atlas.disease.render import _dedup_gencc
    rows = [
        # off-disease but stronger classification — must NOT outrank the on-disease one
        {"symbol": "SUFU", "gencc_classification": "Definitive", "mondo_disease": "Joubert syndrome"},
        {"symbol": "SUFU", "gencc_classification": "Moderate", "mondo_disease": "medulloblastoma"},
        # gene with only an off-disease record — kept, sorted after on-disease genes
        {"symbol": "BRCA2", "gencc_classification": "Definitive",
         "mondo_disease": "Fanconi anemia complementation group D1"},
    ]
    ded = _dedup_gencc(rows, "medulloblastoma")
    by = {best["symbol"]: (best, on) for best, _n, on in ded}
    assert by["SUFU"][0]["mondo_disease"] == "medulloblastoma"     # on-disease preferred
    assert by["SUFU"][1] is True
    assert by["BRCA2"][1] is False
    assert [best["symbol"] for best, _n, _on in ded][0] == "SUFU"  # on-disease sorted first


def test_reverse_mesh_surfaces_incoming_biomarker_genes():
    """A gene→drug (CIViC) edge becomes a 'Biomarker genes' link on the drug."""
    links.reset()
    links._REVERSE = {"/atlas/drug/venetoclax/": [
        ["TP53", "/atlas/gene/TP53/", "gene", "Drugs"],
        ["BCL2", "/atlas/gene/BCL2/", "gene", "Drugs"],
    ]}
    blk = links.related_block("drug", {"2": {}, "4": {}, "7": {}, "10": {}},
                              slug="venetoclax")
    assert "**Biomarker genes:**" in blk
    assert "[TP53](/atlas/gene/TP53/)" in blk
    links.reset()


def test_reverse_mesh_dedups_against_forward_edges():
    """A gene already listed as a forward target is not repeated under reverse."""
    links.reset()
    links._MANIFEST = {"gene": {"BCL2": "bcl2"}, "disease": {}, "drug": {}}  # value = slug
    links._REVERSE = {"/atlas/drug/venetoclax/": [
        ["BCL2", "/atlas/gene/bcl2/", "gene", "Drugs"],   # same url as the forward target
        ["TP53", "/atlas/gene/TP53/", "gene", "Drugs"],
    ]}
    bundle = {"2": {"primary_targets": [{"gene_symbol": "BCL2", "source": "gtopdb"}]},
              "4": {}, "7": {}, "10": {}}
    blk = links.related_block("drug", bundle, slug="venetoclax")
    bio = [l for l in blk.splitlines() if l.startswith("- **Biomarker genes:**")]
    assert bio and "TP53" in bio[0] and "bcl2" not in bio[0].lower()  # BCL2 deduped
    links.reset()


def test_drug_gene_mesh_is_gtopdb_curated_only():
    """#3 extension: only GtoPdb-curated drug targets become gene links — ChEMBL
    bioactivity assay hits (Salmeterol 'targets' KDM4E/TP53) are dropped."""
    links.reset()
    links._MANIFEST = {"gene": {"ABL1": "abl1", "KDM4E": "kdm4e"}, "disease": {}, "drug": {}}
    bundle = {"2": {"primary_targets": [
        {"gene_symbol": "ABL1", "source": "gtopdb"},
        {"gene_symbol": "KDM4E", "source": "chembl"}]}, "4": {}, "7": {}, "10": {}}
    genes = [g for g, _u in links.related_targets("drug", bundle)["Genes"]]
    assert genes == ["ABL1"]                       # chembl bioactivity dropped
    links.reset()


def test_disease_drug_mesh_excludes_bioactivity():
    """#3 extension: disease→drug = title-validated trials + CIViC; NOT
    b10.drugs (bioactivity hits on cohort targets → off-target junk)."""
    links.reset()
    links._MANIFEST = {"gene": {}, "drug": {"Cisplatin": "cisplatin",
                                            "Clotrimazole": "clotrimazole"}, "disease": {}}
    bundle = {"4": {}, "5": {}, "13": {"trial_drugs": [{"name": "Cisplatin"}]},
              "10": {"drugs": [{"name": "Clotrimazole"}]}}   # bioactivity junk
    drugs = [d for d, _u in links.related_targets("disease", bundle)["Drugs"]]
    assert "Cisplatin" in drugs and "Clotrimazole" not in drugs
    links.reset()


def test_reverse_mesh_absent_without_slug_or_index():
    links.reset()
    blk = links.related_block("drug", {"2": {}, "4": {}, "7": {}, "10": {}})  # no slug
    assert "Biomarker genes" not in blk
    links.reset()


def test_related_block_always_emits_with_placeholder():
    """PAGE_CONTRACT: #related is always present (placeholder when no links)."""
    links.reset()
    blk = links.related_block("gene", {}, slug="NOPE")
    assert "## Related Atlas pages {#related}" in blk
    assert "No linked Atlas pages yet" in blk
    links.reset()


def test_yaml_escape_strips_control_chars():
    from atlas.pipeline import _yaml_escape
    assert _yaml_escape("ola\x7fparib") == "olaparib"      # 0x7f stripped
    assert _yaml_escape('a "b"') == 'a \\"b\\"'            # quotes still escaped
    assert _yaml_escape("keep\ttab") == "keep\ttab"        # tab kept


def test_canonical_label_resolves_destination(monkeypatch):
    links.reset()
    links._MANIFEST = {"gene": {}, "drug": {},
                       "disease": {"schizoaffective disorder": "schizophrenia",
                                   "schizophrenia": "schizophrenia"}}
    links._CANON = {"gene": {}, "drug": {},
                    "disease": {"schizophrenia": "schizophrenia"}}
    assert links.canonical_label("/atlas/disease/schizophrenia/") == "schizophrenia"
    assert links.canonical_label("/atlas/disease/unknown/") is None
    assert links.canonical_label("not-a-url") is None
    # a synonym reference relabels to the page it actually opens
    bundle = {"2": {}, "4": {}, "7": {},
              "10": {"civic_evidence": [{"disease": "schizoaffective disorder",
                                         "therapy": "metformin"}]}}
    diseases = links.related_targets("drug", bundle)["Diseases"]
    assert diseases == [("schizophrenia", "/atlas/disease/schizophrenia/")]
    links.reset()
