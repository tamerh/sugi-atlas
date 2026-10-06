"""Molecular data pickups #3 (GtoPdb) / #5 (BindingDB) / #2 (BRENDA) — pure
helpers (collectors hit biobtree and are exercised live, not here)."""
from atlas.gene.sections.s10_drugs import _affinity_nm, _clean_ligand


# ── #5 BindingDB affinity normalization ──────────────────────────────────────
def test_affinity_nm_units():
    assert _affinity_nm("1.0 nM") == 1.0
    assert _affinity_nm("0.5 µM") == 500.0
    assert _affinity_nm("0.5 uM") == 500.0
    assert _affinity_nm("2 mM") == 2_000_000.0
    assert _affinity_nm("100 pM") == 0.1
    assert _affinity_nm(">5 nM") == 5.0          # qualifier stripped


def test_affinity_nm_rejects_junk_and_zero():
    assert _affinity_nm("") is None
    assert _affinity_nm(None) is None
    assert _affinity_nm("n/a") is None
    assert _affinity_nm("0 nM") is None          # zero is not a real affinity
    assert _affinity_nm("5 furlongs") is None    # unknown unit


# ── #5 BindingDB ligand-name cleanup ──────────────────────────────────────────
def test_clean_ligand_picks_readable_name():
    assert _clean_ligand("Imatinib::CHEMBL941::STI-571") == "Imatinib"
    # all-CHEMBL → falls back to first segment
    assert _clean_ligand("CHEMBL941::CHEMBL2") == "CHEMBL941"
    assert _clean_ligand("") == ""
    assert _clean_ligand(None) is None


def test_clean_ligand_returns_full_name():
    # No baked-in truncation: a long IUPAC name comes back whole (the frontend
    # clamps for display). Web-team report: names were arriving cut with "…".
    long = "x" * 100
    out = _clean_ligand(long)
    assert out == long and "…" not in out


# ── §5 ortholog organism recovery (biobtree #37) ─────────────────────────────
def test_ortholog_organism_from_id():
    from atlas.gene.sections.s05_orthologs import _organism_from_id as f
    # biobtree leaves WormBase ortholog genome empty → recover from the id prefix
    assert f("WBGENE00000264") == "caenorhabditis_elegans"
    assert f("WBGene00018327") == "caenorhabditis_elegans"   # case-insensitive
    assert f("FBgn0000546") == "drosophila_melanogaster"
    # Ensembl-namespaced orthologs are populated upstream → no fallback
    assert f("ENSMUSG00000017146") == ""
    assert f("") == "" and f(None) == ""


def test_pharmgkb_summary_annotation_token(monkeypatch):
    # PharmGKB renamed clinical annotations → "SummaryAnnotation"; the old token-only
    # filter rendered clinical annotations on 0 drug pages.
    from atlas.drug.sections import s09_pharmacogenomics as S
    monkeypatch.setattr(S, "entry", lambda *a, **k: {"Attributes": {"Pharmgkb": {"related_genes": [
        {"gene_symbol": "VKORC1", "evidence_type": "SummaryAnnotation,VariantAnnotation"},
        {"gene_symbol": "XYZ", "evidence_type": "VariantAnnotation"}]}}})
    monkeypatch.setattr(S, "map_all", lambda sym, chain, **k: [
        {"variant": "rs9923231", "chemicals": "warfarin", "type": "Dosage",
         "level_of_evidence": "1A", "phenotypes": ""}] if sym == "VKORC1" else [])
    rows = S._clinical_annotations("PA451906", "warfarin")
    assert [(r["gene"], r["level"]) for r in rows] == [("VKORC1", "1A")]


def test_curated_target_genes_order():
    # ChEMBL curated MOA genes lead; GtoPdb primary targets follow, strongest first.
    from atlas.page.drug_declarative import _curated_target_genes
    b2 = {"mechanism_genes": [{"gene_symbol": "ABL1"}, {"gene_symbol": "KIT"}],
          "primary_targets": [{"gene_symbol": "TRHR", "affinity": "5.1"},
                              {"gene_symbol": "ABL1", "affinity": "8.0"},
                              {"gene_symbol": "DDR1", "affinity": "8.5"}]}
    assert _curated_target_genes(b2) == ["ABL1", "KIT", "DDR1", "TRHR"]
    assert _curated_target_genes({"primary_targets": [{"gene_symbol": "TRHR", "affinity": "5.15"},
                                                      {"gene_symbol": "GABRA1", "affinity": "7.79"}]}) == ["GABRA1", "TRHR"]


def test_drug_mechanisms_salt_children_genes_and_organisms(monkeypatch):
    from atlas.drug.sections import s02_targets as S
    def fake(root, chain, **k):
        if chain == ">>chembl_molecule>>chembl_mechanism":
            return [{"mechanism_of_action": "X inhibitor", "action_type": "INHIBITOR",
                     "target_name": "X", "target_type": "SINGLE PROTEIN"}] if root == "SALT" else []
        if chain == ">>chembl_molecule>>chembl_mechanism>>chembl_target":
            return [{"id": "T1", "type": "SINGLE PROTEIN"}, {"id": "T2", "type": "SINGLE PROTEIN"},
                    {"id": "T3", "type": "PROTEIN COMPLEX GROUP"}] if root == "SALT" else []
        if chain == ">>chembl_target>>uniprot>>hgnc":
            return [{"id": "HGNC:76"}] if root == "T1" else []
        if chain == ">>chembl_target>>taxonomy":
            return [{"id": "1773", "name": "Mycobacterium tuberculosis"}] if root == "T2" else []
        return []
    monkeypatch.setattr(S, "map_all", fake)
    monkeypatch.setattr(S, "_hgnc_symbol", lambda h: {"HGNC:76": "ABL1"}.get(h))
    moa, genes, orgs = S._mechanisms("PARENT", ("SALT",))
    assert len(moa) == 1                                   # found via the salt form
    assert [g["gene_symbol"] for g in genes] == ["ABL1"]   # single protein → gene; complex not expanded
    assert orgs == ["Mycobacterium tuberculosis"]


def test_drug_mechanisms_complex_target_organism(monkeypatch):
    # telaprevir: HCV NS3/NS4A is a PROTEIN COMPLEX — organism still reported; a human
    # complex (Homo sapiens) is not, and complexes are never expanded into genes.
    from atlas.drug.sections import s02_targets as S
    def fake(root, chain, **k):
        if chain == ">>chembl_molecule>>chembl_mechanism>>chembl_target":
            return [{"id": "CPX_HCV", "type": "PROTEIN COMPLEX"},
                    {"id": "CPX_HUMAN", "type": "PROTEIN COMPLEX GROUP"}]
        if chain == ">>chembl_target>>taxonomy":
            return [{"name": {"CPX_HCV": "Orthohepacivirus hominis",
                              "CPX_HUMAN": "Homo sapiens"}[root]}]
        if chain == ">>chembl_target>>uniprot>>hgnc":
            raise AssertionError("complexes must not be expanded into genes")
        return []
    monkeypatch.setattr(S, "map_all", fake)
    moa, genes, orgs = S._mechanisms("CHEMBL231813")
    assert genes == [] and orgs == ["Orthohepacivirus hominis"]
