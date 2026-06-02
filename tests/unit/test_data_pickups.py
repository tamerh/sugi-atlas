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


def test_clean_ligand_truncates():
    long = "x" * 100
    out = _clean_ligand(long)
    assert out.endswith("…") and len(out) == 61
