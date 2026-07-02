"""Sugi Predict cross-link URL builders (atlas.predict). Gating against the
bundled 4,884-target manifest + SMILES query encoding."""
from atlas.predict import target_url, smiles_url, _covered


def test_manifest_loads():
    symbols, accs = _covered()
    assert len(symbols) > 4000 and len(accs) > 4000     # ~4,884 covered targets


def test_covered_target_links_by_symbol():
    # EGFR / P00533 is a known covered target.
    assert target_url(symbol="EGFR") == "https://sugi.bio/predict/target/EGFR/"


def test_covered_target_case_insensitive_and_acc_fallback():
    assert target_url(symbol="egfr") == "https://sugi.bio/predict/target/egfr/"
    # symbol unknown but accession covered → still links (via the acc), URL by acc
    assert target_url(symbol="", uniprot="P00533") == "https://sugi.bio/predict/target/P00533/"


def test_non_target_returns_none():
    # MALAT1 (lncRNA) is not a Predict target → must NOT link (would 404).
    assert target_url(symbol="MALAT1") is None
    assert target_url(symbol="NOT_A_GENE", uniprot="X99999") is None
    assert target_url(symbol="", uniprot="") is None


def test_smiles_url_encodes_special_chars():
    # imatinib SMILES — '=', '(', ')', '#' etc. must be percent-encoded.
    url = smiles_url("CC(=O)OC1=CC=CC=C1C(=O)O")  # aspirin
    assert url.startswith("https://sugi.bio/predict/?q=")
    assert "=" not in url.split("?q=", 1)[1]        # raw '=' encoded, not left bare
    assert "%3D" in url                              # '=' → %3D


def test_smiles_url_empty():
    assert smiles_url("") is None
    assert smiles_url(None) is None
