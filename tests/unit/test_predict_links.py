"""Sugi Predict cross-link URL builders (atlas.predict). Gating against the
bundled 4,884-target manifest + SMILES query encoding."""
import pytest

from atlas import predict
from atlas.predict import target_url, compound_url, resolve_schembl, PredictResolverError, _covered


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


def test_compound_url():
    assert compound_url("SCHEMBL10883") == "https://sugi.bio/predict/compound/SCHEMBL10883"
    assert compound_url("") is None
    assert compound_url(None) is None


def test_resolve_schembl_empty_is_none_without_calling():
    assert resolve_schembl("") is None
    assert resolve_schembl(None) is None


def test_resolve_schembl_raises_when_resolver_down(monkeypatch):
    # A down resolver must FAIL LOUDLY (no silent fallback). Point it at a closed
    # port so the connection is refused immediately.
    monkeypatch.setattr(predict, "_RESOLVER", "http://127.0.0.1:1/predict")
    with pytest.raises(PredictResolverError):
        resolve_schembl("CC(=O)OC1=CC=CC=C1C(=O)O")
