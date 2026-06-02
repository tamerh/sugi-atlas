"""Direct-Go vs gate transport routing in atlas.biobtree.client (hermetic —
_POOL.request is stubbed; no network)."""
import json
import pytest

import atlas.biobtree.client as C


class _Resp:
    def __init__(self, body, status=200):
        self.status = status
        self.data = json.dumps(body).encode()


def _capture(monkeypatch, body):
    cap = {}

    def fake(method, url, fields=None):
        cap["url"], cap["fields"] = url, dict(fields or {})
        return _Resp(body)
    monkeypatch.setattr(C._POOL, "request", fake)
    return cap


def test_direct_map_uses_ws_and_injects_mode_lite(monkeypatch):
    monkeypatch.setattr(C, "TRANSPORT", "direct")
    cap = _capture(monkeypatch, {"mappings": [], "schema": ""})
    C.bbmap("TP53", ">>hgnc>>uniprot")
    assert cap["url"].endswith("/ws/map/")
    assert cap["fields"]["mode"] == "lite"
    assert cap["fields"]["m"] == ">>hgnc>>uniprot"
    assert cap["fields"]["i"] == "TP53"


def test_direct_search_uses_ws_root_with_mode(monkeypatch):
    monkeypatch.setattr(C, "TRANSPORT", "direct")
    cap = _capture(monkeypatch, {"data": [], "schema": ""})
    C.search("TP53")
    assert cap["url"].endswith("/ws/")
    assert cap["fields"]["mode"] == "lite"


def test_direct_entry_takes_no_mode(monkeypatch):
    monkeypatch.setattr(C, "TRANSPORT", "direct")
    cap = _capture(monkeypatch, {"xrefs": {}})
    C.entry("HGNC:11998", "hgnc")
    assert cap["url"].endswith("/ws/entry/")
    assert "mode" not in cap["fields"]            # Go hardwires entryLite


def test_gate_uses_api_paths_without_mode(monkeypatch):
    monkeypatch.setattr(C, "TRANSPORT", "gate")
    cap = _capture(monkeypatch, {"mappings": []})
    C.bbmap("TP53", ">>x")
    assert cap["url"].endswith("/api/map")
    assert "mode" not in cap["fields"]


def test_inline_err_body_raises(monkeypatch):
    """Direct Go returns query errors inline with HTTP 200 — must still raise."""
    monkeypatch.setattr(C, "TRANSPORT", "direct")
    _capture(monkeypatch, {"Err": "unknown dataset: 'notadataset'"})
    with pytest.raises(C.BiobtreeError):
        C.bbmap("TP53", ">>hgnc>>notadataset")
