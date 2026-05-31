"""Cross-transport parity: urllib and urllib_pool must return identical
search/entry/bbmap dicts. If a future transport (e.g. gRPC) silently
diverges in field-naming or coercion, this test breaks before it reaches
a section collector.

Requires biobtree REST running on 127.0.0.1:8000 — skips otherwise.
"""
import importlib
import json
import socket
import sys

import pytest


def _biobtree_up() -> bool:
    try:
        with socket.create_connection(("127.0.0.1", 8000), timeout=0.5):
            return True
    except OSError:
        return False


pytestmark = pytest.mark.skipif(
    not _biobtree_up(),
    reason="biobtree REST not reachable on 127.0.0.1:8000",
)


def _fresh(module_path: str):
    # Drop the cached transport modules so each import resolves freshly
    # for the current value of ATLAS_BIOBTREE_TRANSPORT.
    for k in list(sys.modules):
        if k.startswith("atlas.biobtree"):
            del sys.modules[k]
    return importlib.import_module(module_path)


def _load(transport: str, monkeypatch):
    monkeypatch.setenv("ATLAS_BIOBTREE_TRANSPORT", transport)
    return _fresh("atlas.biobtree")


@pytest.mark.parametrize("probe", [
    ("search", ("TP53", "hgnc")),
    ("entry",  ("HGNC:11998", "hgnc")),
    ("bbmap",  ("HGNC:11998", ">>hgnc>>uniprot")),
])
def test_urllib_urllib_pool_parity(probe, monkeypatch):
    """urllib + urllib_pool must agree byte-for-byte (same wire protocol)."""
    name, args = probe
    a = _load("urllib", monkeypatch)
    out_a = getattr(a, name)(*args)
    b = _load("urllib_pool", monkeypatch)
    out_b = getattr(b, name)(*args)
    for blob in (out_a, out_b):
        blob.get("pagination", {}).pop("next_token", None)
    assert json.dumps(out_a, sort_keys=True) == json.dumps(out_b, sort_keys=True), (
        f"{name}{args} disagrees across urllib transports")


def _grpc_up() -> bool:
    try:
        with socket.create_connection(("127.0.0.1", 7776), timeout=0.5):
            return True
    except OSError:
        return False


@pytest.mark.skipif(not _grpc_up(), reason="biobtree gRPC not reachable on :7776")
def test_grpc_search_parity(monkeypatch):
    """gRPC produces REST-shape search output via the adapter.

    Search schema is fixed (`id|dataset|name|xref_count`), so byte parity
    is achievable. Data rows must be a set-equal match."""
    a = _load("urllib", monkeypatch)
    out_a = a.search("TP53", "hgnc")
    b = _load("grpc", monkeypatch)
    out_b = b.search("TP53", "hgnc")
    assert out_a["schema"] == out_b["schema"]
    assert sorted(out_a["data"]) == sorted(out_b["data"])


@pytest.mark.skipif(not _grpc_up(), reason="biobtree gRPC not reachable on :7776")
def test_grpc_entry_parity(monkeypatch):
    """Entry: Attributes PascalCase keys + xrefs schema|data must align."""
    a = _load("urllib", monkeypatch)
    out_a = a.entry("HGNC:11998", "hgnc")
    b = _load("grpc", monkeypatch)
    out_b = b.entry("HGNC:11998", "hgnc")
    assert set(out_a.get("Attributes", {}).keys()) == set(out_b.get("Attributes", {}).keys())
    # xrefs are a set of "dataset|count" rows; order doesn't matter
    assert sorted(out_a["xrefs"]["data"]) == sorted(out_b["xrefs"]["data"])


@pytest.mark.skipif(not _grpc_up(), reason="biobtree gRPC not reachable on :7776")
def test_grpc_map_parity(monkeypatch):
    """Mapping: targets (pipe-encoded) set-equality + matching schema."""
    a = _load("urllib", monkeypatch)
    out_a = a.bbmap("HGNC:11998", ">>hgnc>>uniprot")
    b = _load("grpc", monkeypatch)
    out_b = b.bbmap("HGNC:11998", ">>hgnc>>uniprot")
    assert out_a["schema"] == out_b["schema"]
    rest_targets = sorted(t for m in (out_a.get("mappings") or []) for t in (m.get("targets") or []))
    grpc_targets = sorted(t for m in (out_b.get("mappings") or []) for t in (m.get("targets") or []))
    assert rest_targets == grpc_targets
