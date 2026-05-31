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
def test_transport_parity(probe, monkeypatch):
    name, args = probe
    a = _load("urllib", monkeypatch)
    out_a = getattr(a, name)(*args)
    b = _load("urllib_pool", monkeypatch)
    out_b = getattr(b, name)(*args)
    # Pagination cursor is opaque + per-request token; ignore in equality.
    for blob in (out_a, out_b):
        blob.get("pagination", {}).pop("next_token", None)
    assert json.dumps(out_a, sort_keys=True) == json.dumps(out_b, sort_keys=True), (
        f"{name}{args} disagrees across transports")
