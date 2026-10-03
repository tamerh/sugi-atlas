"""resolve_mondo: when biobtree has no Mondo name for a MONDO:id (deprecated/odd
term), recover the disease name from a GenCC curation (disease_title) rather than
leaving the page titled with the raw id (audit: 7 content-rich pages, e.g.
MONDO:0011271 → 'rigid spine muscular dystrophy 1')."""
from atlas.disease import anchors


def test_mondo_name_recovered_from_gencc(monkeypatch):
    # biobtree returns an empty Mondo entry (no name) and the search index has no
    # matching row — the nameless-id case.
    monkeypatch.setattr(anchors, "entry", lambda *a, **k: {"Attributes": {"Mondo": {}}})
    monkeypatch.setattr(anchors, "search", lambda *a, **k: {})
    monkeypatch.setattr(anchors, "rows", lambda resp: [])
    monkeypatch.setattr(anchors, "map_all",
                        lambda ident, chain, **k: [{"disease_title": "rigid spine muscular dystrophy 1"}]
                        if chain == ">>mondo>>gencc" else [])
    mid, _en, canonical = anchors.resolve_mondo("MONDO:0011271")
    assert mid == "MONDO:0011271"
    assert canonical == "rigid spine muscular dystrophy 1"   # not the raw id


def test_mondo_no_recovery_leaves_none(monkeypatch):
    # No Mondo name and no GenCC curation → canonical stays None (caller falls back
    # to the id, as before — we only improved the recoverable case).
    monkeypatch.setattr(anchors, "entry", lambda *a, **k: {"Attributes": {"Mondo": {}}})
    monkeypatch.setattr(anchors, "search", lambda *a, **k: {})
    monkeypatch.setattr(anchors, "rows", lambda resp: [])
    monkeypatch.setattr(anchors, "map_all", lambda *a, **k: [])
    _mid, _en, canonical = anchors.resolve_mondo("MONDO:9999999")
    assert canonical is None
