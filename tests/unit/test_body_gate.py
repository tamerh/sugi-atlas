"""Unit tests for atlas.validation.body_gate.

Covers the pure logic that drives the publish-gating verdict:
  - _size()        — leaf-size scoring
  - _flat_sizes()  — bundle flattening (skips meta keys)
  - diff()         — bundle-vs-snapshot structural diff
  - verdict()      — clean / drift / regression / first_run classification

No biobtree, no filesystem (snapshot file I/O is exercised separately in the
integration entry point).
"""
import pytest

from atlas.validation import body_gate


# ───── _size ──────────────────────────────────────────────────────────

@pytest.mark.parametrize("v,expected", [
    (None,        0),
    ("",          0),
    ([],          0),
    ({},          0),
    ([1, 2, 3],   3),
    ({"a": 1},    1),
    ("x",         1),
    (42,          1),
    ("nonempty",  1),
])
def test_size(v, expected):
    assert body_gate._size(v) == expected


# ───── _flat_sizes ────────────────────────────────────────────────────

def test_flat_sizes_skips_meta_keys():
    bundle = {
        "1": {"section": "01_gene_ids", "symbol": "TP53",
              "hgnc_id": "HGNC:11998", "mim": ["191170"]},
    }
    flat = body_gate._flat_sizes(bundle)
    # section/symbol must NOT appear
    assert ("1", "section") not in flat
    assert ("1", "symbol") not in flat
    # real keys present with correct sizes
    assert flat[("1", "hgnc_id")] == 1
    assert flat[("1", "mim")] == 1   # list of length 1


def test_flat_sizes_handles_non_dict_section():
    # defensive: a non-dict section value should not blow up
    bundle = {"1": "garbage", "2": {"x": [1, 2]}}
    flat = body_gate._flat_sizes(bundle)
    assert flat == {("2", "x"): 2}


# ───── diff ───────────────────────────────────────────────────────────

def test_diff_first_run_when_no_snapshot():
    new = {"1": {"x": [1, 2]}}
    d = body_gate.diff(None, new)
    assert d == [{"kind": "first_run"}]


def test_diff_clean_when_identical():
    b = {"1": {"x": [1, 2, 3]}}
    assert body_gate.diff(b, b) == []
    assert body_gate.diff(b, {"1": {"x": [1, 2, 3]}}) == []


def test_diff_detects_added_key():
    old = {"1": {"a": [1]}}
    new = {"1": {"a": [1], "b": [2, 3]}}
    [entry] = body_gate.diff(old, new)
    assert entry == {"section": "1", "key": "b", "old": 0, "new": 2, "kind": "added"}


def test_diff_detects_removed_key():
    old = {"1": {"a": [1], "b": [2, 3]}}
    new = {"1": {"a": [1]}}
    [entry] = body_gate.diff(old, new)
    assert entry == {"section": "1", "key": "b", "old": 2, "new": 0, "kind": "removed"}


def test_diff_detects_count_changed():
    old = {"1": {"a": [1, 2, 3]}}
    new = {"1": {"a": [1, 2, 3, 4, 5]}}
    [entry] = body_gate.diff(old, new)
    assert entry["kind"] == "changed"
    assert entry["old"] == 3 and entry["new"] == 5


def test_diff_handles_multiple_sections_sorted():
    old = {"1": {"a": [1]}, "2": {"b": [1, 2]}}
    new = {"1": {"a": [1, 2]}, "2": {"b": [1]}}
    out = body_gate.diff(old, new)
    # results are sorted by (section, key) so ordering is stable
    assert [(d["section"], d["key"]) for d in out] == [("1", "a"), ("2", "b")]


# ───── verdict ────────────────────────────────────────────────────────

def test_verdict_clean_for_empty_diff():
    assert body_gate.verdict([]) == "clean"


def test_verdict_first_run_short_circuits():
    diffs = [{"kind": "first_run"}]
    assert body_gate.verdict(diffs) == "first_run"


def test_verdict_first_run_wins_over_other_diffs():
    # implementation detail: if first_run marker is present, verdict is first_run
    # regardless of any other diff entries
    diffs = [{"kind": "first_run"},
             {"section": "1", "key": "x", "old": 10, "new": 1, "kind": "changed"}]
    assert body_gate.verdict(diffs) == "first_run"


def test_verdict_regression_on_any_removed():
    diffs = [{"section": "1", "key": "x", "old": 5, "new": 0, "kind": "removed"}]
    assert body_gate.verdict(diffs) == "regression"


def test_verdict_regression_when_count_halves():
    # > 50% drop → regression
    diffs = [{"section": "1", "key": "x", "old": 100, "new": 49, "kind": "changed"}]
    assert body_gate.verdict(diffs) == "regression"


def test_verdict_drift_for_minor_count_change():
    # within 50% → drift
    diffs = [{"section": "1", "key": "x", "old": 100, "new": 110, "kind": "changed"}]
    assert body_gate.verdict(diffs) == "drift"


def test_verdict_drift_for_small_drop():
    # 100 -> 75 = 25% drop, within tolerance → drift, not regression
    diffs = [{"section": "1", "key": "x", "old": 100, "new": 75, "kind": "changed"}]
    assert body_gate.verdict(diffs) == "drift"


def test_verdict_drift_on_added_only():
    diffs = [{"section": "1", "key": "x", "old": 0, "new": 10, "kind": "added"}]
    assert body_gate.verdict(diffs) == "drift"


def test_verdict_regression_count_changed_to_zero_is_treated_as_removed_threshold():
    # old=10, new=0 with kind=removed -> regression
    diffs = [{"section": "1", "key": "x", "old": 10, "new": 0, "kind": "removed"}]
    assert body_gate.verdict(diffs) == "regression"


# ───── check() — composed behavior, no I/O via monkeypatching ─────────

def test_check_composes_diff_verdict_and_summary(monkeypatch):
    """check() should: load snapshot (we stub None), compute diff, compute
    verdict, and report a human-readable summary string."""
    bundle = {"1": {"a": [1]}}
    # stub the loader to return None (first run scenario)
    monkeypatch.setattr(body_gate, "load_snapshot", lambda snap_dir, sym: None)
    r = body_gate.check("FAKE", bundle, snap_dir="/nowhere")
    assert r["verdict"] == "first_run"
    assert "no snapshot" in r["summary"].lower()


def test_check_clean_against_identical(monkeypatch):
    bundle = {"1": {"a": [1, 2]}}
    monkeypatch.setattr(body_gate, "load_snapshot", lambda snap_dir, sym: bundle)
    r = body_gate.check("FAKE", bundle, snap_dir="/nowhere")
    assert r["verdict"] == "clean"
    assert "byte-identical" in r["summary"]


def test_check_summary_counts_changes(monkeypatch):
    """Summary string should report n_changed / n_added / n_removed."""
    old = {"1": {"a": [1], "b": [1, 2]}}
    new = {"1": {"a": [1, 2, 3], "c": [9]}}  # 'b' removed, 'c' added, 'a' changed
    monkeypatch.setattr(body_gate, "load_snapshot", lambda snap_dir, sym: old)
    r = body_gate.check("FAKE", new, snap_dir="/nowhere")
    # 'b' removed -> regression verdict
    assert r["verdict"] == "regression"
    assert "1 changed" in r["summary"]
    assert "1 added" in r["summary"]
    assert "1 removed" in r["summary"]
