"""Prominence scoring — component extraction, composite, percentile, lookup."""
from atlas.page import evidence as P


def test_components_have_stable_keys_and_int_values():
    # Every entity of a type gets the SAME component keys (0 when absent).
    full = {"6": {"clinvar_total": 100}, "10": {"molecule_count": 5,
            "civic_evidence_total": 3, "disease_trial_count": 2},
            "8": {"string_count": 900}, "12": {"gwas_total": 7}}
    empty = {}
    cf, ce = P.components("gene", full), P.components("gene", empty)
    assert set(cf) == set(ce)                       # identical key set
    assert all(isinstance(v, int) for v in cf.values())
    assert cf["variant_count"] == 100 and ce["variant_count"] == 0
    assert cf["civic_count"] == 3


def test_raw_signal_monotonic():
    lo = P.raw_signal("gene", P.components("gene", {"6": {"clinvar_total": 1}}))
    hi = P.raw_signal("gene", P.components("gene", {"6": {"clinvar_total": 9999},
                                                    "10": {"civic_evidence_total": 50}}))
    assert hi > lo > 0


def test_percentiles_span_0_to_100():
    raw = {"a": 0.0, "b": 1.0, "c": 5.0, "d": 50.0}
    pct = P.percentiles(raw)
    assert pct["a"] == 0 and pct["d"] == 100
    assert all(0 <= v <= 100 for v in pct.values())
    assert pct["a"] <= pct["b"] <= pct["c"] <= pct["d"]


def test_lookup_uses_loaded_scores_then_distribution():
    P.reset()
    P._PROM = {"gene": {"TP53": 99}}
    P._DIST = {"gene": [0.0, 1.0, 2.0, 3.0, 4.0]}          # 5-point frozen dist
    assert P.lookup("gene", "TP53", 0.0) == 99             # known slug wins
    # unknown slug → percentile against the frozen distribution
    assert P.lookup("gene", "NEW", 4.0) == 100
    assert P.lookup("gene", "NEW", 0.0) == 0
    assert P.lookup("drug", "X", 1.0) is None              # no dist → omit
    P.reset()
