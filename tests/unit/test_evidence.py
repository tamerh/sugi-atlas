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


def test_component_percentile_and_rank_clause():
    P.reset()
    # frozen per-component distribution: 100 genes, drug_count mostly 0 with a
    # heavy right tail (the real shape — most genes aren't drug targets).
    P._DIST_COMP = {"gene": {"drug_count": sorted([0] * 90 + list(range(1, 11)))}}
    # a count of 10 is the max → top of the distribution
    assert P.component_percentile("gene", "drug_count", 10) == 100
    assert P.component_percentile("gene", "drug_count", 0) == 0
    # rank_clause: count clears the floor (10) AND lands in the top decile
    assert P.rank_clause("gene", "drug_count", 10) == " (top 1% of genes corpus-wide)"
    # below the absolute floor → no clause even if percentile is high
    assert P.rank_clause("gene", "drug_count", 3) == ""
    # an un-ranked component → no clause
    assert P.rank_clause("gene", "trial_count", 999) == ""
    # no distribution loaded → graceful empty (unit-test / no-batch path)
    P.reset()
    assert P.component_percentile("gene", "drug_count", 10) is None
    assert P.rank_clause("gene", "drug_count", 10) == ""


def test_components_for_reads_per_slug_map():
    P.reset()
    P._COMP = {"gene": {"FER": {"gwas_count": 33, "drug_count": 56}}}
    assert P.components_for("gene", "FER") == {"gwas_count": 33, "drug_count": 56}
    assert P.components_for("gene", "NOPE") == {}     # unknown slug
    assert P.components_for("drug", "X") == {}        # unknown type
    P.reset()
