"""ORA statistics — hypergeometric tail + BH-FDR + fold (atlas.ora)."""
import math
from atlas.ora import hypergeom_sf, bh_fdr, enrich


def test_hypergeom_sf_exact_small():
    # N=4 population, K=2 successes, n=2 draws: P(X=0)=1/6, P(X=1)=4/6, P(X=2)=1/6
    assert abs(hypergeom_sf(0, 4, 2, 2) - 1.0) < 1e-12
    assert abs(hypergeom_sf(1, 4, 2, 2) - 5 / 6) < 1e-12      # P(>=1)
    assert abs(hypergeom_sf(2, 4, 2, 2) - 1 / 6) < 1e-12      # P(>=2)
    assert hypergeom_sf(3, 4, 2, 2) == 0.0                    # impossible


def test_hypergeom_sf_stable_large():
    # Large universe shouldn't overflow; a strongly over-represented hit is tiny p
    p = hypergeom_sf(9, 11000, 12, 68)        # MITF-M-like: 9 of 68 in a 12-gene pathway
    assert 0.0 < p < 1e-9


def test_bh_fdr_monotone():
    adj = bh_fdr([0.01, 0.02, 0.03, 0.04])
    assert all(abs(a - 0.04) < 1e-12 for a in adj)            # all p*4/rank == 0.04
    # monotone non-decreasing in p-order
    assert bh_fdr([0.0001, 0.5]) == sorted(bh_fdr([0.0001, 0.5]))
    assert bh_fdr([]) == []


def test_enrich_keeps_counts_adds_stats():
    items = [{"id": "big", "k": 17, "K": 1179},     # umbrella: many genes, low fold
             {"id": "mitf", "k": 9, "K": 12},        # specific: high fold, tiny fdr
             {"id": "none", "k": 1, "K": 0}]         # untestable → passthrough
    out = enrich(items, cohort_n=68, universe_n=11000)
    by = {o["id"]: o for o in out}
    assert by["mitf"]["fold"] > by["big"]["fold"]            # specific beats umbrella
    assert by["mitf"]["fdr"] < by["big"]["fdr"]
    assert by["none"]["fold"] is None and by["none"]["fdr"] is None
    assert all("k" in o for o in out)                        # raw count preserved
