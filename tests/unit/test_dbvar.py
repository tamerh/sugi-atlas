"""dbVar structural-variant section — pure helpers (no network).

The collector's clustering is what makes the section readable: dbVar is
per-submission, so 1,066 Duchenne records are 685 overlapping deletions of one
locus. These tests pin the merge rules, the significance majority, and the
insertion-size fallback.
"""
from types import SimpleNamespace

from atlas.disease.sections.s03_variant_details import (
    _dbvar_rank, _dbvar_cluster, _collect_dbvar, dbvar_size,
)


def _rec(chrom="7", start=100, end=200, svtype="DEL", sig="Pathogenic", svlen=None):
    return {"id": f"{chrom}:{start}-{end}:{svtype}", "chrom": chrom,
            "start": start, "end": end, "svtype": svtype, "sig": sig,
            "svlen": svlen, "clinvar": None}


def test_rank_orders_clinical_significance():
    assert _dbvar_rank("Pathogenic") == 0
    assert _dbvar_rank("Likely pathogenic") == 1
    assert _dbvar_rank("Conflicting interpretations") == 3
    assert _dbvar_rank("Uncertain significance") == 4
    assert _dbvar_rank("Benign") == 5
    assert _dbvar_rank("not provided") == 6
    assert _dbvar_rank(None) == 4


def test_size_formatting():
    assert dbvar_size(0) == "—"
    assert dbvar_size(500) == "500 bp"
    assert dbvar_size(62_000) == "62 kb"
    assert dbvar_size(1_100_000) == "1.1 Mb"


def test_cluster_merges_overlapping_same_type():
    loci = _dbvar_cluster([_rec(start=100, end=200), _rec(start=150, end=300)])
    assert len(loci) == 1
    assert loci[0]["n"] == 2
    assert (loci[0]["start"], loci[0]["end"]) == (100, 300)


def test_cluster_never_merges_across_chromosome_or_type():
    recs = [_rec(chrom="7", start=100, end=200),
            _rec(chrom="8", start=100, end=200),          # different chromosome
            _rec(chrom="7", start=100, end=200, svtype="DUP")]  # different type
    assert len(_dbvar_cluster(recs)) == 3


def test_cluster_requires_reciprocal_overlap():
    # 1 kb vs 10 kb overlapping by 0.5 kb → 0.5/1 = 50% (merges at the boundary)
    a = _rec(start=0, end=1_000)
    b = _rec(start=500, end=10_500)
    assert len(_dbvar_cluster([a, b])) == 1
    # overlap 0.4 kb / 1 kb = 40% → stays separate
    c = _rec(start=600, end=10_600)
    assert len(_dbvar_cluster([a, c])) == 2


def test_cluster_majority_significance_or_mixed():
    majority = _dbvar_cluster([_rec(start=0, end=100, sig="Pathogenic"),
                               _rec(start=10, end=110, sig="Pathogenic"),
                               _rec(start=20, end=120, sig="Uncertain significance")])
    assert majority[0]["sig"] == "Pathogenic"
    mixed = _dbvar_cluster([_rec(start=0, end=100, sig="Pathogenic"),
                            _rec(start=10, end=110, sig="Uncertain significance"),
                            _rec(start=20, end=120, sig="Uncertain significance")])
    assert mixed[0]["sig"].startswith("mixed (1/3")


def test_insertion_uses_svlen_when_start_equals_end():
    loci = _dbvar_cluster([_rec(start=500, end=500, svtype="INS", svlen="115")])
    assert loci[0]["size"] == 115


def test_no_dbvar_xref_makes_no_call():
    # xref_counts says 0 → the collector returns empty without touching biobtree
    a = SimpleNamespace(mondo_id="MONDO:0000001", xref_counts={})
    out = _collect_dbvar(a)
    assert out["dbvar_total"] == 0 and out["dbvar_loci"] == []


def test_render_block_and_empty_placeholder():
    from atlas.disease.render import r_variant_details
    bundle = {"dbvar_total": 2, "dbvar_retrieved": 2, "dbvar_plp": 1,
              "dbvar_svtype_counts": {"DEL": 2},
              "dbvar_loci": [{"chrom": "7", "start": 1, "end": 100, "svtype": "DEL",
                              "n": 2, "sig": "Pathogenic", "clinvar": "RCV1",
                              "size": 99}]}
    md = r_variant_details(bundle)
    assert "### Structural variants (dbVar) {#structural-variants}" in md
    assert "| DEL | 7:1-100 | 99 bp | Pathogenic | 2 | RCV1 |" in md
    # no dbVar and no other variant data → placeholder mentions all three, no block
    empty = r_variant_details({})
    assert "structural-variants" not in empty
    assert "structural variants" in empty
