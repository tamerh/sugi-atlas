"""ClinVar variant lists rank by review status (gold stars), not biobtree order —
v1.11.7 BRCA1 led its top-pathogenic table with 'no assertion criteria' rows."""
from types import SimpleNamespace

from atlas.gene.render import r_variants
from atlas.gene.sections import s06_variants as S
from atlas.render_common import clinvar_stars


def test_clinvar_stars():
    assert clinvar_stars("practice guideline") == 4
    assert clinvar_stars("Reviewed by expert panel") == 3
    assert clinvar_stars("criteria provided, multiple submitters, no conflicts") == 2
    assert clinvar_stars("criteria provided, single submitter") == 1
    assert clinvar_stars("no assertion criteria provided") == 0
    assert clinvar_stars(None) == 0


def test_gene_top_pathogenic_ranked_and_capped(monkeypatch):
    def row(i, cls, rs):
        return {"id": str(i), "name": f"v{i}", "germline_classification": cls, "review_status": rs}
    P = [row(1, "Pathogenic", "no assertion criteria provided"),
         row(2, "Pathogenic", "criteria provided, single submitter")] + \
        [row(100 + i, "Pathogenic", "no assertion criteria provided") for i in range(20)]
    LP = [row(3, "Likely pathogenic", "reviewed by expert panel"),
          row(4, "Likely pathogenic", "criteria provided, single submitter")]

    def fake_map_all(root, chain, **k):
        if '"Pathogenic"' in chain:
            return P
        if '"Likely pathogenic"' in chain:
            return LP
        return []
    monkeypatch.setattr(S, "map_all", fake_map_all)
    monkeypatch.setattr(S, "xref_counts", lambda e: {})
    b = S.collect(SimpleNamespace(symbol="BRCA1", hgnc_id="HGNC:1100", hgnc_entry={},
                                  canonical_transcript=None))
    ids = [v["id"] for v in b["top_pathogenic"]]
    assert ids[:3] == ["3", "2", "4"]            # expert panel, then 1★ P before 1★ LP
    assert len(ids) == S.TOP_PATHOGENIC_CAP and b["top_pathogenic_total"] == 24
    md = r_variants(b)
    assert "the 15 best-reviewed of 24" in md and "| Review |" in md
