"""v1.11.8 small page cleanups: MANE 'None', '0 indication records', BioGRID
self-partners, HPA Brain Atlas sub-regions flooding the tissue table."""
from types import SimpleNamespace

from atlas.drug.render import r_indications
from atlas.gene import render as R
from atlas.gene.sections import s08_interactions as S8


def test_mane_none_not_rendered():
    md = R.r_transcripts({"refseq_mrna_count": 0, "mane_select_refseq": "None"})
    assert "MANE Select" not in md
    md = R.r_transcripts({"refseq_mrna_count": 1, "mane_select_refseq": "NM_000546.6",
                          "refseq_mrna": ["NM_000546.6"]})
    assert "MANE Select: `NM_000546.6`" in md


def test_zero_indications_plain_sentence():
    md = r_indications({"indication_count": 0, "indications": []})
    assert "0 indication record" not in md
    assert "No ChEMBL drug-indication records" in md


def test_biogrid_drops_self_and_dedupes(monkeypatch):
    rows = [{"interactor_b_symbol": s, "experimental_system": m} for s, m in [
        ("TP53", "Affinity Capture-Western"), ("MDM2", "Affinity Capture-Western"),
        ("tp53", "Reconstituted Complex"), ("MDM2", "Two-hybrid"), ("RCHY1", "Two-hybrid")]]
    monkeypatch.setattr(S8, "map_all",
                        lambda r, c, **k: rows if "biogrid" in c else [])
    monkeypatch.setattr(S8, "entry", lambda *a, **k: {}, raising=False)
    b = S8.collect(SimpleNamespace(symbol="TP53", canonical_uniprot="P04637", hgnc_id="HGNC:11998"))
    assert [x["partner"] for x in b["biogrid"]] == ["MDM2", "RCHY1"]
    assert b["biogrid"][0]["method"] == "Affinity Capture-Western"


def test_hpa_brain_regions_collapsed():
    exp = ([{"entity": f"Nucleus {i}", "axis": "tissue", "ntpm": str(20000 - i), "brain_region": True}
            for i in range(80)]
           + [{"entity": "Liver", "axis": "tissue", "ntpm": "5.0"},
              {"entity": "Astrocytes", "axis": "cell", "ntpm": "900.0"}])
    exp.sort(key=lambda r: -float(r["ntpm"]))
    md = R.r_hpa_expression({"13": {"hpa_expression": exp, "hpa": {"rna_tissue_specificity": "Tissue enriched"}}})
    assert "| Liver | tissue | 5.0 |" in md and "| Astrocytes | cell | 900.0 |" in md
    assert "| Nucleus 0 |" not in md
    assert "**Brain sub-regions (HPA Brain Atlas):** detected in 80; highest nTPM — Nucleus 0 (20000)" in md
