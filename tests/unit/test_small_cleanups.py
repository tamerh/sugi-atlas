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


def test_biogrid_human_partners_either_side_ranked(monkeypatch):
    # Projection has only interactor_b_symbol; the full entry has both sides + organisms.
    recs = {
        "1": ("TP53", "P04637", 9606, "MDM2", "Q00987", 9606, "Affinity Capture-Western"),
        "2": ("RCHY1", "Q96PM5", 9606, "TP53", "P04637", 9606, "Two-hybrid"),     # TP53 is B
        "3": ("TP53", "P04637", 9606, "MDM2", "Q00987", 9606, "Reconstituted Complex"),
        "4": ("TP53", "P04637", 9606, "XRS2", "P33301", 559292, "Synthetic Rescue"),  # yeast
        "5": ("TP53", "P04637", 9606, "TP53", "P04637", 9606, "Two-hybrid"),     # self
    }
    rows = [{"id": k, "interactor_b_symbol": v[3]} for k, v in recs.items()]
    def fake_entry(i, ds):
        a = recs[i]
        return {"Attributes": {"BiogridInteraction": {
            "interactor_a_symbol": a[0], "interactor_a_id": a[1], "interactor_a_organism": a[2],
            "interactor_b_symbol": a[3], "interactor_b_id": a[4], "interactor_b_organism": a[5],
            "experimental_system": a[6]}}}
    monkeypatch.setattr(S8, "entry", fake_entry)
    out = S8._biogrid_partners(rows, "TP53", "P04637")
    assert [(x["partner"], x["records"]) for x in out] == [("MDM2", 2), ("RCHY1", 1)]
    md = R.r_interactions({"biogrid": out, "biogrid_count": 6159})
    assert "BioGRID (6,159 interactions; top human partners by supporting records)" in md
    assert "MDM2 (Affinity Capture-Western; 2 records)" in md and "RCHY1 (Two-hybrid)" in md
    md = R.r_interactions({"biogrid": out, "biogrid_count": 6159, "biogrid_sampled": 600})
    assert "by supporting records in a sample of 600 records)" in md


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


def test_no_none_for_missing_transcript_or_uniprot():
    md = R.r_transcripts({"refseq_mrna_count": 0, "canonical_transcript": None})
    assert "`None`" not in md and "{#canonical-exons}" not in md
    md = R.r_protein_ids({"canonical_uniprot": None, "reviewed_uniprot": []})
    assert "`None`" not in md and "No reviewed (Swiss-Prot) UniProt entry" in md
