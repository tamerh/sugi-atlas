"""ClinVar per-transcript consequence + protein-alphaMissense join (pure)."""
from atlas.gene.sections.s06_variants import (
    _protein_substitution, _mane_transcript, _enrich_clinvar_rows,
)


def test_protein_substitution_three_letter_to_one():
    assert _protein_substitution("NP_000537.3:p.Glu294Lys") == ("NP_000537", "E294K")
    assert _protein_substitution("p.Arg175His") == (None, "R175H")
    assert _protein_substitution("p.(Arg175His)") == (None, "R175H")


def test_protein_substitution_rejects_non_missense():
    assert _protein_substitution("NP_000537.3:p.Arg175del") == (None, None)
    assert _protein_substitution("NP_000537.3:p.Arg175fs") == (None, None)
    assert _protein_substitution("p.=") == (None, None)
    assert _protein_substitution(None) == (None, None)


def test_mane_transcript_prefers_mane_select():
    attrs = {"hgvs_transcripts": [
        {"transcript": "NM_1", "consequence": "intron variant"},
        {"transcript": "NM_2", "consequence": "missense variant", "mane_select": True},
    ]}
    assert _mane_transcript(attrs)["transcript"] == "NM_2"
    assert _mane_transcript({}) is None


def test_enrich_clinvar_rows_joins_consequence_and_alphamissense(monkeypatch):
    from atlas.gene.sections import s06_variants as m

    def fake_entry(ident, dataset):
        if dataset == "clinvar":
            return {"Attributes": {"Clinvar": {"hgvs_transcripts": [
                {"transcript": "NM_000546.6", "consequence": "missense variant",
                 "hgvs_p": "NP_000537.3:p.Glu294Lys", "mane_select": True}]}}}
        if dataset == "alphamissense_aa":
            assert ident == "P04637:E294K"
            return {"Attributes": {"Alphamissense": {
                "am_pathogenicity": 0.98, "am_class": "pathogenic"}}}
        return {"Attributes": {}}

    def fake_map(ids, chain, cap=60):
        assert chain == ">>refseq>>uniprot" and ids == "NP_000537"
        return [{"id": "P04637"}]

    monkeypatch.setattr(m, "entry", fake_entry)
    monkeypatch.setattr(m, "map_all", fake_map)
    rows = [{"id": "123"}]
    m._enrich_clinvar_rows(rows)
    assert rows[0]["consequence"] == "missense variant"
    assert rows[0]["am_score"] == 0.98
    assert rows[0]["am_class"] == "pathogenic"


def test_variants_table_has_new_columns():
    from atlas.gene.render import r_variants
    md = r_variants({"symbol": "TP53", "top_pathogenic_total": 1, "top_pathogenic": [{
        "id": "123", "hgvs": "NM_000546.6(TP53):c.880G>A (p.Glu294Lys)",
        "consequence": "missense variant", "am_score": 0.98, "am_class": "pathogenic",
        "classification": "Pathogenic", "review_status": "reviewed by expert panel",
    }]})
    assert "Consequence" in md and "AlphaMissense" in md
    assert "missense variant" in md and "0.98 (pathogenic)" in md


def test_transcripts_render_mane_version():
    from atlas.gene.render import r_transcripts
    md = r_transcripts({"refseq_mrna_count": 1, "refseq_mrna": ["NM_000546"],
                        "mane_select_refseq": "NM_000546",
                        "mane_select_refseq_version": "NM_000546.6",
                        "mane_select_refseq_ensembl": "ENST00000269305.9"})
    assert "NM_000546.6" in md and "ENST00000269305.9" in md
