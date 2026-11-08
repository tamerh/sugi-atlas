"""ChIP-Atlas TF targets — pure parsing, genome fallback, and the render block."""
from atlas.gene.sections.s09_tf_regulation import _parse_chip_targets, _int


def test_parse_chip_targets_preserves_source_order():
    raw = ["RPS27L|680.816", "CDKN1A|678.176", "DDB2|559.112"]
    out = _parse_chip_targets(raw)
    assert [t["symbol"] for t in out] == ["RPS27L", "CDKN1A", "DDB2"]
    assert out[0]["score"] == "680.816"


def test_parse_chip_targets_skips_empty_rows():
    assert _parse_chip_targets([]) == []
    assert _parse_chip_targets(["|1.0", "  ", None]) == []


def test_int_tolerates_garbage():
    assert _int("392") == 392
    assert _int(None) == 0
    assert _int(".") == 0


def test_chip_targets_prefers_hg38(monkeypatch):
    from atlas.gene.sections import s09_tf_regulation as m
    def fake_entry(ident, dataset):
        assert ident == "CHIPATLAS:HG38:TP53:1"
        return {"Attributes": {"ChipAtlasTfbs": {
            "tf": "TP53", "genome": "hg38", "top_targets": ["CDKN1A|678.176"],
            "target_count": "13822", "experiment_count": "392", "cell_type_count": "68"}}}
    monkeypatch.setattr(m, "entry", fake_entry)
    out = m._chip_atlas_targets("TP53")
    assert out["chip_atlas_genome"] == "hg38"
    assert out["chip_atlas_experiment_count"] == 392
    assert out["chip_atlas_targets"] == [{"symbol": "CDKN1A", "score": "678.176"}]


def test_chip_targets_falls_back_to_hg19(monkeypatch):
    from atlas.gene.sections import s09_tf_regulation as m
    seen = []
    def fake_entry(ident, dataset):
        seen.append(ident)
        if ident.startswith("CHIPATLAS:HG38"):
            raise RuntimeError("no hg38 entry")
        return {"Attributes": {"ChipAtlasTfbs": {
            "tf": "TP53", "genome": "hg19", "top_targets": ["BAX|351.778"],
            "target_count": "10", "experiment_count": "5", "cell_type_count": "3"}}}
    monkeypatch.setattr(m, "entry", fake_entry)
    out = m._chip_atlas_targets("TP53")
    assert seen == ["CHIPATLAS:HG38:TP53:1", "CHIPATLAS:HG19:TP53:1"]
    assert out["chip_atlas_genome"] == "hg19"


def test_chip_targets_absent_returns_empty(monkeypatch):
    from atlas.gene.sections import s09_tf_regulation as m
    monkeypatch.setattr(m, "entry", lambda ident, dataset: {"Attributes": {}})
    assert m._chip_atlas_targets("NOTATF") == {}


def test_render_block_and_elide():
    from atlas.gene.render import r_tf_regulation
    md = r_tf_regulation({
        "symbol": "TP53", "is_transcription_factor": True,
        "chip_atlas_genome": "hg38", "chip_atlas_target_count": 13822,
        "chip_atlas_experiment_count": 392, "chip_atlas_cell_type_count": 68,
        "chip_atlas_targets": [{"symbol": "CDKN1A", "score": "678.176"}],
    })
    assert "### ChIP-Atlas experimental targets {#chip-atlas-targets}" in md
    assert "CDKN1A" in md and "13,822" in md
    # no targets → no block
    assert "chip-atlas-targets" not in r_tf_regulation({"symbol": "X"})
