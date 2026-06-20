"""GO ECO evidence + Reactome evidence/disease tiering (gene §7). Pure helpers;
the live collect() hits biobtree and is exercised by the corpus build."""
from atlas.gene.sections.s07_pathways import (
    _ECO_EXPERIMENTAL, _ECO_LABEL, _eco_tier, _merge_reactome,
)


# ── the invariant that would have caught the v1.9.0 raw-ECO leak ──────────────
def test_every_experimental_eco_has_a_display_label():
    # An ECO code we count as "experimental" (so it's shown) MUST have a readable
    # label, else the render falls back to the raw "ECO:000…" id (the leak bug).
    missing = [eco for eco in _ECO_EXPERIMENTAL if eco not in _ECO_LABEL]
    assert not missing, f"experimental ECO codes missing a label: {missing}"


def test_eco_label_values_are_not_raw_codes():
    for eco, label in _ECO_LABEL.items():
        assert not label.startswith("ECO:"), f"{eco} maps to a raw-looking label"


# ── evidence tiering ─────────────────────────────────────────────────────────
def test_eco_tier_buckets():
    assert _eco_tier("ECO:0000314") == "experimental"   # direct assay
    assert _eco_tier("ECO:0007005") == "experimental"   # high-throughput assay
    assert _eco_tier("ECO:0000501") == "electronic"     # IEA
    assert _eco_tier("ECO:0000318") == "computational"  # phylogenetic
    assert _eco_tier("ECO:0000303") == "author/curator"
    assert _eco_tier(None) == "other"


# ── Reactome evidence merge: TAS (curated) beats IEA (electronic) ─────────────
def test_reactome_merge_tas_beats_iea_regardless_of_order():
    for first, second in (("IEA", "TAS"), ("TAS", "IEA")):
        rx = {}
        _merge_reactome(rx, {"id": "R-HSA-1", "name": "P", "evidence": first})
        _merge_reactome(rx, {"id": "R-HSA-1", "name": "P", "evidence": second})
        assert rx["R-HSA-1"]["evidence"] == "curated"


def test_reactome_merge_maps_codes_and_unions_disease_flag():
    rx = {}
    _merge_reactome(rx, {"id": "R-HSA-2", "name": "Q", "evidence": "IEA",
                         "is_disease_pathway": "false"})
    assert rx["R-HSA-2"]["evidence"] == "electronic"
    assert rx["R-HSA-2"]["is_disease"] is False
    # a later disease-flagged row flips it on (union)
    _merge_reactome(rx, {"id": "R-HSA-2", "name": "Q", "evidence": "TAS",
                         "is_disease_pathway": "true"})
    assert rx["R-HSA-2"]["is_disease"] is True
    assert rx["R-HSA-2"]["evidence"] == "curated"


def test_reactome_merge_skips_idless_rows():
    rx = {}
    _merge_reactome(rx, {"name": "no id", "evidence": "TAS"})
    assert rx == {}
