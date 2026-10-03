"""Gene §drugs: CIViC AMP/ASCO/CAP tier assertions — net-new vs the per-item
civic_evidence (the amp_category tier is the gold-standard actionability call)."""
from atlas.gene.render import r_drugs


def test_civic_assertions_block_renders():
    b = {"symbol": "BRAF", "molecules": [], "molecule_count": 0,
         "civic_assertions": [
             {"profile": "BRAF V600E", "disease": "Melanoma", "type": "Predictive",
              "tier": "Tier I - Level A", "significance": "Sensitivity/Response"},
         ]}
    md = r_drugs(b)
    assert "Clinical-actionability assertions (CIViC AMP/ASCO/CAP) {#civic-assertions}" in md
    assert "BRAF V600E" in md and "Tier I - Level A" in md and "Sensitivity/Response" in md


def test_civic_assertions_elides_when_empty():
    assert "civic-assertions" not in r_drugs(
        {"symbol": "X", "civic_assertions": [], "molecules": [], "molecule_count": 0})
