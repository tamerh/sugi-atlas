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


def test_civic_assertions_therapy_and_direction():
    b = {"symbol": "EGFR", "molecules": [], "molecule_count": 0,
         "civic_assertions": [
             {"profile": "EGFR L858R", "disease": "NSCLC", "type": "Predictive",
              "tier": "Tier I - Level A", "significance": "Sensitivity/Response",
              "therapies": ["Gefitinib"], "direction": "Supports", "companion_test": True},
             {"profile": "EGFR T790M", "disease": "NSCLC", "type": "Predictive",
              "tier": "Tier I - Level A", "significance": "Resistance",
              "therapies": ["Erlotinib"], "direction": "Does Not Support", "companion_test": False},
         ]}
    md = r_drugs(b)
    assert "| Gefitinib † |" in md and "companion diagnostic" in md
    assert "Resistance (does not support)" in md


def test_civic_wrong_gene_assertion_dropped(monkeypatch):
    # biobtree's hgnc→civic link pointed FDXR at CIViC gene AR — must be filtered.
    from types import SimpleNamespace
    from atlas.gene.sections import s10_drugs as S

    def fake_map_all(root, chain, **k):
        if chain == ">>hgnc>>civic>>civic_assertion":
            return [{"id": "1", "molecular_profile": "AR AR-V7", "disease": "Prostate",
                     "assertion_type": "Predictive", "amp_category": "Tier I - Level A",
                     "significance": "Resistance"}]
        return []
    monkeypatch.setattr(S, "map_all", fake_map_all)
    monkeypatch.setattr(S, "entry", lambda *a, **k: {})
    b = S.collect(SimpleNamespace(symbol="FDXR", canonical_uniprot=None, hgnc_id="HGNC:3642"))
    assert b["civic_assertions"] == []
    # …but a fusion profile still matches its gene
    assert S.re.compile(r"(?<![A-Za-z0-9-])ALK(?![A-Za-z0-9-])").search("EML4::ALK Fusion")
