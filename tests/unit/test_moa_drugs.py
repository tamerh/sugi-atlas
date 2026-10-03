"""Gene §drugs: the curated mechanism-of-action subsection surfaces the modality gap
(antibodies / ADCs / oligonucleotides the bioactivity table misses, e.g. cetuximab→
EGFR, inclisiran→PCSK9) — and is net-new, not a duplicate of the bioactivity list."""
from atlas.gene.render import r_drugs
from atlas.gene.sections import s10_drugs as S


def test_moa_subsection_renders_modality_and_links():
    b = {"symbol": "EGFR", "molecules": [], "molecule_count": 0,
         "moa_drugs": [
             {"id": "CHEMBL1201585", "name": "CETUXIMAB", "type": "Antibody", "phase": "4"},
             {"id": "CHEMBL3990033", "name": "INCLISIRAN", "type": "Oligonucleotide", "phase": "4"},
         ]}
    md = r_drugs(b)
    assert "Targeted drugs by mechanism of action (ChEMBL) {#chembl-moa}" in md
    assert "CETUXIMAB" in md and "Antibody" in md
    assert "INCLISIRAN" in md and "Oligonucleotide" in md
    assert "does not surface" in md                    # honest "net-new vs bioactivity" framing


def test_moa_subsection_elides_when_empty():
    assert "chembl-moa" not in r_drugs({"symbol": "X", "moa_drugs": [], "molecules": [], "molecule_count": 0})


def test_collect_moa_dedups_and_filters(monkeypatch):
    # Drive collect with canned biobtree rows; everything not listed returns [].
    def fake_map_all(root, chain, **k):
        if chain == ">>uniprot>>chembl_target":
            return [{"id": "CT1", "title": "EGFR", "type": "SINGLE PROTEIN"}]
        if "chembl_target>>chembl_molecule[highestDevelopmentPhase" in chain:
            return [{"id": "CHEMBL939", "name": "GEFITINIB", "type": "Small molecule",
                     "highestDevelopmentPhase": "4"}]                      # already in bioactivity list
        if chain == ">>uniprot>>chembl_target>>chembl_mechanism>>chembl_molecule":
            return [
                {"id": "CHEMBL939", "name": "GEFITINIB", "type": "Small molecule", "highestDevelopmentPhase": "4"},   # dup → drop
                {"id": "CHEMBL1201585", "name": "CETUXIMAB", "type": "Antibody", "highestDevelopmentPhase": "4"},     # new → keep
                {"id": "CX", "name": None, "type": "Antibody", "highestDevelopmentPhase": "3"},                        # unnamed → drop
                {"id": "CY", "name": "CHEMBL999", "type": "Unknown", "highestDevelopmentPhase": "1"},                  # bare-id name → drop
            ]
        if chain == ">>hgnc>>chembl_mechanism>>chembl_molecule":
            return [{"id": "CHEMBL3990033", "name": "INCLISIRAN", "type": "Oligonucleotide",
                     "highestDevelopmentPhase": "3"}]                      # oligo catch → keep
        return []
    monkeypatch.setattr(S, "map_all", fake_map_all)
    from types import SimpleNamespace
    b = S.collect(SimpleNamespace(symbol="EGFR", canonical_uniprot="P00533", hgnc_id="HGNC:3236"))
    names = {d["name"] for d in b["moa_drugs"]}
    assert names == {"CETUXIMAB", "INCLISIRAN"}        # dup + unnamed + bare-id excluded
