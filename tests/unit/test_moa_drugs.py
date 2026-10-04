"""Gene §drugs #chembl-moa: curated ChEMBL mechanism-of-action drugs whose MOA
target is the gene's OWN protein. Guards the v1.11.7 QA findings: family/complex
targets must not leak a whole drug class onto every member gene; curated small
molecules must not vanish (no dedupe against the capped bioactivity table); salt
forms collapse; phase -1 renders as "—"."""
from types import SimpleNamespace

from atlas.gene.render import r_drugs
from atlas.gene.sections import s10_drugs as S


def _mol(i, name, typ="Small molecule", ph="4"):
    return {"id": i, "name": name, "type": typ, "highestDevelopmentPhase": ph}


def _fake(monkeypatch, targets, per_target, direct=(), bioactivity=()):
    def fake_map_all(root, chain, **k):
        if chain == ">>uniprot>>chembl_target":
            return targets
        if "chembl_target>>chembl_molecule[highestDevelopmentPhase" in chain:
            return list(bioactivity) if root == targets[0]["id"] else []
        if chain == ">>chembl_target>>chembl_mechanism>>chembl_molecule":
            return per_target.get(root, [])
        if chain == ">>hgnc>>chembl_mechanism>>chembl_molecule":
            return list(direct)
        return []
    monkeypatch.setattr(S, "map_all", fake_map_all)
    return S.collect(SimpleNamespace(symbol="EGFR", canonical_uniprot="P00533", hgnc_id="HGNC:3236"))


def test_family_and_complex_targets_excluded(monkeypatch):
    targets = [{"id": "T_SP", "title": "EGFR", "type": "SINGLE PROTEIN"},
               {"id": "T_FAM", "title": "ErbB family", "type": "PROTEIN FAMILY"},
               {"id": "T_CPX", "title": "Some complex", "type": "PROTEIN COMPLEX GROUP"}]
    b = _fake(monkeypatch, targets, {
        "T_SP": [_mol("C1", "GEFITINIB"), _mol("C2", "CETUXIMAB", "Antibody")],
        "T_FAM": [_mol("C9", "FAMILYDRUG")],
        "T_CPX": [_mol("C8", "COMPLEXDRUG")],
    })
    names = {d["name"] for d in b["moa_drugs"]}
    assert names == {"GEFITINIB", "CETUXIMAB"}       # family/complex drugs excluded


def test_not_deduped_against_bioactivity_and_salts_collapse(monkeypatch):
    targets = [{"id": "T_SP", "title": "EGFR", "type": "SINGLE PROTEIN"}]
    b = _fake(monkeypatch, targets, {"T_SP": [
        _mol("C1", "GEFITINIB"),
        _mol("C3", "ABIVERTINIB", ph="3"), _mol("C4", "ABIVERTINIB MALEATE", ph="3"),
        _mol("C5", None), _mol("C6", "CHEMBL999"),
    ]}, bioactivity=[_mol("C1", "GEFITINIB")])
    names = [d["name"] for d in b["moa_drugs"]]
    assert "GEFITINIB" in names                       # kept even though it's in the bioactivity list
    assert names.count("ABIVERTINIB") == 1 and "ABIVERTINIB MALEATE" not in names
    assert "CHEMBL999" not in names and None not in names


def test_direct_route_only_adds_oligonucleotides(monkeypatch):
    targets = [{"id": "T_SP", "title": "PCSK9", "type": "SINGLE PROTEIN"}]
    b = _fake(monkeypatch, targets, {"T_SP": [_mol("A1", "EVOLOCUMAB", "Antibody")]},
              direct=[_mol("O1", "INCLISIRAN", "Oligonucleotide", "3"),
                      _mol("X1", "NOTANOLIGO", "Small molecule")])
    names = {d["name"] for d in b["moa_drugs"]}
    assert names == {"EVOLOCUMAB", "INCLISIRAN"}


def test_salt_base():
    assert S._salt_base("ELACESTRANT HYDROCHLORIDE") == "elacestrant"
    assert S._salt_base("PENTOBARBITAL SODIUM") == "pentobarbital"
    assert S._salt_base("SODIUM CHLORIDE") == "sodium"     # never strips the last token
    assert S._salt_base("CERTOLIZUMAB PEGOL") == "certolizumab pegol"


def test_render_moa_block_first_with_dash_phase():
    b = {"symbol": "EGFR", "molecule_count": 1,
         "molecules": [{"id": "C1", "name": "GEFITINIB", "phase": "4"}],
         "moa_drugs": [{"id": "C2", "name": "CETUXIMAB", "type": "Antibody", "phase": "4"},
                       {"id": "C7", "name": "OLDTHING", "type": "Small molecule", "phase": "-1"}]}
    md = r_drugs(b)
    assert "Drugs by curated mechanism of action (ChEMBL) {#chembl-moa}" in md
    assert "single-protein targets only" in md
    assert md.index("{#chembl-moa}") < md.index("{#chembl-molecules}")   # MOA leads
    assert "| OLDTHING | Small molecule | — |" in md                   # phase -1 → —


def test_moa_block_elides_when_empty():
    assert "chembl-moa" not in r_drugs({"symbol": "X", "moa_drugs": [], "molecules": [], "molecule_count": 0})
