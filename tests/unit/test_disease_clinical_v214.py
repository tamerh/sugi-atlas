"""Disease Clinical-features zone (biobtree v2.14.0): Mondo/DO definitions,
Wikidata common-disease symptoms (labelled, separate from HPO), and no disease
listed as its own HPO feature."""
from atlas.disease import render as DR


def test_mondo_definition_shown_but_not_when_copy_of_orphanet():
    md = DR.r_clinical_description({"mondo_definition": "A progressive, neurodegenerative disease X.",
                                    "mesh_scope_note": "A degenerative disease of the BRAIN."})
    assert "**Mondo:** A progressive" in md and "**MeSH:**" in md
    orph = "Rare syndrome characterized by short stature and intellectual disability in children."
    md = DR.r_clinical_description({"orphanet_definition": orph, "mondo_definition": orph + " (Orphanet)"})
    assert md.count("short stature") == 1


def test_doid_definition_only_as_fallback():
    assert "Disease Ontology" not in DR.r_clinical_description(
        {"mondo_definition": "A disease.", "doid_definition": "A tauopathy characterized by memory lapses."})
    assert "**Disease Ontology:** A tauopathy" in DR.r_clinical_description(
        {"doid_definition": "A tauopathy characterized by memory lapses."})


def test_wikidata_symptoms_block_separate_and_labelled():
    w = [{"qid": "Q86", "name": "headache"}, {"qid": "Q186889", "name": "nausea"}]
    md = DR.r_symptoms({"wikidata_symptoms": w})
    assert "{#wikidata-symptoms}" in md and "{#hpo-features}" not in md
    assert "[headache](https://www.wikidata.org/wiki/Q86)" in md
    assert "crowd-curated" in md and "2 symptoms and signs" in md
    assert "1 symptom or sign listed" in DR.r_symptoms({"wikidata_symptoms": w[:1]})
    md = DR.r_symptoms({"phenotypes": [{"hpo_id": "HP:0002315", "hpo_term": "Headache"}],
                        "wikidata_symptoms": w})
    assert md.index("{#hpo-features}") < md.index("{#wikidata-symptoms}")
    assert DR.r_symptoms({}) == ""


def test_collector_symptoms_via_mondo_and_doid_and_self_guard(monkeypatch):
    from types import SimpleNamespace
    from atlas.disease.sections import s01_disease_ids as S
    def fake_map_all(root, chain, **k):
        if chain == ">>mondo>>doid":
            return [{"id": "DOID:9352"}]
        if chain == ">>mondo>>wikidata_symptom":
            return [{"id": "Q1", "name": "fatigue"}]
        if chain == ">>doid>>wikidata_symptom":
            return [{"id": "Q1", "name": "fatigue"}, {"id": "Q2", "name": "polyuria"}]
        if chain == ">>mondo>>hpo":
            return [{"id": "HP:0005978", "name": "Type II diabetes mellitus"},
                    {"id": "HP:0003074", "name": "Hyperglycemia"}]
        return []
    def fake_entry(i, ds):
        if ds == "doid":
            return {"Attributes": {"Ontology": {"definition": "A diabetes characterized by ..."}}}
        return {}
    monkeypatch.setattr(S, "map_all", fake_map_all)
    monkeypatch.setattr(S, "entry", fake_entry)
    a = SimpleNamespace(mondo_id="MONDO:0005148", name="type 2 diabetes mellitus",
                        canonical_name="type 2 diabetes mellitus", synonyms=(), efo_id=None,
                        mesh_ids=(), omim_ids=(), orphanet_ids=(), orphanet_attrs={}, obo_xrefs={},
                        anatomy_uberon_ids=(), xref_counts={}, is_cancer=False,
                        mondo_entry={"Attributes": {"Ontology": {"definition": "A type of diabetes."}}})
    b = S.collect(a)
    assert [x["name"] for x in b["wikidata_symptoms"]] == ["fatigue", "polyuria"]   # deduped
    assert b["mondo_definition"] == "A type of diabetes." and b["doid_definition"].startswith("A diabetes")
    assert [p["hpo_term"] for p in b["phenotypes"]] == ["Hyperglycemia"]          # self-row dropped


def test_biogrid_uses_projection_without_entry_calls(monkeypatch):
    from atlas.gene.sections import s08_interactions as S8
    monkeypatch.setattr(S8, "entry", lambda *a, **k: (_ for _ in ()).throw(AssertionError("no entry calls")))
    rows = [{"id": "1", "interactor_a_symbol": "RCHY1", "interactor_a_organism": "9606",
             "interactor_b_symbol": "TP53", "interactor_b_organism": "9606", "experimental_system": "Two-hybrid"},
            {"id": "2", "interactor_a_symbol": "TP53", "interactor_a_organism": "9606",
             "interactor_b_symbol": "XRS2", "interactor_b_organism": "559292", "experimental_system": "Synthetic Rescue"}]
    out = S8._biogrid_partners(rows, "TP53", "P04637")
    assert [x["partner"] for x in out] == ["RCHY1"]


def test_same_text_tolerates_spelling_variants():
    o = "Lissencephaly syndrome, Norman-Roberts type is characterised by the association of lissencephaly type I with craniofacial anomalies."
    m = o.replace("characterised", "characterized")
    assert DR._same_text(m, o)
    assert not DR._same_text("A progressive neurodegenerative disease of the brain.",
                             "A degenerative disease of the BRAIN characterized by dementia.")
