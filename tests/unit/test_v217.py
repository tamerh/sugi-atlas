"""biobtree v2.17 integration: DrugCentral label indications + ATC drug class,
HPO inheritance fallback, MIR targets, MIM names."""


def test_label_match_only_qualifiers_may_differ():
    from atlas.drug.anchors import label_match
    assert label_match("chronic myeloid leukemia", ["Chronic phase chronic myeloid leukemia"])
    assert label_match("type 2 diabetes mellitus", ["Diabetes mellitus type 2"])
    assert not label_match("hypertension", ["Pulmonary hypertension"])
    assert not label_match("cancer", ["Breast cancer"])
    assert not label_match("bacterial infectious disease", ["Bacterial infections"])


def test_label_rows_marked_and_added(monkeypatch):
    from atlas.drug import anchors as A
    monkeypatch.setattr(A, "map_all", lambda i, chain, cap=None:
                        [{"id": "MONDO:1", "name": "chronic myeloid leukemia"}] if i == "EFO:1" else [])
    label = {"names": ["Chronic phase chronic myeloid leukemia", "Systemic mast cell disease"],
             "mondo": {"MONDO:9": "systemic mastocytosis"}, "atc_classes": []}
    recs = A._resolve_indications({"indications": [{"efo": "EFO:1", "highestDevelopmentPhase": 3}]},
                                  (), label)
    by = {r.name: r for r in recs}
    assert by["chronic myeloid leukemia"].label and by["chronic myeloid leukemia"].max_phase == 3
    assert by["systemic mastocytosis"].label and by["systemic mastocytosis"].max_phase == 0


def test_label_rows_are_approved_and_rendered_with_source():
    from atlas.indication import indication_tiers
    from atlas.drug.render import r_indications
    b4 = {"label_indications": ["Chronic phase chronic myeloid leukemia"], "indications": [
        {"name": "chronic myeloid leukemia", "mondo_id": "MONDO:1", "max_phase": 3,
         "approved": True, "label": True, "late": False},
        {"name": "stroke", "mondo_id": "MONDO:3", "max_phase": 2, "approved": False, "late": False}]}
    appr, late, trials = indication_tiers(b4)
    assert [i["name"] for i in appr] == ["chronic myeloid leukemia"] and not late
    md = r_indications(b4)
    assert "**On-label indications (DrugCentral):** Chronic phase chronic myeloid leukemia." in md
    assert "| DrugCentral label |" in md and "(DrugCentral on-label use or ChEMBL phase 4)" in md


def test_atc_drug_class():
    from atlas.drug.roles import drug_class
    split = lambda s: s.split(";")
    assert drug_class(["non-steroidal anti-inflammatory drug"],
                      split("ANTIPROTOZOALS;ANTIMALARIALS;Aminoquinolines")) == "antimalarial"
    assert drug_class(["tyrosine kinase inhibitor"],
                      split("ANTINEOPLASTIC AGENTS;OTHER ANTINEOPLASTIC AGENTS;Protein kinase inhibitors")) \
        == "tyrosine kinase inhibitor"
    assert drug_class(["P450 inhibitor"], split(
        "ANTIFUNGALS FOR DERMATOLOGICAL USE;ANTIFUNGALS FOR TOPICAL USE;Imidazole and triazole derivatives;"
        "ANTIMYCOTICS FOR SYSTEMIC USE;Triazole derivatives")) == "antimycotic"
    assert drug_class([], split("LIPID MODIFYING AGENTS;LIPID MODIFYING AGENTS, PLAIN;"
                                "HMG CoA reductase inhibitors")) == \
        "lipid modifying agent (HMG CoA reductase inhibitor)"
    assert drug_class([], split("AGENTS ACTING ON THE RENIN-ANGIOTENSIN SYSTEM;ANGIOTENSIN II RECEPTOR "
                                "BLOCKERS (ARBs), PLAIN;Angiotensin II receptor blockers (ARBs), plain")) == \
        "angiotensin II receptor blocker"
    assert drug_class(["anticonvulsant"], [], "drug") == "anticonvulsant"     # no ATC: ChEBI


def test_hpo_inheritance_fallback_when_no_orphanet():
    from atlas.page.disease_declarative import inheritance_modes, declarative_sentence
    b = {"1": {"canonical_name": "cystic fibrosis", "mondo_id": "MONDO:1",
               "hpo_inheritance": ["Autosomal recessive inheritance"]},
         "2": {}, "4": {}, "5": {}}
    assert inheritance_modes(b) == (["Autosomal recessive inheritance"], "HPO")
    assert "inherited in an autosomal recessive pattern" in declarative_sentence(b)


def test_disease_ids_show_mim_names():
    from atlas.disease.render import r_disease_ids
    md = r_disease_ids({"canonical_name": "cystic fibrosis", "omim_ids": ["219700"],
                        "omim_names": {"219700": "Cystic fibrosis"}})
    assert "219700 (Cystic fibrosis)" in md


def test_mir_gene_targets_block():
    from atlas.gene.render import r_tf_regulation
    md = r_tf_regulation({"symbol": "MIR21", "mirna_products": [
        {"id": "hsa-miR-21-5p", "target_count": 1035, "targets": [{"symbol": "PDCD4", "score": "99.0"}]}]})
    assert "{#mirdb-targets}" in md and "**hsa-miR-21-5p** (1,035 predicted targets): PDCD4 (99.0)" in md


def test_atc_class_edge_cases():
    from atlas.drug.roles import drug_class
    split = lambda s: s.split(";")
    assert drug_class(["bronchodilator agent"], split(
        "DRUGS FOR OBSTRUCTIVE AIRWAY DISEASES;OTHER DRUGS FOR OBSTRUCTIVE AIRWAY DISEASES, INHALANTS;"
        "Anticholinergics")) == "bronchodilator agent"                     # generic ATC → ChEBI
    assert drug_class(["non-narcotic analgesic"], split(
        "ANTINEOPLASTIC AGENTS;OTHER ANTINEOPLASTIC AGENTS;Other antineoplastic agents;"
        "ANTIINFLAMMATORY AND ANTIRHEUMATIC PRODUCTS;ANTIINFLAMMATORY AND ANTIRHEUMATIC PRODUCTS, "
        "NON-STEROIDS;Coxibs")) == "non-steroidal anti-inflammatory drug"
    assert drug_class(["antilipemic drug"], split(
        "OTHER ALIMENTARY TRACT AND METABOLISM PRODUCTS;Amino acids and derivatives")) == "antilipemic drug"
    assert drug_class(["anti-allergic agent"], split(
        "STOMATOLOGICAL PREPARATIONS;Other agents for local oral treatment")) == "anti-allergic agent"
    assert drug_class(["antineoplastic agent"], split(
        "ENDOCRINE THERAPY;HORMONE ANTAGONISTS AND RELATED AGENTS;Anti-androgens")) == "hormone antagonist"


def test_hpo_inheritance_not_used_for_complex_or_umbrella():
    from atlas.page.disease_declarative import inheritance_modes
    b = {"1": {"canonical_name": "type 2 diabetes mellitus", "mondo_id": "MONDO:1",
               "hpo_inheritance": ["Autosomal dominant inheritance"]},
         "2": {"assoc_total": 8000}, "4": {}, "5": {}}
    assert inheritance_modes(b) == ([], "")
