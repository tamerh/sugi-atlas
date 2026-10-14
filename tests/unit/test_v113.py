"""v1.13: three-tier indications, molecule-level approval in FAQ, and the
precision fixes from the v1.12.2 audit."""


def test_three_tier_indications_lead_table_faq():
    from atlas.page.drug_declarative import _indications_clause
    from atlas.drug.render import r_indications
    from atlas.page.faq_intents import intent_qas
    b4 = {"indication_count": 3, "indications": [
        {"name": "chronic myeloid leukemia", "mondo_id": "MONDO:1", "max_phase": 3, "approved": False, "late": True},
        {"name": "gastrointestinal stromal tumor", "mondo_id": "MONDO:2", "max_phase": 3, "approved": False, "late": True},
        {"name": "stroke", "mondo_id": "MONDO:3", "max_phase": 2, "approved": False, "late": False}]}
    s = _indications_clause(b4)
    assert s == ("; its most advanced recorded uses (phase 3) include chronic myeloid leukemia "
                 "and gastrointestinal stromal tumor, with clinical trials for 1 other condition")
    md = r_indications(b4)
    assert "**2 conditions with phase-3 evidence** for this approved drug" in md
    assert "approved indication" not in md.split("phase-3 evidence")[0]
    qas = dict(intent_qas("drug", {"1": {"is_fda_approved": True, "max_phase": 4}, "4": b4}, "Imatinib"))
    assert qas["Is Imatinib an approved drug?"] == "Yes. Imatinib status: FDA-approved."
    assert "ChEMBL records no phase-4 indication" in qas["What is Imatinib used for?"]


def test_faq_fda_approved_drug_never_answered_not_approved():
    from atlas.page.faq_intents import intent_qas
    qas = dict(intent_qas("drug", {"1": {"is_fda_approved": True, "max_phase": 3},
                                   "4": {"indications": [{"name": "eczema", "max_phase": 3,
                                                          "approved": False, "late": True}]}},
                          "Clobetasol"))
    assert qas["Is Clobetasol an approved drug?"].startswith("Yes.")


def test_drug_glance_indications_tiers_and_targets_follow_mechanism():
    from atlas.page.drug_at_a_glance import at_a_glance
    md = at_a_glance({"1": {"max_phase": 4}, "2": {
        "mechanisms": [{"mechanism_of_action": "Sodium channel alpha subunit blocker"}],
        "mechanism_genes": [{"gene_symbol": "SCN1A"}],
        "primary_targets": [{"gene_symbol": "FZD8", "source": "gtopdb"}]},
        "4": {"indication_count": 9, "indications": [
            {"name": "epilepsy", "max_phase": 4, "approved": True},
            {"name": "pain", "max_phase": 3, "approved": False, "late": True}]}})
    assert "**Indications (ChEMBL):** 1 approved · 1 with phase-3 evidence" in md
    assert "SCN1A" in md and "FZD8" not in md


def _dz(name, orph_name=None, orph_inh=(), gencc=(), **b1):
    return {"1": {"canonical_name": name, "mondo_id": "MONDO:1", "orphanet_name": orph_name,
                  "orphanet_inheritance": list(orph_inh), **b1},
            "2": {}, "4": {"disease_gencc": [{"symbol": s, "classification": c, "moi": list(m)}
                                             for s, c, m in gencc]}, "5": {}}


def test_parent_orphanet_group_not_used_for_subtype():
    from atlas.disease.summary import orphanet_is_parent, definition
    from atlas.page.disease_declarative import inheritance_modes, declarative_sentence
    b = _dz("long QT syndrome 2", "Congenital long QT syndrome", ["Autosomal dominant", "Autosomal recessive"],
            [("KCNH2", "Definitive", ["Autosomal dominant"])],
            orphanet_definition="A rare group of genetic, cardiac rhythm diseases characterized by QT prolongation.",
            mondo_definition="A long QT syndrome caused by mutation in the KCNH2 gene.")
    assert orphanet_is_parent(b["1"])
    assert inheritance_modes(b) == (["Autosomal dominant"], "GenCC")
    assert definition(b["1"])[1] == "Mondo"
    s = declarative_sentence(b)
    assert "rare group" not in s and "autosomal recessive" not in s
    assert not orphanet_is_parent({"canonical_name": "Stickler syndrome", "orphanet_name": "Stickler syndrome"})
    assert orphanet_is_parent({"canonical_name": "Fanconi anemia complementation group C",
                               "orphanet_name": "Fanconi anemia"})


def test_orphanet_modes_narrowed_to_causal_gene():
    from atlas.page.disease_declarative import inheritance_modes
    b = _dz("cystinosis", "Cystinosis", ["Autosomal recessive", "X-linked recessive"],
            [("CTNS", "Definitive", ["Autosomal recessive"])])
    assert inheritance_modes(b) == (["Autosomal recessive"], "Orphanet, consistent with GenCC")


def test_clinvar_tier_needs_two_records_and_mendelian_anchor():
    from atlas.disease.cohort import causal_genes
    def b(plp, **b1):
        return {"1": {"canonical_name": "sialidosis", "orphanet_ids": ["Orphanet:1"], **b1},
                "2": {}, "4": {}, "3": {"clinvar_plp_by_gene": plp},
                "5": {"genes": [{"symbol": s, "evidence": {"clinvar": True}}
                                for s in ("NEU1", "SGCB", "MIR181A1HG")]}}
    assert causal_genes(b({"NEU1": 7, "SGCB": 1, "MIR181A1HG": 3})) == [("NEU1", "ClinVar-linked")]
    assert causal_genes(b({"NEU1": 7}, orphanet_ids=[])) == []          # no Mendelian anchor
    assert causal_genes(b({"NEU1": 7}, is_cancer=True)) == []


def test_faq_inherited_for_complex_disease_is_not_a_flat_yes():
    from atlas.page.faq_intents import intent_qas
    b = _dz("amyotrophic lateral sclerosis", "Amyotrophic lateral sclerosis",
            ["Autosomal dominant", "Autosomal recessive"])
    b["2"] = {"assoc_total": 169}
    a = dict(intent_qas("disease", b, "Amyotrophic lateral sclerosis"))["Is Amyotrophic lateral sclerosis inherited?"]
    assert a.startswith("Mostly not:") and "familial forms are inherited in an autosomal dominant or autosomal recessive pattern" in a


def test_subtype_genes_cancer_wording():
    from atlas.page.disease_declarative import _subtype_gene_clause
    b = {"1": {"child_count": 9, "is_cancer": True},
         "4": {"subtype_genes": [{"symbol": "EGFR"}], "subtype_gene_total": 1}}
    assert _subtype_gene_clause(b) == " whose hereditary forms are linked to EGFR"


def test_definition_text_cleanup():
    from atlas.disease.summary import clean_text
    assert clean_text("A disease that has_material_basis_in mutation of the &amp;lt;i&amp;gt;CFTR&amp;lt;/i&amp;gt; "
                      "gene and results_in fatigue , fatal .") == \
        "A disease that is caused by mutation of the CFTR gene and results in fatigue, fatal."


def test_article_agreement_in_definition_predicate():
    from atlas.disease.summary import definition_predicate
    assert definition_predicate("An disease or disorder caused by a helminth.", "x") == \
        "a disease or disorder caused by a helminth."
    assert definition_predicate("A inherited disorder of the eye.", "x") == "an inherited disorder of the eye."


def test_umbrella_subtype_approved_drugs():
    from atlas.disease.render import r_drugs_indicated as r_indicated_drugs
    from atlas.page.faq_intents import intent_qas
    from atlas.page.disease_declarative import _drugs_clause
    b = _dz("breast cancer")
    b["_indicated_drugs"] = [{"name": "ASPIRIN", "url": "/atlas/drug/aspirin/", "max_phase": 3,
                              "approved": False, "late": False}]
    b["_subtype_approved_drugs"] = [{"name": "PALBOCICLIB", "url": "/atlas/drug/palbociclib/",
                                     "subtypes": ["breast carcinoma"]}]
    md = r_indicated_drugs(b)
    assert "**1 drug approved for its subtypes**" in md and "breast carcinoma" in md
    a = dict(intent_qas("disease", b, "Breast cancer"))["What drugs are approved for Breast cancer?"]
    assert a.startswith("No drug is recorded as approved for Breast cancer as a whole")
    assert _drugs_clause({}, b["_indicated_drugs"], b["_subtype_approved_drugs"]) == \
        " Drugs approved for its subtypes include palbociclib."


def test_drug_pgx_adverse_events_recruiting_faq():
    from atlas.page.faq_intents import intent_qas
    b = {"1": {"is_fda_approved": True},
         "9": {"drug_labels": [{"source": "FDA", "level": "Actionable PGx", "genes": ["CYP2C19"]},
                               {"source": "EMA", "level": "Informative PGx", "genes": ["CYP2C19"]}],
               "guidelines": [{"source": "CPIC", "genes": "CYP2C19"}, {"source": "DPWG", "genes": "CYP2C19"}]},
         "13": {"total_reports": 1000, "most_reported": [
             {"reaction": "drug ineffective"}, {"reaction": "dyspnoea"}, {"reaction": "off label use"},
             {"reaction": "myocardial infarction"}]},
         "5": {"status_counts": {"RECRUITING": 3, "COMPLETED": 9},
               "top_trials": [{"id": "NCT1", "status": "RECRUITING"}]}}
    q = dict(intent_qas("drug", b, "Clopidogrel"))
    assert q["Is genetic testing recommended before taking Clopidogrel?"] == (
        "Pharmacogenomic guidance exists. Drug labels with pharmacogenomic content (PharmGKB): "
        "FDA — Actionable PGx (CYP2C19). Dosing guidelines: CPIC, DPWG (CYP2C19).")
    assert q["What adverse events are reported for Clopidogrel?"].endswith(": dyspnoea and myocardial infarction.")
    assert q["Are there recruiting clinical trials for Clopidogrel?"] == (
        "Yes — at least 3 registered trials are recruiting or about to recruit (ClinicalTrials.gov), e.g. NCT1.")


def test_gene_expression_disease_mouse_faq_and_glance():
    from atlas.page.faq_intents import intent_qas
    from atlas.page.at_a_glance import at_a_glance
    b = {"13": {"hpa": {"rna_tissue_specificity": "Tissue enriched", "top_tissues": ["Liver|198523.8", "Pancreas|1480"],
                        "subcellular_main": ["endoplasmic reticulum"]}},
         "12": {"gencc": [{"disease": "analbuminemia", "classification": "Definitive"},
                          {"disease": "x", "classification": "Limited"}],
                "mondo": [{"name": "familial dysalbuminemic hyperthyroxinemia"}]},
         "5": {"mouse_phenotypes": [{"statement": "decreased circulating serum albumin level"}],
               "mouse_phenotype_total": 4}}
    q = dict(intent_qas("gene", b, "ALB"))
    assert q["Where is ALB expressed?"] == ("Human Protein Atlas RNA tissue specificity: Tissue enriched; "
                                           "highest in liver and pancreas. Main subcellular location: endoplasmic reticulum.")
    assert "analbuminemia" in q["What diseases is ALB associated with?"] and " x" not in q["What diseases is ALB associated with?"]
    assert q["Is there a mouse model for ALB?"].startswith("Yes — MGI records 4 phenotypes")
    assert "**Expression (HPA):** Tissue enriched — highest in liver, pancreas" in at_a_glance(b)


def test_chemistry_brackets_not_links():
    from atlas.render_common import _cell
    assert _cell("{[3-(dimethylamino)propyl](methyl)amino}") == "{[3-(dimethylamino)propyl\\](methyl)amino}"
    assert _cell("[EGFR](/atlas/gene/EGFR/)") == "[EGFR](/atlas/gene/EGFR/)"


def test_anhydrous_dropped_from_display_names():
    from atlas.render_common import display_name
    assert display_name("DOCETAXEL ANHYDROUS") == "Docetaxel"
    assert display_name("Lisinopril anhydrous") == "Lisinopril"
    assert display_name("ABT199") == "ABT199"


def test_lastmod_from_content_hash():
    from atlas.lastmod import stamp, content_hash
    md = '---\ntitle: "X"\ngenerated_at: "2026-10-13T10:00:00+00:00"\natlas_version: "v1.13.0"\n---\nbody\n'
    md2 = md.replace("2026-10-13T10:00:00", "2026-11-01T09:00:00").replace("v1.13.0", "v1.14.0")
    assert content_hash(md) == content_hash(md2)                 # stamps ignored
    out, e = stamp(md2, "gene/X", {"gene/X": [content_hash(md), "2026-10-13"]}, "2026-11-01")
    assert 'lastmod: "2026-10-13"' in out and e[1] == "2026-10-13"
    out, e = stamp(md2.replace("body", "new body"), "gene/X", {"gene/X": [content_hash(md), "2026-10-13"]}, "2026-11-01")
    assert 'lastmod: "2026-11-01"' in out


def test_complex_disease_lead_inheritance_attaches_to_familial_forms():
    from atlas.page.disease_declarative import declarative_sentence
    b = _dz("amyotrophic lateral sclerosis", "Amyotrophic lateral sclerosis",
            ["Autosomal dominant", "Autosomal recessive"], [("FUS", "Definitive", [])])
    b["2"] = {"assoc_total": 169}
    s = declarative_sentence(b)
    assert ("whose monogenic forms are linked to FUS; its familial forms are inherited in an "
            "autosomal dominant or autosomal recessive pattern.") in s


def test_anhydrous_with_comma_separator():
    from atlas.render_common import display_name
    assert display_name("SODIUM PHOSPHATE, DIBASIC, ANHYDROUS") == "Sodium Phosphate, Dibasic"
