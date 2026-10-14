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
    assert "ChEMBL records its uses only up to phase 3" in qas["What is Imatinib used for?"]


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
