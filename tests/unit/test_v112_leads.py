"""v1.12 lead / claim corrections: IntAct own-gene rows, drug approved-vs-trial
indications, ClinVar-only causal claims, disease definition-first leads,
inheritance fallback gating, intent FAQ questions."""


# ── IntAct ──────────────────────────────────────────────────────────────────

def test_intact_keeps_only_this_genes_pairs_one_row_per_partner():
    from atlas.gene.sections.s08_interactions import _intact_partners
    rows = [
        {"protein_a_gene": "CDK4", "protein_b_gene": "CCND3",          # foreign pair
         "interaction_type": 'psi-mi:"MI:0914"(association)', "confidence_score": "0.98"},
        {"protein_a_gene": "TMEM126A", "protein_b_gene": "IPO7",
         "interaction_type": 'psi-mi:"MI:0914"(association)', "confidence_score": "0.35"},
        {"protein_a_gene": "IPO7", "protein_b_gene": "TMEM126A",       # same partner again
         "interaction_type": 'psi-mi:"MI:0407"(direct interaction)', "confidence_score": "0.60"},
        {"protein_a_gene": "Tuba3a", "protein_b_gene": "TMEM126A",     # mouse symbol
         "interaction_type": 'psi-mi:"MI:0914"(association)', "confidence_score": "0.90"},
        {"protein_a_gene": "TMEM126A", "protein_b_gene": "TMEM126A",   # self
         "interaction_type": 'psi-mi:"MI:0914"(association)', "confidence_score": "0.99"},
        {"protein_a_gene": "TMEM126A", "protein_b_gene": "C1orf159",
         "interaction_type": 'psi-mi:"MI:0915"(physical association)', "confidence_score": "0.40"},
    ]
    out = _intact_partners(rows, "TMEM126A")
    assert [p["partner"] for p in out] == ["IPO7", "C1orf159"]
    ipo7 = out[0]
    assert ipo7["records"] == 2 and ipo7["score"] == "0.60" and ipo7["type"] == "direct interaction"


def test_intact_render_partner_table():
    from atlas.gene.render import r_interactions
    md = r_interactions({"symbol": "TMEM126A", "intact_count": 87, "intact_sampled": 87,
                         "intact_partner_count": 1,
                         "intact": [{"partner": "IPO7", "type": "association",
                                     "score": "0.350", "records": 1}]})
    assert "| Partner | Interaction type | MIscore | Records |" in md
    assert "| IPO7 | association | 0.350 | 1 |" in md
    assert "CDK4" not in md


# ── drug lead / JSON-LD ─────────────────────────────────────────────────────

def _inds(*rows):
    return {"indications": [{"name": n, "approved": a} for n, a in rows],
            "indication_count": len(rows)}


def test_drug_lead_separates_approved_from_trials():
    from atlas.page.drug_declarative import _indications_clause
    s = _indications_clause(_inds(("diabetes mellitus", True), ("obesity", False),
                                  ("pcos", False)))
    assert s == "; approved for 1 indication (diabetes mellitus), and in clinical trials for 2 other conditions"
    s = _indications_clause(_inds(("non-small cell lung carcinoma", False), ("breast cancer", False)))
    assert s.startswith("; investigated in clinical trials for 2 conditions including")
    assert "approved" not in s and "indicated" not in s


def test_drug_lead_names_the_mechanism():
    from atlas.page.drug_declarative import _targets_clause
    b2 = {"mechanisms": [{"mechanism_of_action": "Vascular endothelial growth factor receptor inhibitor"},
                         {"mechanism_of_action": "HMG-CoA reductase inhibitor"},
                         {"mechanism_of_action": "Tubulin inhibitor"}],
          "mechanism_genes": [{"gene_symbol": "KDR"}]}
    s = _targets_clause(b2)
    assert s == (" that acts as a vascular endothelial growth factor receptor inhibitor and "
                 "HMG-CoA reductase inhibitor (+1 more curated mechanism), targeting KDR")


def test_drug_jsonld_treats_only_approved():
    from atlas.page.drug_jsonld import _treats
    b4 = {"indications": [{"name": "asthma", "approved": True},
                          {"name": "copd", "approved": False}]}
    assert [t["name"] for t in _treats(b4)] == ["asthma"]


# ── disease causal claims ───────────────────────────────────────────────────

def _bundle(name, gwas=0, children=0, gencc=(), clinvar=(), orph_inh=(), **b1):
    genes = [{"symbol": s, "evidence": {"clinvar": True}} for s in clinvar]
    return {"1": {"canonical_name": name, "mondo_id": "MONDO:0000001",
                  "child_count": children, "orphanet_inheritance": list(orph_inh), **b1},
            "2": {"assoc_total": gwas},
            "4": {"disease_gencc": [{"symbol": s, "classification": c, "moi": list(m)}
                                    for s, c, m in gencc]},
            "5": {"genes": genes, "gene_count": len(genes)}}


def test_chromosomal_syndrome_gets_no_clinvar_causal_gene():
    from atlas.disease.cohort import causal_genes
    from atlas.page.disease_declarative import declarative_sentence
    b = _bundle("2q13 microdeletion syndrome", clinvar=("RGPD8",))
    assert causal_genes(b) == []
    assert "caused by" not in declarative_sentence(b)


def test_common_disease_gets_no_clinvar_monogenic_claim():
    from atlas.page.disease_declarative import declarative_sentence
    b = _bundle("asthma", gwas=3000, clinvar=("CARD11", "TBXA2R"))
    s = declarative_sentence(b)
    assert "CARD11" not in s and "monogenic" not in s


def test_clinvar_only_rare_disease_reports_evidence_not_cause():
    from atlas.page.disease_declarative import declarative_sentence, causal_kind
    b = _bundle("Harel-Tora syndrome", clinvar=("ANKLE2",))
    assert causal_kind(b)[0] == "clinvar"
    s = declarative_sentence(b)
    assert "with pathogenic variants reported in ANKLE2 (ClinVar)" in s
    assert "caused by" not in s


def test_gencc_inheritance_only_for_a_caused_disease():
    from atlas.page.disease_declarative import inheritance_modes
    caused = _bundle("Some syndrome", gencc=[("ABC1", "Definitive", ["Autosomal recessive"])])
    assert inheritance_modes(caused) == (["Autosomal recessive"], "GenCC")
    t2d = _bundle("type 2 diabetes mellitus", gwas=8000,
                  gencc=[("GCK", "Strong", ["Autosomal dominant"])])
    assert inheritance_modes(t2d) == ([], "")
    orph = _bundle("cystic fibrosis", orph_inh=["Autosomal recessive"])
    assert inheritance_modes(orph) == (["Autosomal recessive"], "Orphanet")


# ── disease definition-first lead ───────────────────────────────────────────

def test_disease_lead_definition_inheritance_features():
    from atlas.page.disease_declarative import declarative_sentence
    b = _bundle("cystic fibrosis", gencc=[("CFTR", "Definitive", ["Autosomal recessive"])],
                orph_inh=["Autosomal recessive"],
                orphanet_definition="A rare, genetic pulmonary disorder characterized by thick mucus.",
                phenotypes=[{"hpo_term": "Bronchiectasis", "frequency_value": 0.9},
                            {"hpo_term": "Malabsorption", "frequency_value": 0.95},
                            {"hpo_term": "Abnormal EEG", "frequency_value": 0.1}])
    s = declarative_sentence(b)
    assert s.startswith("**Cystic fibrosis** (MONDO:0000001) is a disease caused by CFTR "
                        "(GenCC Definitive), inherited in an autosomal recessive pattern. "
                        "It is a rare, genetic pulmonary disorder characterized by thick mucus (Orphanet).")
    assert "include malabsorption, bronchiectasis, and Abnormal EEG." in s
    assert "cohort gene" not in s


def test_disease_lead_definition_as_predicate():
    from atlas.page.disease_declarative import declarative_sentence
    b = _bundle("asthma", gwas=3000,
                mondo_definition="A bronchial disease that is characterized by chronic inflammation.")
    assert declarative_sentence(b).startswith(
        "**Asthma** (MONDO:0000001) is a bronchial disease that is characterized by chronic "
        "inflammation (Mondo).")


def test_summary_helpers():
    from atlas.disease import summary as S
    assert S.inheritance_phrase(["Autosomal dominant", "Autosomal recessive"]) == \
        "inherited in an autosomal dominant or autosomal recessive pattern"
    assert S.inheritance_phrase(["Mitochondrial inheritance"]) == "with mitochondrial inheritance"
    assert S._uncap_mesh("A degenerative disease of the BRAIN marked by DEMENTIA in HIV and AIDS.") == \
        "A degenerative disease of the brain marked by dementia in HIV and AIDS."
    assert S._uncap_mesh("caused by the CYSTIC FIBROSIS TRANSMEMBRANE CONDUCTANCE REGULATOR") == \
        "caused by the cystic fibrosis transmembrane conductance regulator"
    assert S.definition_predicate("Asthma is a chronic disease.", "asthma") == "a chronic disease."
    assert S.definition_predicate("Persistently high blood pressure.", "hypertension") == ""


def test_meta_description_drops_accession():
    from atlas.page.meta_facts import _meta_description
    b = _bundle("asthma", mondo_definition="A bronchial disease of the airways.")
    assert _meta_description("disease", b).startswith("Asthma is a bronchial disease")


def test_capped_cohort_count_is_not_ranked():
    from atlas.page.disease_at_a_glance import at_a_glance
    b = _bundle("asthma")
    b["5"]["gene_count"] = 75
    md = at_a_glance(b)
    assert "**Cohort genes:** 75 (strongest-evidence selection; the list is capped)" in md
    assert "top 1%" not in md


# ── intent FAQ ──────────────────────────────────────────────────────────────

def test_faq_polygenic_disease_not_caused_or_inherited():
    from atlas.page.faq_intents import intent_qas
    t2d = _bundle("type 2 diabetes mellitus", gwas=8000,
                  gencc=[("GCK", "Strong", ["Autosomal dominant"])])
    qas = dict(intent_qas("disease", t2d, "Type 2 diabetes mellitus"))
    assert "monogenic forms are linked to GCK" in qas["What causes Type 2 diabetes mellitus?"]
    assert "Is Type 2 diabetes mellitus inherited?" not in qas


def test_faq_drug_unapproved_and_gene_exons():
    from atlas.page.faq_intents import intent_qas
    drug = {"1": {"max_phase": 3}, "2": {"mechanisms": [{"mechanism_of_action": "EGFR inhibitor"}]},
            "4": _inds(("lung carcinoma", False))}
    qas = dict(intent_qas("drug", drug, "Poziotinib"))
    assert qas["Is Poziotinib an approved drug?"].startswith(
        "No approved indication is recorded for Poziotinib in ChEMBL; its highest development phase is 3.")
    assert "EGFR inhibitor" in qas["What is the mechanism of action of Poziotinib?"]
    gene = {"2": {"canonical_exon_count": 11, "canonical_transcript": "ENST00000269305",
                  "ensembl_transcript_count": 40}}
    assert dict(intent_qas("gene", gene, "TP53"))["How many exons does TP53 have?"] == (
        "The canonical Ensembl transcript of TP53 (ENST00000269305) has 11 exons. "
        "TP53 has 40 annotated Ensembl transcripts.")
