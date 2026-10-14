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
    return {"indications": [{"name": n, "approved": a, "max_phase": 4 if a else 3} for n, a in rows],
            "indication_count": len(rows)}


def test_drug_lead_separates_approved_from_trials():
    from atlas.page.drug_declarative import _indications_clause
    s = _indications_clause(_inds(("diabetes mellitus", True), ("obesity", False),
                                  ("pcos", False)))
    assert s == "; approved for 1 indication (diabetes mellitus), with clinical trials for 2 other conditions"
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
    b4 = {"indications": [{"name": "asthma", "approved": True, "max_phase": 4},
                          {"name": "copd", "approved": False, "max_phase": 3}]}
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
    b = _bundle("Harel-Tora syndrome", clinvar=("ANKLE2",), omim_ids=["MIM:1"])
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
    assert S.inheritance_phrase(["X-linked recessive"]) == "inherited in an X-linked recessive pattern"
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
        "No. Poziotinib status: Max clinical phase 3 (not approved).")
    assert "EGFR inhibitor" in qas["What is the mechanism of action of Poziotinib?"]
    gene = {"2": {"canonical_exon_count": 11, "canonical_transcript": "ENST00000269305",
                  "ensembl_transcript_count": 40}}
    assert dict(intent_qas("gene", gene, "TP53"))["How many exons does TP53 have?"] == (
        "The canonical Ensembl transcript of TP53 (ENST00000269305) has 11 exons. "
        "TP53 has 40 annotated Ensembl transcripts.")


# ── page hygiene ────────────────────────────────────────────────────────────

def test_prune_empty_headings_and_zero_labels():
    from atlas.pipeline import prune_empty
    md = "\n".join([
        "## Disease & clinical {#disease}", "",
        "### HPO {#hpo}", "",
        "### SIGNOR {#signor}", "",
        "#### Sub {#sub}", "",
        "### GWAS {#gwas-assoc}", "", "| a |", "",
        "**Orphanet (0):**", "**Mondo (2):** x, y", "",
        "### Linked {#linked}", "",
        "## Next {#next}", "", "see [here](#linked)"])
    out = prune_empty(md)
    assert "{#hpo}" not in out and "{#signor}" not in out and "{#sub}" not in out
    assert "{#disease}" in out and "{#gwas-assoc}" in out and "{#linked}" in out
    assert "Orphanet (0)" not in out and "**Mondo (2):** x, y" in out


def test_umbrella_lead_names_subtype_genes():
    from atlas.page.disease_declarative import declarative_sentence
    b = _bundle("inherited retinal dystrophy", children=40)
    b["4"]["subtype_genes"] = [{"symbol": s, "subtypes": n, "examples": []}
                               for s, n in (("PRPH2", 10), ("BEST1", 10), ("RHO", 10), ("RPE65", 8))]
    b["4"]["subtype_gene_total"] = 150
    s = declarative_sentence(b)
    assert "whose subtypes are caused by variants in PRPH2, BEST1, RHO, and 147 other genes" in s


def test_redirects_from_slug_history():
    from atlas.redirects import compute
    hist = {"disease": {"old-name": "MONDO:1", "same": "MONDO:2", "gone": "MONDO:3"},
            "gene": {"C1orf159": "HGNC:26062"}}
    built = [("disease", "new-name", "MONDO:1"), ("disease", "same", "MONDO:2"),
             ("gene", "C1orf159", "HGNC:26062")]
    rd, upd = compute(hist, built)
    assert rd == {"disease": {"old-name": "new-name"}}
    assert upd["disease"]["new-name"] == "MONDO:1" and upd["disease"]["old-name"] == "MONDO:1"
    # an old slug re-used by a live page is not redirected
    rd2, _ = compute({"disease": {"x": "MONDO:9"}}, [("disease", "x", "MONDO:8"), ("disease", "y", "MONDO:9")])
    assert rd2 == {}


def test_mondo_templated_definition_becomes_class():
    from atlas.page.disease_declarative import declarative_sentence
    b = _bundle("Focal segmental glomerulosclerosis 5",
                gencc=[("INF2", "Strong", ["Autosomal dominant"])],
                mondo_definition="Any focal segmental glomerulosclerosis in which the cause of the "
                                 "disease is a mutation in the INF2 gene.")
    s = declarative_sentence(b)
    assert s.startswith("**Focal segmental glomerulosclerosis 5** (MONDO:0000001) is a form of focal "
                        "segmental glomerulosclerosis caused by INF2 (GenCC Strong), inherited in an "
                        "autosomal dominant pattern.")
    assert "in which the cause" not in s


def test_drug_indications_render_matches_lead_counts():
    from atlas.drug.render import r_indications
    from atlas.page.drug_declarative import _indications_clause
    b4 = {"indication_count": 5, "indications": [
        {"name": "diabetes mellitus", "mondo_id": "MONDO:1", "max_phase": 4, "approved": True},
        {"name": "diabetes mellitus", "mondo_id": "MONDO:1", "max_phase": 4, "approved": True},
        {"name": "obesity", "mondo_id": "MONDO:2", "max_phase": 3, "approved": False},
        {"name": "EFO:0000001", "max_phase": 2, "approved": False},
        {"name": None, "max_phase": 4, "approved": True}]}
    md = r_indications(b4)
    assert "**1 approved indication** (ChEMBL phase 4)" in md and "**1 disease in clinical trials**" in md
    assert _indications_clause(b4) == ("; approved for 1 indication (diabetes mellitus), "
                                       "with clinical trials for 1 other condition")


# ── biobtree v2.16 data ─────────────────────────────────────────────────────

def test_rmc_render_and_glance():
    from atlas.gene.render import r_functional_genomics, rmc_summary
    from atlas.page.at_a_glance import at_a_glance
    rmc = {"transcript": "ENST1", "regions": [
        {"chrom": "2", "start": 100, "end": 199, "obs": 37, "exp": 154.9, "oe": 0.24, "p": 2.7e-21, "constrained": True},
        {"chrom": "2", "start": 200, "end": 299, "obs": 141, "exp": 185.4, "oe": 0.76, "p": None, "constrained": False}]}
    assert rmc_summary(rmc)[:2] == (2, 1)
    md = r_functional_genomics({"gnomad_rmc": rmc})
    assert "{#regional-constraint}" in md and "| 2:100–199 | 37 / 154.9 | 0.24 | yes (p 2.7e-21) |" in md
    g = at_a_glance({"3": {"gnomad_rmc": rmc}})
    assert "**Regional missense constraint (gnomAD):** 1 of 2 regions constrained (≈46% of expected missense variation)" in g
    one = {"transcript": "ENST2", "regions": [{"chrom": "17", "start": 1, "end": 9, "oe": 0.885, "constrained": False}]}
    assert "no regional variation along `ENST2` — one region, missense o/e 0.89." in r_functional_genomics({"gnomad_rmc": one})


def test_gtex_egene_render_and_glance():
    from atlas.gene.render import r_expression
    from atlas.page.at_a_glance import at_a_glance
    g = {"tested": 95, "egenes": [
        {"tissue": "Skin Sun Exposed Lower leg", "type": "sqtl", "q": 6.1e-7, "slope": "0.581", "lead_variant": "17:7673183:G:A"},
        {"tissue": "Lung", "type": "eqtl", "q": 0.004, "slope": "-1.17", "lead_variant": "17:1:A:C"}]}
    md = r_expression({"symbol": "TP53", "gtex_egene": g})
    assert "(eQTL) in **1** tissue and its splicing (sQTL) in **1**" in md
    assert "| Skin Sun Exposed Lower leg | sqtl | 6.1e-07 | +0.58 | 17:7673183:G:A |" in md
    assert "**GTEx eGene:** eQTL in 1 tissue, sQTL in 1 (q < 0.05)" in at_a_glance({"11": {"gtex_egene": g}})


def test_transcripts_mane_tsl():
    from atlas.gene.render import r_transcripts
    md = r_transcripts({"ensembl_transcript_count": 3, "ensembl_transcripts": [
        {"id": "ENST2", "biotype": "protein_coding", "tsl": "2"},
        {"id": "ENST1", "biotype": "protein_coding", "mane_select": True, "tsl": "1"},
        {"id": "ENST3", "biotype": "protein_coding", "mane_plus_clinical": True, "tsl": "NA"}]})
    assert "**MANE Select:** `ENST1` · **MANE Plus Clinical:** `ENST3`" in md
    assert "TSL 1: 1, TSL 2: 1" in md
    assert "`ENST1`, `ENST3`, `ENST2`" in md


def test_disease_jsonld_wikidata_sameas():
    from atlas.page.disease_jsonld import same_as_urls
    assert "https://www.wikidata.org/wiki/Q11081" in same_as_urls(
        {"1": {"mondo_id": "MONDO:0004975", "wikidata_qid": "Q11081"}})


def test_cp1252_c1_chars_do_not_break_frontmatter():
    import yaml
    from atlas.pipeline import _yaml_escape
    v = _yaml_escape("Shprintzen\x96Goldberg omphalocele syndrome \x85")
    assert v == "Shprintzen–Goldberg omphalocele syndrome …"
    assert yaml.safe_load(f'd: "{v}"')["d"] == v


def test_drug_lead_unknown_mechanism_and_phase3_with_phase4_rows():
    from atlas.page.drug_declarative import declarative_sentence, _targets_clause
    assert _targets_clause({"mechanisms": [{"mechanism_of_action": "Unknown"}]}) == ""
    b = {"1": {"canonical_name": "ACETYLCARNITINE", "chembl_id": "CHEMBL1", "max_phase": 3,
               "molecule_type": "Small molecule"},
         "2": {"mechanisms": [{"mechanism_of_action": "Unknown"}]},
         "4": {"indications": [{"name": "adhd", "approved": True, "max_phase": 4}]}}
    s = declarative_sentence(b)
    assert "unknown" not in s.lower()
    assert "phase-3 clinical-stage" in s and "approved for" not in s
    assert "with phase-4 (approval-stage) ChEMBL records for 1 indication (adhd)" in s
