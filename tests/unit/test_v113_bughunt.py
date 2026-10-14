"""v1.13 bug-hunt fixes: salt-form indications, one inheritance scope for lead
and FAQ, approved-only drug answers, GenCC-refuted links, wording."""


def _dz(name, orph_inh=(), gencc=(), gwas=0, **b1):
    return {"1": {"canonical_name": name, "mondo_id": "MONDO:1", "orphanet_name": name,
                  "orphanet_inheritance": list(orph_inh), **b1},
            "2": {"assoc_total": gwas},
            "4": {"disease_gencc": [{"symbol": s, "classification": c, "moi": ["Autosomal dominant"]}
                                    for s, c in gencc]}, "5": {}}


def test_salt_form_indications_merged_with_max_phase(monkeypatch):
    from atlas.drug import anchors as A
    monkeypatch.setattr(A, "map_all", lambda i, chain, cap=None:
                        [{"id": "MONDO:9", "name": "chronic myeloid leukemia"}] if i == "EFO:1" else [])
    mol = {"indications": [{"efo": "EFO:1", "highestDevelopmentPhase": 3}]}
    salt = [{"efo": "EFO:1", "highestDevelopmentPhase": 4}]
    (rec,) = A._resolve_indications(mol, salt)
    assert rec.max_phase == 4 and rec.mondo_id == "MONDO:9"


def test_indication_without_mondo_keeps_mesh_label(monkeypatch):
    from atlas.drug import anchors as A
    monkeypatch.setattr(A, "map_all", lambda *a, **k: [])
    monkeypatch.setattr(A, "entry", lambda i, src: {"Attributes": {"Mesh": {"descriptor_name": "Pain"}}}
                        if src == "mesh" else {"Attributes": {"Ontology": {"name": "obsolete_pain"}}})
    (rec,) = A._resolve_indications({"indications": [{"efo": "EFO:3", "mesh": "D010146",
                                                       "highestDevelopmentPhase": 4}]})
    assert rec.name == "Pain" and rec.slug is None and rec.mondo_id is None


def test_hereditary_syndrome_faq_matches_lead():
    from atlas.page.faq_intents import intent_qas
    from atlas.page.disease_declarative import declarative_sentence, inheritance_scope
    b = _dz("Hereditary breast ovarian cancer syndrome", ["Autosomal dominant"], is_cancer=True)
    assert inheritance_scope(b) == "disease"
    a = dict(intent_qas("disease", b, "HBOC"))["Is HBOC inherited?"]
    assert a.startswith("Yes.") and "autosomal dominant" in a
    assert "inherited in an autosomal dominant pattern" in declarative_sentence(b)


def test_umbrella_hereditary_cancer_syndrome_is_subtypes_not_sporadic():
    from atlas.page.faq_intents import intent_qas
    from atlas.page.disease_declarative import declarative_sentence, inheritance_scope
    b = _dz("Lynch syndrome", ["Autosomal dominant"], [("MLH1", "Definitive")],
            is_cancer=True, child_count=5)
    assert inheritance_scope(b) == "subtypes"
    a = dict(intent_qas("disease", b, "Lynch syndrome"))["Is Lynch syndrome inherited?"]
    assert "Mostly not" not in a and "its subtypes are inherited in an autosomal dominant" in a
    s = declarative_sentence(b)
    assert "hereditary forms" not in s and "its subtypes are inherited" in s


def test_sporadic_cancer_familial_forms_have_multifactorial():
    from atlas.page.faq_intents import intent_qas
    from atlas.page.disease_declarative import declarative_sentence
    b = _dz("glioblastoma", ["Multigenic/multifactorial"], is_cancer=True)
    a = dict(intent_qas("disease", b, "Glioblastoma"))["Is Glioblastoma inherited?"]
    assert a.startswith("Mostly not:") and "familial forms have multigenic/multifactorial inheritance" in a
    s = declarative_sentence(b)
    assert "are with" not in s and "whose familial forms have multigenic/multifactorial inheritance" in s


def test_disease_faq_approved_only_phase4_late_separate():
    from atlas.page.faq_intents import intent_qas
    b = _dz("idiopathic pulmonary fibrosis")
    b["_indicated_drugs"] = [{"name": "PIRFENIDONE", "approved": True, "late": False},
                             {"name": "AMBRISENTAN", "approved": False, "late": True}]
    q = dict(intent_qas("disease", b, "IPF"))
    assert "Ambrisentan" not in q["What drugs are approved for IPF?"]
    assert "Ambrisentan" in q["Which approved drugs have been tested in phase 3 for IPF?"]
    b["_indicated_drugs"] = b["_indicated_drugs"][1:]
    assert "What drugs are approved for IPF?" not in dict(intent_qas("disease", b, "IPF"))


def test_disease_jsonld_drug_is_phase4_only():
    from atlas.page.disease_jsonld import _drugs
    rows = [{"name": "AMBRISENTAN", "url": "/atlas/drug/ambrisentan/", "approved": False, "late": True}]
    assert not _drugs({"_indicated_drugs": rows})


def test_phase4_drug_not_fda_is_still_approved():
    from atlas.page.faq_intents import intent_qas
    b = {"1": {"is_fda_approved": False, "max_phase": 4},
         "4": {"indications": [{"name": "bacterial infectious disease", "max_phase": 4, "approved": True}]}}
    a = dict(intent_qas("drug", b, "Cefamandole"))["Is Cefamandole an approved drug?"]
    assert a.startswith("Yes") and "no FDA approval" in a


def test_gene_faq_ignores_refuted_gencc_links():
    from atlas.page.faq_intents import intent_qas
    b = {"12": {"gencc": [{"disease": "hypertrophic cardiomyopathy 3", "classification": "Definitive"}],
                "mondo": [{"id": "MONDO:1", "name": "arrhythmogenic right ventricular cardiomyopathy",
                           "gencc": "Limited", "clinvar_records": 0},
                          {"id": "MONDO:2", "name": "dilated cardiomyopathy 1Y", "gencc": None,
                           "clinvar_records": 4}]}}
    a = dict(intent_qas("gene", b, "TPM1"))["What diseases is TPM1 associated with?"]
    assert "arrhythmogenic" not in a and "dilated cardiomyopathy 1Y" in a


def test_glance_causal_label_follows_causal_kind():
    from atlas.page.disease_at_a_glance import at_a_glance
    b = _dz("type 2 diabetes mellitus", gencc=[("GCK", "Strong")], gwas=500)
    md = at_a_glance(b)
    assert "**Gene of its monogenic forms:** GCK" in md and "Causal gene" not in md


def test_article_only_fixes_clear_errors():
    from atlas.disease.summary import _article
    assert _article("a urea cycle disorder") == "a urea cycle disorder"
    assert _article("an X-linked disorder") == "an X-linked disorder"
    assert _article("an RNA virus") == "an RNA virus"
    assert _article("an disease") == "a disease"
    assert _article("a inherited disorder") == "an inherited disorder"


def test_noncoding_filter_spares_protein_coding_hg_genes():
    from atlas.disease.cohort import _NONCODING_NEIGHBOR as R
    assert not R.search("PLEKHG5") and not R.search("PGGHG")
    assert R.search("MIR17HG") and R.search("CFTR-AS1") and R.search("SNHG5")


def test_fake_links_escaped_in_prose_not_scripts():
    from atlas.pipeline import escape_fake_links
    md = 'raises [Ca(2+)](i) via [MLH1](/atlas/gene/MLH1/)\n<script>{"a":"x](y)"}</script>'
    out = escape_fake_links(md)
    assert "[Ca(2+)\\](i)" in out and "[MLH1](/atlas/gene/MLH1/)" in out and '"x](y)"' in out
    assert escape_fake_links(out) == out                       # idempotent


def test_capped_table_singular_caption():
    from atlas.render_common import capped_table
    assert capped_table(["Drug"], [("x",)], 10, noun="drugs via a complex or family").startswith(
        "1 drug via a complex or family:")


def test_cancer_flag_ignores_definition_text(monkeypatch):
    from atlas.disease import anchors as A
    monkeypatch.setattr(A, "_under_neoplasm", lambda m, **k: False)
    en = {"Attributes": {"Ontology": {"definition": "Exposure increases the risk of cancer."}}}
    assert not A._is_cancer("MONDO:1", en, "heavy metal poisoning")
    assert A._is_cancer("MONDO:1", en, "Wilms tumor 1")
