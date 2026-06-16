"""Gene page canonical H2 taxonomy + the functional-residue map
(docs/PAGE_CONTRACT.md, docs/MOLECULAR_ENRICHMENT.md)."""
import atlas.pipeline as P
from atlas.gene import render as R

# The FROZEN gene H2 set (order matters) — Summary/#related are added by
# assemble_page/related_block, so render_all emits these six.
GENE_H2 = [
    "## Identifiers {#identifiers}",
    "## Gene structure {#gene-structure}",
    "## Protein {#protein}",
    "## Function {#function}",
    "## Disease & clinical {#disease}",
    "## Drugs & pharmacology {#drugs}",
]


# ── zoning ──────────────────────────────────────────────────────────────────
def _stub_renderers(monkeypatch):
    monkeypatch.setattr(R, "RENDER",
                        {s: (lambda b, _s=s: f"## Section {_s}\n\nbody{_s}")
                         for s in (str(i) for i in range(1, 13))})
    monkeypatch.setattr(R, "r_functional_genomics", lambda b: "## Functional genomics\n\nfg")
    monkeypatch.setattr(R, "r_generifs", lambda b: "## GeneRIF\n\ngr")
    monkeypatch.setattr(R, "r_residue_map", lambda b: "## Functional residue map\n\nrm")


def _h2(md):
    return [l for l in md.splitlines() if l.startswith("## ")]


def test_canonical_h2_set_and_order_coding(monkeypatch):
    _stub_renderers(monkeypatch)
    bundle = {s: {} for s in (str(i) for i in range(1, 13))}
    bundle["3"] = {"reviewed_uniprot": ["P1"], "canonical_uniprot": "P1"}
    md = P.render_all(bundle)
    assert _h2(md) == GENE_H2                      # exact set + frozen order
    assert '<a id="protein-P1"></a>' in md         # JSON-LD @id linkage kept
    # sub-sections demoted to H3 under their canonical H2
    assert "### Section 1" in md and "\n## Section 1" not in md
    prot = md.split("{#protein}")[1].split("{#function}")[0]
    assert "### Functional residue map" in prot
    gs = md.split("{#gene-structure}")[1].split("{#protein}")[0]
    assert "### Functional genomics" in gs
    assert "### Section 10" in md.split("{#drugs}")[1]   # drugs


def test_ncrna_emits_all_sections_with_placeholders(monkeypatch):
    """Emit-even-if-empty: the TOC is identical across genes; ncRNA shows the
    same H2 set, with informative placeholders for the protein-layer sections."""
    _stub_renderers(monkeypatch)
    bundle = {s: {} for s in (str(i) for i in range(1, 13))}
    bundle["3"] = {"reviewed_uniprot": [], "canonical_uniprot": None}
    bundle["_noncoding"] = "lncRNA"
    md = P.render_all(bundle)
    assert _h2(md) == GENE_H2                      # same set, not collapsed
    assert "*Non-coding RNA — no protein product" in md.split("{#protein}")[1]
    assert 'id="protein-' not in md                # no JSON-LD @id anchor for ncRNA


DISEASE_H2 = [
    "## Clinical features {#clinical}",
    "## Identifiers {#identifiers}",
    "## Disease family {#family}",
    "## Genetics & variants {#genetics}",
    "## Genes & proteins {#genes}",
    "## Function {#function}",
    "## Therapeutics {#drugs}",
    "## Clinical trials & evidence {#trials}",
]
DRUG_H2 = [
    "## Identifiers {#identifiers}",
    "## Targets {#targets}",
    "## Indications & clinical {#indications}",
    "## Pharmacology {#pharmacology}",
    "## Related molecules {#related-molecules}",
]


def test_clinical_description_dual_source():
    """Clinical description shows Orphanet + MeSH source-labelled (both when both
    exist, Orphanet first) plus an inheritance/onset line; '' when nothing."""
    from atlas.disease import render as DR
    md = DR.r_clinical_description({"mesh_scope_note": "Persistently high BLOOD PRESSURE."})
    assert "## Clinical description" in md
    assert "**MeSH:** Persistently high BLOOD PRESSURE." in md
    assert "**Orphanet:**" not in md
    md2 = DR.r_clinical_description({
        "orphanet_definition": "Brachytelephalangic chondrodysplasia punctata is …",
        "mesh_scope_note": "A heterogeneous disorder …",
        "orphanet_inheritance": ["X-linked recessive"],
        "orphanet_onset": ["Antenatal", "Neonatal"]})
    assert md2.index("**Orphanet:**") < md2.index("**MeSH:**")     # Orphanet leads
    assert "**Inheritance:** X-linked recessive" in md2
    assert "**Onset:** Antenatal, Neonatal" in md2
    assert DR.r_clinical_description({}) == ""


def test_disease_canonical_h2_set_and_order():
    """Every canonical section always emits, in the frozen order
    (docs/PAGE_CONTRACT.md) — sub-renderers carry their own empty-state text, so
    the H2 set is identical across every disease page even with no data."""
    from atlas.disease import render as DR
    md = DR.render_all({str(i): {} for i in range(1, 15)})
    assert _h2(md) == DISEASE_H2


def test_drug_canonical_h2_set_and_order():
    from atlas.drug import render as DRR
    md = DRR.render_all({str(i): {} for i in range(1, 13)})
    assert _h2(md) == DRUG_H2


def test_demote_bumps_headings():
    assert P._demote("## A\ntext\n### B") == "### A\ntext\n#### B"
    assert P._demote("no heading") == "no heading"


# ── residue map ───────────────────────────────────────────────────────────────
def _feat(u, t, b, e=None, d=None):
    return {"uniprot": u, "type": t, "begin": str(b),
            "end": str(e if e is not None else b), "description": d}


def test_residue_map_categorizes_and_orders():
    b = {"reviewed_uniprot": ["P0"], "canonical_uniprot": "P0", "ufeatures": [
        _feat("P0", "active site", 837, d="proton acceptor"),
        _feat("P0", "binding site", 718, 726),
        _feat("P0", "modified residue", 693),
        _feat("P0", "disulfide bond", 31, 58),
        _feat("P0", "glycosylation site", 56),
        _feat("P0", "mutagenesis site", 275, d="reduced activity"),
        _feat("P0", "strand", 1, 5),          # structural — excluded
    ]}
    md = R.r_residue_map(b)
    assert "## Functional residue map" in md
    # each residue category is now its own H4 sub-section (count kept in heading),
    # with a residue table underneath (web renderer wraps each into a collapsible)
    assert "### Catalytic / active sites (1)" in md
    assert "| 837 | proton acceptor |" in md    # perres: Position | Role table
    assert "### Ligand- & substrate-binding residues (1)" in md
    assert "| 718–726 |" in md
    assert "### Post-translational modifications (1)" in md
    assert "### Disulfide bonds (1)" in md
    assert "### Mutagenesis-validated functional residues (1)" in md
    assert "strand" not in md                  # secondary structure not surfaced


def test_residue_map_labels_modifications_by_type():
    """PTM / glycosylation positions are grouped by modification identity (from
    the UniProt feature description) into a `Type | Positions` table instead of a
    bare, meaningless number list. Unlabeled features (most intrachain disulfide
    bonds) collapse onto one row tagged with the group's singular noun."""
    b = {"reviewed_uniprot": ["P0"], "canonical_uniprot": "P0", "ufeatures": [
        _feat("P0", "modified residue", 229, d="phosphoserine"),
        _feat("P0", "modified residue", 695, d="phosphoserine"),
        _feat("P0", "modified residue", 869, d="phosphotyrosine; by src"),
        _feat("P0", "glycosylation site", 56, d="n-linked (glcnac...) asparagine"),
        _feat("P0", "disulfide bond", 31, 58),     # no description
    ]}
    md = R.r_residue_map(b)
    assert "| Phosphoserine | 229, 695 |" in md    # grouped, ordered, no kinase tail
    assert "| Phosphotyrosine | 869 |" in md        # '; by src' tail dropped
    assert "| N-linked (glcnac...) asparagine | 56 |" in md
    assert "### Disulfide bonds (1)" in md
    # the disulfide bond carries no description → falls onto a 'Disulfide bond' row
    assert "| Disulfide bond | 31–58 |" in md


def test_residue_map_binding_site_shows_ligand():
    """Binding-site rows show the ligand (ATP, Mg(2+)…) from the clean `ligand`
    field, not the lowercased description; active sites keep their description."""
    b = {"reviewed_uniprot": ["P0"], "canonical_uniprot": "P0", "ufeatures": [
        {"uniprot": "P0", "type": "binding site", "begin": "718", "end": "726",
         "description": "atp", "ligand": "ATP"},
        {"uniprot": "P0", "type": "active site", "begin": "837", "end": "837",
         "description": "proton acceptor", "ligand": None},
    ]}
    md = R.r_residue_map(b)
    assert "| 718–726 | ATP |" in md            # clean ligand casing, not 'atp'
    assert "| 837 | proton acceptor |" in md     # active site falls back to description


def test_hpa_protein_shows_antibody_reliability():
    """HPA protein block surfaces the antibody-staining reliability tier (IH/IF),
    omitting whichever field is absent."""
    md = R.r_hpa_protein({"13": {"hpa": {
        "subcellular_main": ["Nucleoplasm"],
        "reliability_ih": "enhanced", "reliability_if": "supported"}}})
    assert "**Antibody reliability (HPA):** IH enhanced · IF supported" in md
    # IF absent → only IH shown, no dangling separator
    md2 = R.r_hpa_protein({"13": {"hpa": {
        "subcellular_main": ["Cytosol"], "reliability_ih": "enhanced"}}})
    assert "**Antibody reliability (HPA):** IH enhanced" in md2
    assert "·" not in md2.split("Antibody reliability")[1]


def test_residue_map_per_product_skips_empty():
    b = {"reviewed_uniprot": ["P0", "P1", "P2"], "canonical_uniprot": "P0",
         "ufeatures": [
            _feat("P0", "active site", 10),
            _feat("P1", "binding site", 20, 25),
            _feat("P2", "chain", 1, 99),        # only a non-residue feature
         ]}
    md = R.r_residue_map(b)
    # two products with residues → each category H4 carries a product suffix
    # (no separate per-product bold header — keeps the heading hierarchy flat)
    assert "### Catalytic / active sites (1) — P0 (canonical)" in md
    assert "### Ligand- & substrate-binding residues (1) — P1" in md
    assert "P2" not in md                        # residue-less product gets no heading


def test_residue_map_single_product_no_header():
    """One product with residues → no redundant per-product header."""
    b = {"reviewed_uniprot": ["P0", "P9"], "canonical_uniprot": "P0", "ufeatures": [
        _feat("P0", "active site", 10),
        _feat("P9", "chain", 1, 50),
    ]}
    md = R.r_residue_map(b)
    assert "**P0 (canonical)**" not in md
    assert "Catalytic / active sites (1)" in md


def test_residue_map_empty_when_no_features():
    assert R.r_residue_map({"reviewed_uniprot": ["P0"], "ufeatures": []}) == ""


def test_sense_gene_orientation_bullet(monkeypatch):
    """An antisense lncRNA orients to its sense gene, quoting the SENSE gene's
    biology (attributed), manifest-gated."""
    from atlas.page import at_a_glance as AAG
    from atlas.page import evidence, links
    monkeypatch.setattr(links, "gene_url",
                        lambda symbol=None, **k: f"/atlas/gene/{symbol}/" if symbol == "FER" else None)
    evidence._COMP = {"gene": {"FER": {"gwas_count": 33, "drug_count": 56}}}
    b = AAG._sense_gene_bullet("FER-AS1")
    assert "antisense to" in b and "[FER](/atlas/gene/FER/)" in b
    assert "33 GWAS associations" in b and "56 ChEMBL molecules" in b
    assert "not this transcript" in b                       # honest attribution
    assert AAG._sense_gene_bullet("TP53") == ""             # not antisense
    assert AAG._sense_gene_bullet("XXX-AS1") == ""          # sense gene not built
    evidence.reset()


def test_parent_evidence_clause(monkeypatch):
    """A thin subtype quantifies its parent's evidence from the parent's frozen
    components — attributed, never folded into this page's score."""
    from atlas.disease import render as DR
    from atlas.page import evidence, links
    monkeypatch.setattr(links, "_lookup", lambda et, *a: "heavy-metal-poisoning")
    evidence._COMP = {"disease": {"heavy-metal-poisoning":
                                  {"trial_count": 4, "gwas_count": 1, "drug_count": 1}}}
    c = DR._parent_evidence_clause({"id": "MONDO:1", "name": "heavy metal poisoning"})
    assert "4 clinical trials" in c and "1 GWAS association" in c and "1 approved drug" in c
    evidence.reset()


def test_mechanism_synthesis(monkeypatch):
    """Cross-entity join: indicated drugs that also target a cohort gene, matched
    by manifest URL (never by fuzzy name)."""
    from atlas.disease import render as DR
    from atlas.page import links
    monkeypatch.setattr(links, "drug_url",
                        lambda chembl_id=None, name=None: f"/atlas/drug/{(chembl_id or name or '').lower()}/")
    bundles = {
        "5": {"gene_count": 8},
        "10": {"drugs": [{"id": "CHEMBL1", "name": "Vemurafenib", "gene_targets": ["BRAF"]},
                         {"id": "CHEMBL2", "name": "OtherDrug", "gene_targets": ["XYZ"]}]},
        "_indicated_drugs": [{"name": "Vemurafenib", "url": "/atlas/drug/chembl1/"},
                             {"name": "Aspirin", "url": "/atlas/drug/aspirin/"}],
    }
    md = DR.r_mechanism_synthesis(bundles)
    assert "Mechanistic alignment" in md
    assert "1 of the 2 drugs indicated" in md
    assert "Vemurafenib" in md and "BRAF" in md
    assert "Aspirin" not in md                     # not a cohort drug → not aligned
    assert DR.r_mechanism_synthesis({"5": {}, "10": {}, "_indicated_drugs": []}) == ""


def test_clingen_disease_name_dedup():
    """ClinGen condition labels differing only by a trailing numeric subtype
    collapse to the shortest representative (Li-Fraumeni syndrome subsumes
    'Li-Fraumeni syndrome 1'); genuinely distinct names are kept."""
    from atlas.gene.sections.s06_variants import _dedup_disease_names as dd
    assert dd(["Li-Fraumeni syndrome", "Li-Fraumeni syndrome 1"]) == ["Li-Fraumeni syndrome"]
    assert dd(["Li-Fraumeni syndrome 1"]) == ["Li-Fraumeni syndrome 1"]   # nothing to subsume
    assert dd(["BRCA1-related cancer", "Lynch syndrome"]) == ["BRCA1-related cancer", "Lynch syndrome"]
    assert dd(["", None, "  "]) == []


def test_drug_target_landscape(monkeypatch):
    """Reverse-index reads over a drug's curated targets: disease reach + per-target
    competitor count (incl. self)."""
    from atlas.page import links
    monkeypatch.setattr(links, "gene_url",
                        lambda symbol=None, hgnc_id=None: f"/atlas/gene/{symbol}/")
    links._REVERSE = {
        "/atlas/gene/ABL1/": [["m", "/atlas/disease/cml/", "disease", "Genes"],
                              ["m", "/atlas/drug/dasatinib/", "drug", "Genes"],
                              ["m", "/atlas/drug/imatinib/", "drug", "Genes"]],
    }
    bundle = {"2": {"primary_targets": [{"gene_symbol": "ABL1", "source": "gtopdb"},
                                        {"gene_symbol": "OFF", "source": "chembl"}]}}
    land = links.drug_target_landscape(bundle)
    assert land["n_targets"] == 1                       # only the GtoPdb target counts
    assert land["reach"] == {"/atlas/disease/cml/"}
    assert land["per_target"] == [("ABL1", 1, 2)]       # 1 cohort disease, 2 drugs (incl self)
    links._REVERSE = {}


def test_causal_genes_clinvar_fallback():
    """When GenCC/OMIM are empty, a small ClinVar-route cohort becomes the causal
    tier; a large (polygenic) cohort does not overclaim."""
    from atlas.disease.cohort import causal_genes
    b = {"4": {}, "5": {"genes": [
        {"symbol": "ARSL", "evidence": {"clinvar": True}},
        {"symbol": "ARSD", "evidence": {"clinvar": True}},
        {"symbol": "G", "evidence": {"gwas": True}}]}}
    assert causal_genes(b) == [("ARSD", "ClinVar-linked"), ("ARSL", "ClinVar-linked")]
    big = {"4": {}, "5": {"genes": [{"symbol": f"G{i}", "evidence": {"clinvar": True}}
                                    for i in range(8)]}}
    assert causal_genes(big) == []                 # >5 → no causal overclaim


def test_molecular_basis_block():
    from atlas.disease import render as DR
    b = {"4": {}, "5": {
        "genes": [{"symbol": "ARSL", "evidence": {"clinvar": True}}],
        "molecular_basis": [{"symbol": "ARSL", "protein_name": "Arylsulfatase L",
            "uniprot": "P51690", "family": "Belongs to the sulfatase family.",
            "function": "Exhibits arylsulfatase activity.",
            "subcellular": "Golgi apparatus.", "cofactor": "Binds 1 Ca(2+) ion.",
            "disease": "Chondrodysplasia punctata 1 (CDPX1)."}]}}
    md = DR.r_molecular_basis(b)
    assert "## Molecular basis" in md and "ARSL — Arylsulfatase L" in md
    assert "sulfatase family" in md and "Chondrodysplasia punctata 1 (CDPX1)" in md
    assert "Golgi apparatus; Binds 1 Ca(2+) ion" in md     # cleaned join, no stray period
    assert DR.r_molecular_basis({"4": {}, "5": {}}) == ""   # no causal genes → elides


def test_name_overlap_gate():
    from atlas.disease.anchors import _name_overlap
    assert _name_overlap("Chondrodysplasia punctata, brachytelephalangic, autosomal",
                         "Brachytelephalangic chondrodysplasia punctata") == 1.0
    assert _name_overlap("Chondrodysplasia punctata", "Cystic fibrosis") == 0.0


def test_pdb_table_shows_title():
    """The PDB experimental-structures table carries the structure Title column."""
    from atlas.gene import render as R
    b4 = {"reviewed_uniprot": ["P0"], "pdb_count": 1, "pdb": [
        {"id": "1M17", "title": "EGFR tyrosine kinase domain with erlotinib",
         "method": "X-RAY DIFFRACTION", "resolution": "2.6"}]}
    md = R.r_structure(b4)
    assert "| PDB | Title | Method | Resolution (Å) |" in md
    assert "1M17" in md and "EGFR tyrosine kinase domain with erlotinib" in md


def test_cross_species_homologs_render():
    """Diamond cross-species homologs render as an organism-labelled table (beyond
    Compara), with UniProt links and % sequence identity."""
    from atlas.gene import render as R
    b = {"ortholog_count": 1,
         "orthologs": [{"organism": "mus_musculus", "symbol": "Egfr", "id": "ENSMUSG..."}],
         "cross_species_homologs": [
             {"organism": "Pongo abelii", "accession": "Q5RB22", "similarity": 0.98, "source": "Diamond"},
             {"organism": "Bos taurus", "accession": "P04412", "similarity": 0.541, "source": "Diamond"}]}
    md = R.r_orthologs(b)
    assert "### Additional cross-species homologs {#cross-species-homologs}" in md
    # 1-decimal identity (not rounded), Diamond-only framing, and a "% sequence
    # identity" header — ESM2 (which over-called, e.g. EGFR→bovine FKBP9) is gone.
    assert "Pongo abelii" in md and "Q5RB22" in md and "98.0%" in md and "54.1%" in md
    assert "beyond Ensembl Compara" in md
    assert "Diamond" in md and "identity" in md and "ESM2" not in md


def test_ncrna_layer_renders_disease_interaction_drug_function():
    """The §14 non-coding RNA layer renders its four blocks from a bundle."""
    from atlas.gene import render as R
    b14 = {
        "rfam": [{"rfam_id": "RF00658", "rfam_description": "microRNA mir-21", "rna_type": "miRNA"}],
        "go": [{"id": "GO:0016442", "type": "cellular_component", "name": "RISC complex"}],
        "diseases": [{"disease_name": "Stomach Neoplasms", "causality": "Yes",
                      "validated_method": "qRT-PCR//Western Blot", "category": "LncRNA"}],
        "disease_total": 246,
        "interactions": [{"partner_name": "hsa-miR-412", "partner_type": "miRNA",
                          "level": "RNA-RNA", "datasource": "Literature mining"}],
        "interaction_total": 2968,
        "drugs": [{"drug_name": "Cisplatin", "relation": "drug_resistance",
                   "effect": "resistant", "condition": "cell line"}],
        "drug_total": 81,
    }
    fn = R.r_ncrna_function(b14)
    assert "RF00658" in fn and "RISC complex" in fn and "microRNA mir-21" in fn
    dis = R.r_ncrna_disease(b14)
    assert "Stomach Neoplasms" in dis and "246" in dis and "qRT-PCR, Western Blot" in dis
    inter = R.r_ncrna_interactions(b14)
    assert "hsa-miR-412" in inter and "2,968" in inter
    drg = R.r_ncrna_drugs(b14)
    assert "Cisplatin" in drg and "drug_resistance" in drg
    # all elide on an empty bundle
    assert R.r_ncrna_function({}) == "" and R.r_ncrna_disease({}) == ""


def test_cellphonedb_ligand_receptor_render():
    """§8 renders CellPhoneDB ligand–receptor pairs as a distinct subsection,
    with the partner gene and this gene's role (ligand/receptor)."""
    from atlas.gene import render as R
    b = {"cellphonedb": [{"partner": "EREG", "role": "receptor",
                          "classification": "Signaling by Epidermal growth factor"},
                         {"partner": "TGFA", "role": "receptor",
                          "classification": "Signaling by Transforming growth factor"}],
         "cellphonedb_count": 7}
    md = R.r_interactions(b)
    assert "Ligand–receptor pairs (CellPhoneDB)" in md and "{#cellphonedb}" in md
    assert "EREG" in md and "receptor" in md
    assert R.r_interactions({}).find("CellPhoneDB") == -1   # elides with no data


def test_drug_moa_targets_render_and_lead():
    """ChEMBL mechanism-of-action gives RNA therapeutics a target (#49): MOA block
    + 'targeting X' in the lead, even with no GtoPdb/bioactivity target."""
    from atlas.drug import render as DR
    from atlas.page.drug_declarative import _targets_clause
    b2 = {"primary_targets": [], "bioactivity_target_count": 0,
          "mechanisms": [{"mechanism_of_action": "PCSK9 mRNA RNAi inhibitor",
                          "action_type": "RNAI INHIBITOR", "target_name": "PCSK9 mRNA",
                          "target_type": "NUCLEIC-ACID"}],
          "mechanism_genes": [{"hgnc_id": "HGNC:20001", "gene_symbol": "PCSK9"}]}
    md = DR.r_targets(b2)
    assert "Mechanism of action (ChEMBL curated)" in md and "PCSK9 mRNA RNAi inhibitor" in md
    assert "No target linkage" not in md            # MOA counts as linkage
    assert _targets_clause(b2) == " targeting PCSK9"


def test_cross_species_homolog_virus_classifier():
    """Viral taxa (v-erbB carriers) are detected so they de-prioritise below the
    cellular orthologs; no cellular organism carries 'virus' in its name."""
    from atlas.gene.sections.s05_orthologs import _is_virus
    assert _is_virus("Avian leukosis virus")
    assert _is_virus("Avian erythroblastosis virus (strain ES4)")
    assert not _is_virus("Macaca mulatta")
    assert not _is_virus("Bos taurus") and not _is_virus("")


def test_clinical_trials_sponsor_and_intervention_drugs():
    """Top-trials table carries the sponsor; intervention drugs not in ChEMBL
    (e.g. daraxonrasib) surface in their own block."""
    from atlas.disease import render as DR
    b = {"trial_count": 1,
         "top_trials": [{"id": "NCT1", "title": "Study of X", "phase": "PHASE3",
                         "status": "RECRUITING", "sponsor": "Revolution Medicines, Inc."}],
         "trial_drugs": [],
         "trial_intervention_drugs": [{"name": "Daraxonrasib", "max_phase": 3, "phase": "PHASE3"}]}
    md = DR.r_clinical_trials(b)
    assert "| NCT | Phase | Status | Sponsor | Title |" in md
    assert "Revolution Medicines, Inc." in md
    assert "### Other drugs named in trials {#trial-intervention-drugs}" in md
    assert "Daraxonrasib" in md


def test_disease_no_about_block_in_body():
    """The medical disclaimer lives in the web theme's footer, NOT a body block —
    disease pages must not emit a duplicate '**About this page**'."""
    import atlas.pipeline as P
    md = P.assemble_page("m", "", "## Clinical features {#clinical}\n\nb",
                         {"entity_type": "disease"},
                         bundle={"1": {"canonical_name": "Marfan syndrome"}})
    assert "About this page" not in md
    assert "not as medical advice" not in md


def test_noncoding_genesets_function_block():
    """Non-coding genes surface MSigDB gene-set membership in the Function zone."""
    from atlas.gene import render as R
    b7 = {"msigdb": [{"id": "M1", "name": "TFEB_TARGET_GENES"},
                     {"id": "M2", "name": "chr6q26"}], "msigdb_total": 2}
    md = R.r_noncoding_genesets(b7)
    assert "Gene-set membership (MSigDB)" in md
    assert "member of 2 curated MSigDB gene sets" in md
    assert "`TFEB_TARGET_GENES`" in md
    assert R.r_noncoding_genesets({}) == ""          # no sets → elide


def test_noncoding_overlap_gene_bullet(monkeypatch):
    """A non-coding gene names + links the overlapping protein-coding gene whose
    positional variant/disease data it would otherwise inherit."""
    from atlas.page import at_a_glance as AAG, links
    monkeypatch.setattr(links, "gene_url", lambda symbol=None, **k: f"/atlas/gene/{symbol}/")
    md = AAG.at_a_glance({"_noncoding": "lncRNA", "6": {"overlap_genes": ["QKI"]}})
    assert "non-coding (lncRNA)" in md
    assert "overlapping protein-coding gene [QKI](/atlas/gene/QKI/)" in md
