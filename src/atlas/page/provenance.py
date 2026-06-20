"""Per-page provenance sidecar — turn Atlas's internal data trail into a
machine-fetchable artifact so AI agents can cite individual facts with their
upstream source, not just the page.

Schema-wrapped as schema.org Dataset; emitted alongside page.md as
`provenance.json` in each gene's dist directory.

For every section the collector ran, we expose:
  - the biobtree chains used (the queries that produced the facts)
  - the datasets traversed
  - the upstream primary source (NCBI ClinVar, EBI Reactome, etc.) with URL
  - which bundle keys this section produced

Plus the resolved anchors for the entity (HGNC id, Ensembl id, canonical
UniProt, canonical transcript) so an agent can quickly answer
"what did Atlas anchor this gene to?"
"""
from atlas.gene.sections import REGISTRY

BASE_URL = "https://sugi.bio/atlas"

# Map biobtree dataset name -> (human display name, upstream source URL).
# These point at the source-of-record (project home), not at per-id deep links —
# deep links per fact are a follow-up (Path B item #2: per-fact HTML anchors).
UPSTREAM = {
    # identifiers
    "hgnc":              ("HGNC",                                  "https://www.genenames.org/"),
    "ensembl":           ("Ensembl",                               "https://www.ensembl.org/"),
    "entrez":            ("NCBI Gene (Entrez)",                    "https://www.ncbi.nlm.nih.gov/gene/"),
    "neighborentrez":    ("NCBI Gene (neighbors)",                 "https://www.ncbi.nlm.nih.gov/gene/"),
    "orthologentrez":    ("NCBI Gene (orthologs)",                 "https://www.ncbi.nlm.nih.gov/gene/"),
    "mim":               ("OMIM",                                  "https://www.omim.org/"),
    "rnacentral":        ("RNAcentral",                            "https://rnacentral.org/"),
    "ncrna_disease":     ("LncRNADisease / HMDD",                  "http://www.rnanut.net/lncrnadisease/"),
    "ncrna_interaction": ("NPInter v5",                            "http://bigdata.ibp.ac.cn/npinter5/"),
    "ncrna_drug":        ("ncRNADrug",                             "http://www.jianglab.cn/ncRNADrug/"),
    "chembl_mechanism":  ("ChEMBL (mechanism of action)",         "https://www.ebi.ac.uk/chembl/"),
    "cellphonedb":       ("CellPhoneDB",                           "https://www.cellphonedb.org/"),
    # transcripts / exons
    "transcript":        ("Ensembl",                               "https://www.ensembl.org/"),
    "exon":              ("Ensembl",                               "https://www.ensembl.org/"),
    "refseq":            ("NCBI RefSeq",                           "https://www.ncbi.nlm.nih.gov/refseq/"),
    "ccds":              ("CCDS",                                  "https://www.ncbi.nlm.nih.gov/CCDS/"),
    # proteins / domains
    "uniprot":           ("UniProt",                               "https://www.uniprot.org/"),
    "interpro":          ("InterPro",                              "https://www.ebi.ac.uk/interpro/"),
    "pfam":              ("Pfam",                                  "https://www.ebi.ac.uk/interpro/entry/pfam/"),
    "antibody":          ("SAbDab / Thera-SAbDab (Oxford OPIG)",   "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/"),
    "ufeature":          ("UniProt sequence features",             "https://www.uniprot.org/"),
    "brenda":            ("BRENDA enzyme database",                "https://www.brenda-enzymes.org/"),
    "rhea":              ("Rhea reaction database",                "https://www.rhea-db.org/"),
    # structure
    "pdb":               ("RCSB PDB",                              "https://www.rcsb.org/"),
    "alphafold":         ("AlphaFold DB",                          "https://alphafold.ebi.ac.uk/"),
    # orthology
    "ortholog":          ("Ensembl Compara",                       "https://www.ensembl.org/info/genome/compara/"),
    "panelapp_gene":     ("Genomics England PanelApp",             "https://panelapp.genomicsengland.co.uk/"),
    "mgi":               ("MGI (Mouse Genome Informatics)",        "https://www.informatics.jax.org/"),
    "alliance_phenotype": ("Alliance of Genome Resources",         "https://www.alliancegenome.org/"),
    "paralog":           ("Ensembl Compara",                       "https://www.ensembl.org/info/genome/compara/"),
    # variants
    "clinvar":           ("ClinVar",                               "https://www.ncbi.nlm.nih.gov/clinvar/"),
    "spliceai":          ("SpliceAI",                              "https://spliceailookup.broadinstitute.org/"),
    "alphamissense":     ("AlphaMissense",                         "https://alphamissense.hegelab.org/"),
    "dbsnp":             ("dbSNP",                                 "https://www.ncbi.nlm.nih.gov/snp/"),
    # pathways
    "reactome":          ("Reactome",                              "https://reactome.org/"),
    "msigdb":            ("MSigDB",                                "https://www.gsea-msigdb.org/gsea/msigdb/"),
    "go":                ("Gene Ontology",                         "http://geneontology.org/"),
    # interactions
    "string_interaction":("STRING",                                "https://string-db.org/"),
    "intact":            ("IntAct",                                "https://www.ebi.ac.uk/intact/"),
    "biogrid_interaction":("BioGRID",                              "https://thebiogrid.org/"),
    "signor":            ("SIGNOR",                                "https://signor.uniroma2.it/"),
    "esm2_similarity":   ("ESM2 (Meta AI)",                        "https://github.com/facebookresearch/esm"),
    "diamond_similarity":("DIAMOND alignment",                     "https://github.com/bbuchfink/diamond"),
    "taxonomy":          ("UniProt Taxonomy",                      "https://www.uniprot.org/taxonomy"),
    # TF regulation + post-transcriptional regulation
    "collectri":         ("CollecTRI",                             "https://github.com/saezlab/CollecTRI"),
    "jaspar":            ("JASPAR",                                "https://jaspar.genereg.net/"),
    "pubmed":            ("PubMed (NCBI)",                         "https://pubmed.ncbi.nlm.nih.gov/"),
    "mirdb":             ("miRDB",                                 "https://mirdb.org/"),
    # drugs
    "chembl_target":     ("ChEMBL",                                "https://www.ebi.ac.uk/chembl/"),
    "chembl_molecule":   ("ChEMBL",                                "https://www.ebi.ac.uk/chembl/"),
    "chembl_activity":   ("ChEMBL",                                "https://www.ebi.ac.uk/chembl/"),
    "chembl_assay":      ("ChEMBL",                                "https://www.ebi.ac.uk/chembl/"),
    "chembl_document":   ("ChEMBL",                                "https://www.ebi.ac.uk/chembl/"),
    "cellosaurus":       ("Cellosaurus (SIB cell-line index)",     "https://www.cellosaurus.org/"),
    "patent_compound":   ("PubChem patent-compound index",         "https://pubchem.ncbi.nlm.nih.gov/"),
    "pharmgkb_gene":     ("PharmGKB",                              "https://www.pharmgkb.org/"),
    "pharmgkb_clinical": ("PharmGKB Clinical Annotations",         "https://www.pharmgkb.org/clinicalAnnotations"),
    "pharmgkb_variant":  ("PharmGKB Variant Annotations",          "https://www.pharmgkb.org/variantAnnotations"),
    "pharmgkb_var_annotation": ("PharmGKB Variant Annotations",    "https://www.pharmgkb.org/variantAnnotations"),
    "pharmgkb_guideline": ("PharmGKB / CPIC dosing guidelines",    "https://www.pharmgkb.org/guidelineAnnotations"),
    "pharmgkb":          ("PharmGKB",                              "https://www.pharmgkb.org/"),
    "faers":             ("openFDA FAERS",                         "https://open.fda.gov/data/faers/"),
    "faers_reaction":    ("openFDA FAERS",                         "https://open.fda.gov/data/faers/"),
    "drugcentral":       ("DrugCentral",                           "https://drugcentral.org/"),
    "chebi":             ("ChEBI",                                 "https://www.ebi.ac.uk/chebi/"),
    "pubchem":           ("PubChem",                               "https://pubchem.ncbi.nlm.nih.gov/"),
    "gtopdb":            ("Guide to Pharmacology (GtoPdb)",        "https://www.guidetopharmacology.org/"),
    "gtopdb_ligand":     ("Guide to Pharmacology (GtoPdb)",        "https://www.guidetopharmacology.org/"),
    "gtopdb_interaction":("Guide to Pharmacology (GtoPdb)",        "https://www.guidetopharmacology.org/"),
    "bindingdb":         ("BindingDB",                             "https://www.bindingdb.org/"),
    "pubchem_activity":  ("PubChem BioAssay",                      "https://pubchem.ncbi.nlm.nih.gov/"),
    "ctd_gene_interaction":("Comparative Toxicogenomics Database (CTD)", "http://ctdbase.org/"),
    "clinical_trials":   ("ClinicalTrials.gov",                    "https://clinicaltrials.gov/"),
    # diseases / phenotypes
    "efo":               ("EFO (Experimental Factor Ontology)",    "https://www.ebi.ac.uk/efo/"),
    "civic_evidence":    ("CIViC (clinical evidence)",             "https://civicdb.org/"),
    "mondo":             ("Mondo Disease Ontology",                "https://mondo.monarchinitiative.org/"),
    "gencc":             ("GenCC",                                 "https://thegencc.org/"),
    "orphanet":          ("Orphanet",                              "https://www.orpha.net/"),
    "hpo":               ("Human Phenotype Ontology",              "https://hpo.jax.org/"),
    "gwas":              ("GWAS Catalog",                          "https://www.ebi.ac.uk/gwas/"),
    "gwas_study":        ("GWAS Catalog",                          "https://www.ebi.ac.uk/gwas/"),
    "mesh":              ("MeSH (NLM Medical Subject Headings)",   "https://www.ncbi.nlm.nih.gov/mesh/"),
    "intogen":           ("intOGen (cancer driver catalog)",       "https://www.intogen.org/"),
    # Mondo OBO cross-ontology xrefs (biobtree 2026-06-01)
    "doid":              ("Disease Ontology",                      "https://disease-ontology.org/"),
    "alliance_disease":  ("Alliance of Genome Resources",          "https://www.alliancegenome.org/"),
    "sctid":             ("SNOMED CT",                             "https://www.snomed.org/"),
    "umls":              ("UMLS (NLM Unified Medical Language)",   "https://www.nlm.nih.gov/research/umls/"),
    "ncit":              ("NCI Thesaurus",                         "https://ncithesaurus.nci.nih.gov/"),
    "medgen":            ("NCBI MedGen",                           "https://www.ncbi.nlm.nih.gov/medgen/"),
    "icd10cm":           ("ICD-10-CM",                             "https://www.cdc.gov/nchs/icd/icd-10-cm.htm"),
    "icd11":             ("ICD-11",                                "https://icd.who.int/en"),
    "gard":              ("Genetic and Rare Diseases Information Center (GARD)", "https://rarediseases.info.nih.gov/"),
    "meddra":            ("MedDRA",                                "https://www.meddra.org/"),
    "nord":              ("NORD (National Org for Rare Disorders)", "https://rarediseases.org/"),
    "uberon":            ("Uberon (multi-species anatomy ontology)", "https://obophenotype.github.io/uberon/"),
    "civic":             ("CIViC (Clinical Interpretation of Variants in Cancer)", "https://civicdb.org/"),
    "civic_variant":     ("CIViC (curated variants)",              "https://civicdb.org/"),
    "clingen_gene_validity": ("ClinGen Gene-Disease Validity",     "https://search.clinicalgenome.org/kb/gene-validity"),
    "clingen_dosage":    ("ClinGen Gene Dosage Map",               "https://search.clinicalgenome.org/kb/gene-dosage"),
    "clingen_variant":   ("ClinGen Variant Curation",              "https://erepo.clinicalgenome.org/evrepo/"),
    "corum":             ("CORUM Protein Complexes",               "https://mips.helmholtz-muenchen.de/corum/"),
    "depmap":            ("DepMap (Cancer Dependency Map)",        "https://depmap.org/portal/"),
    "depmap_dependency": ("DepMap (Cancer Dependency Map)",        "https://depmap.org/portal/"),
    "gnomad_constraint": ("gnomAD (gene constraint)",              "https://gnomad.broadinstitute.org/"),
    "generif":           ("NCBI GeneRIFs",                         "https://www.ncbi.nlm.nih.gov/gene/about-generif"),
    # expression
    "bgee":              ("Bgee",                                  "https://www.bgee.org/"),
    "bgee_evidence":     ("Bgee",                                  "https://www.bgee.org/"),
    "fantom5_gene":      ("FANTOM5",                               "https://fantom.gsc.riken.jp/5/"),
    "fantom5_promoter":  ("FANTOM5",                               "https://fantom.gsc.riken.jp/5/"),
    "scxa":              ("EBI Single Cell Expression Atlas",      "https://www.ebi.ac.uk/gxa/sc/"),
    "scxa_expression":   ("EBI Single Cell Expression Atlas",      "https://www.ebi.ac.uk/gxa/sc/"),
    "scxa_gene_experiment": ("EBI Single Cell Expression Atlas",   "https://www.ebi.ac.uk/gxa/sc/"),
    "hpa":               ("Human Protein Atlas",                   "https://www.proteinatlas.org/"),
    "hpa_expression":    ("Human Protein Atlas",                   "https://www.proteinatlas.org/"),
    "hpa_pathology":     ("Human Protein Atlas (Pathology)",       "https://www.proteinatlas.org/humanproteome/pathology"),
}

def _section_provenance(sec):
    """Provenance entry for a single Section dataclass."""
    upstreams = {}
    for ds in sec.datasets:
        if ds in UPSTREAM:
            name, url = UPSTREAM[ds]
            upstreams[name] = url
    return {
        "id": sec.id,
        "name": sec.name,
        "description": sec.description,
        "datasets": list(sec.datasets),
        "chains": list(sec.chains),
        "upstream_sources": [{"name": n, "url": u} for n, u in sorted(upstreams.items())],
        "produces": list(sec.produces),
    }

def build_provenance(bundle, meta=None, base_url=BASE_URL):
    """Compose the provenance Dataset blob for a gene page.

    bundle: full collector bundle ({section_id: bundle_dict}).
    meta:   optional dict carrying generated_at / atlas_version / biobtree_version
            (typically the same `meta` dict passed to assemble_page).
    """
    meta = meta or {}
    b1 = bundle.get("1") or {}
    b3 = bundle.get("3") or {}
    b6 = bundle.get("6") or {}
    symbol = b1.get("symbol") or "?"
    page_url = f"{base_url}/gene/{symbol}/"

    sections = [_section_provenance(REGISTRY[sid]) for sid in REGISTRY]

    return {
        "@context": "https://schema.org",
        "@type": "Dataset",
        "name": f"{symbol} — provenance trail",
        "description": (
            "Per-section provenance for Atlas's gene page: datasets touched, "
            "biobtree chains used, upstream-source URLs. Every fact in the "
            "page can be traced back through this trail to its primary source."
        ),
        "isPartOf": page_url,
        "url": f"{page_url}provenance.json",
        "generated_at": meta.get("generated_at"),
        "atlas_version": meta.get("atlas_version"),
        "biobtree_version": meta.get("biobtree_version"),
        "anchors": {
            "symbol": symbol,
            "hgnc_id": b1.get("hgnc_id"),
            "ensembl_id": b1.get("ensembl_id"),
            "canonical_uniprot": b3.get("canonical_uniprot"),
            "canonical_transcript": b6.get("canonical_transcript"),
        },
        "data_access": {
            "biobtree_project": "https://biobtree.org",
            "chain_syntax": "https://biobtree.org",
            "note": "Chains in each section's `chains` field can be replayed against any biobtree instance to verify the source data.",
        },
        "sections": sections,
    }

def as_provenance_string(prov):
    """Pretty-printed JSON for the provenance.json sidecar."""
    import json
    return json.dumps(prov, indent=2) + "\n"
