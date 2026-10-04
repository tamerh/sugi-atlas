"""§9 — pharmacogenomics. Drug-level PGx = CPIC / DPWG genotype-guided dosing
guidelines for THIS drug (drug × metabolizing-gene, e.g. atorvastatin × SLCO1B1
statin myopathy). NOT the drug's *target* gene — PGx is about the patient's
pharmacogene genotype.

There's no direct chembl_molecule→pharmgkb edge, but biobtree is a graph: the
PharmGKB chemical node cross-refs PubChem, so the drug reaches it by ID-join
(no fragile name matching):

    chembl_molecule >> pubchem >> pharmgkb >> pharmgkb_guideline

The intermediate `pharmgkb` chemical node also carries clinical/variant
annotation counts (gene-keyed annotations live on the gene pages), and its
`drug_labels` — PharmGKB's annotations of regulatory drug labels (FDA / EMA /
PMDA / HCSC / Swissmedic) with the label's PGx level ("Testing Required",
"Actionable PGx", ...) and genes. Its curated PK/PD pathways hang off it via
pharmgkb >> pharmgkb_pathway. Empty for drugs with no curated PGx (e.g. newer
targeted agents)."""
from atlas.biobtree import map_all, entry
from atlas.section import Section

_CHEMICAL_CHAIN = ">>chembl_molecule>>pubchem>>pharmgkb"
_GUIDELINE_CHAIN = ">>chembl_molecule>>pubchem>>pharmgkb>>pharmgkb_guideline"
_CLINICAL_CHAIN = ">>hgnc>>pharmgkb_clinical"   # gene → its PharmGKB clinical annotations
_VARANN_CHAIN = ">>hgnc>>pharmgkb_var_annotation"  # gene → per-publication variant annotations
_PATHWAY_CHAIN = ">>pharmgkb>>pharmgkb_pathway"   # PharmGKB chemical → its PK/PD pathways


def _int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def _clinical_annotations(pgx_id, pgx_name):
    """The drug's PharmGKB clinical annotations (variant × gene × type × level ×
    phenotype). There is no direct drug→annotation edge — annotations are
    variant/gene-keyed — so we read the drug's related PGx genes off its chemical
    node, fan each (flagged ClinicalAnnotation) through >>hgnc>>pharmgkb_clinical,
    and keep the rows whose `chemicals` include this drug."""
    if not (pgx_id and pgx_name):
        return []
    try:
        ce = entry(pgx_id, "pharmgkb")
        related = ((ce.get("Attributes") or {}).get("Pharmgkb") or {}).get("related_genes") or []
    except Exception:
        return []
    out, seen = [], set()
    for g in related:
        sym = g.get("gene_symbol")
        ev = g.get("evidence_type") or ""
        # PharmGKB renamed clinical annotations → "SummaryAnnotation"; matching only
        # the old token rendered clinical annotations on 0 of 4,683 drug pages.
        if not sym or not ("ClinicalAnnotation" in ev or "SummaryAnnotation" in ev):
            continue
        try:
            rows = map_all(sym, _CLINICAL_CHAIN)
        except Exception:
            continue
        for r in rows:
            chems = {c.strip().lower() for c in (r.get("chemicals") or "").split(";")}
            if pgx_name not in chems:
                continue
            key = (r.get("variant"), sym, r.get("type"), r.get("phenotypes"))
            if key in seen:
                continue
            seen.add(key)
            out.append({
                "variant": r.get("variant"),
                "gene": sym,
                "type": r.get("type"),
                "level": r.get("level_of_evidence"),
                "phenotypes": r.get("phenotypes"),
            })
    # PharmGKB level of evidence: 1 strongest → 4 weakest; sort that way, then gene.
    out.sort(key=lambda x: (str(x.get("level") or "9"), x.get("gene") or ""))
    return out


def _variant_annotations(pgx_id, pgx_name):
    """The drug's per-publication PharmGKB variant annotations — the raw evidence
    layer BENEATH the clinical annotations: one row per published finding, each
    with a plain-English `sentence` + PMID. Same fan-out as _clinical_annotations
    (read the drug's related PGx genes, fan each through
    >>hgnc>>pharmgkb_var_annotation), keeping only the SIGNIFICANT rows that name
    this drug. Deduped by (variant, gene, pmid). This is the data the old
    'see PharmGKB' tease pointed at — now aggregated."""
    if not (pgx_id and pgx_name):
        return []
    try:
        ce = entry(pgx_id, "pharmgkb")
        related = ((ce.get("Attributes") or {}).get("Pharmgkb") or {}).get("related_genes") or []
    except Exception:
        return []
    out, seen = [], set()
    for g in related:
        sym = g.get("gene_symbol")
        if not sym:
            continue
        try:
            rows = map_all(sym, _VARANN_CHAIN)
        except Exception:
            continue
        for r in rows:
            if r.get("significance") != "yes":          # positive findings only
                continue
            drugs = {d.strip().lower()
                     for d in (r.get("drugs") or "").replace(",", ";").split(";")}
            if pgx_name not in drugs:
                continue
            key = (r.get("variant"), sym, r.get("pmid"))
            if key in seen:
                continue
            seen.add(key)
            out.append({
                "variant": r.get("variant"),
                "gene": sym,
                "category": r.get("phenotype_category"),
                "direction": r.get("direction_of_effect"),
                "sentence": r.get("sentence"),
                "pmid": r.get("pmid"),
            })
    out.sort(key=lambda x: (x.get("gene") or "", x.get("variant") or ""))
    return out


def _drug_labels(pgx_id):
    """The PharmGKB chemical's regulatory drug-label annotations: one row per
    (agency, label annotation) with the PGx level and the genes it names. Rows
    with neither a PGx level nor a gene carry no pharmacogenomic content (an
    annotated label that names no gene) and are dropped. Field names as served
    by biobtree's PharmgkbDrugLabel."""
    if not pgx_id:
        return []
    try:
        ce = entry(pgx_id, "pharmgkb")
        labels = ((ce.get("Attributes") or {}).get("Pharmgkb") or {}).get("drug_labels") or []
    except Exception:
        return []
    out = []
    for l in labels:
        level = (l.get("testing_level") or "").strip()
        genes = [g for g in (l.get("genes") or []) if g]
        if not level and not genes:
            continue
        out.append({
            "label_id": l.get("label_id"),
            "source": (l.get("source") or "").strip(),
            "level": level,
            "genes": genes,
            "has_prescribing_info": bool(l.get("has_prescribing_info")),
            "has_dosing_info": bool(l.get("has_dosing_info")),
            "has_alternate_drug": bool(l.get("has_alternate_drug")),
        })
    return out


def _pathways(pgx_id):
    """PharmGKB's curated pharmacokinetic / pharmacodynamic pathways for the
    chemical (pharmgkb >> pharmgkb_pathway), sorted by name."""
    if not pgx_id:
        return []
    try:
        hits = map_all(pgx_id, _PATHWAY_CHAIN)
    except Exception:
        return []
    out = [{"id": r.get("id"), "name": (r.get("name") or "").strip(),
            "pk": r.get("is_pharmacokinetic") == "true",
            "pd": r.get("is_pharmacodynamic") == "true"}
           for r in hits if r.get("id") and (r.get("name") or "").strip()]
    out.sort(key=lambda p: p["name"].lower())
    return out


def collect(a):
    chem = map_all(a.chembl_id, _CHEMICAL_CHAIN, cap=1)
    c0 = chem[0] if chem else {}
    pgx_name = (c0.get("name") or "").strip().lower()
    guidelines = [{
        "id": r.get("id"),
        "name": r.get("name"),
        "source": r.get("source"),               # CPIC / DPWG / ...
        "genes": r.get("gene_symbols"),
        "chemicals": r.get("chemical_names"),
        "has_dosing": r.get("has_dosing_info") == "true",
        "has_recommendation": r.get("has_recommendation") == "true",
    } for r in map_all(a.chembl_id, _GUIDELINE_CHAIN)]
    return {
        "section": "09_pharmacogenomics",
        "pharmgkb_chemical_id": c0.get("id"),
        "clinical_annotation_count": _int(c0.get("clinical_annotation_count")),
        "variant_annotation_count": _int(c0.get("variant_annotation_count")),
        "clinical_annotations": _clinical_annotations(c0.get("id"), pgx_name),
        "variant_annotations": _variant_annotations(c0.get("id"), pgx_name),
        "guidelines": guidelines,
        "guideline_count": len(guidelines),
        "drug_labels": _drug_labels(c0.get("id")),
        "pathways": _pathways(c0.get("id")),
    }


SECTION = Section(
    id="9", name="pharmacogenomics",
    description=("Drug-level pharmacogenomics: CPIC / DPWG genotype-guided dosing "
                 "guidelines (drug × pharmacogene) via the graph path "
                 "chembl_molecule→pubchem→pharmgkb→pharmgkb_guideline, plus "
                 "PharmGKB regulatory drug-label PGx annotations and PK/PD pathways"),
    needs=("chembl_id",),
    produces=("pharmgkb_chemical_id", "clinical_annotation_count",
              "variant_annotation_count", "clinical_annotations",
              "variant_annotations", "guidelines", "guideline_count",
              "drug_labels", "pathways"),
    datasets=("chembl_molecule", "pubchem", "pharmgkb", "pharmgkb_guideline",
              "hgnc", "pharmgkb_clinical", "pharmgkb_var_annotation",
              "pharmgkb_pathway"),
    chains=(_CHEMICAL_CHAIN, _GUIDELINE_CHAIN, _CLINICAL_CHAIN, _VARANN_CHAIN,
            _PATHWAY_CHAIN),
    collect_fn=collect,
)
