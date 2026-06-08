"""§9 — pharmacogenomics. Drug-level PGx = CPIC / DPWG genotype-guided dosing
guidelines for THIS drug (drug × metabolizing-gene, e.g. atorvastatin × SLCO1B1
statin myopathy). NOT the drug's *target* gene — PGx is about the patient's
pharmacogene genotype.

There's no direct chembl_molecule→pharmgkb edge, but biobtree is a graph: the
PharmGKB chemical node cross-refs PubChem, so the drug reaches it by ID-join
(no fragile name matching):

    chembl_molecule >> pubchem >> pharmgkb >> pharmgkb_guideline

The intermediate `pharmgkb` chemical node also carries clinical/variant
annotation counts (gene-keyed annotations live on the gene pages). Empty for
drugs with no curated PGx (e.g. newer targeted agents)."""
from atlas.biobtree import map_all, entry
from atlas.section import Section

_CHEMICAL_CHAIN = ">>chembl_molecule>>pubchem>>pharmgkb"
_GUIDELINE_CHAIN = ">>chembl_molecule>>pubchem>>pharmgkb>>pharmgkb_guideline"
_CLINICAL_CHAIN = ">>hgnc>>pharmgkb_clinical"   # gene → its PharmGKB clinical annotations


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
        if not sym or "ClinicalAnnotation" not in (g.get("evidence_type") or ""):
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
        "guidelines": guidelines,
        "guideline_count": len(guidelines),
    }


SECTION = Section(
    id="9", name="pharmacogenomics",
    description=("Drug-level pharmacogenomics: CPIC / DPWG genotype-guided dosing "
                 "guidelines (drug × pharmacogene) via the graph path "
                 "chembl_molecule→pubchem→pharmgkb→pharmgkb_guideline"),
    needs=("chembl_id",),
    produces=("pharmgkb_chemical_id", "clinical_annotation_count",
              "variant_annotation_count", "clinical_annotations",
              "guidelines", "guideline_count"),
    datasets=("chembl_molecule", "pubchem", "pharmgkb", "pharmgkb_guideline",
              "hgnc", "pharmgkb_clinical"),
    chains=(_CHEMICAL_CHAIN, _GUIDELINE_CHAIN, _CLINICAL_CHAIN),
    collect_fn=collect,
)
