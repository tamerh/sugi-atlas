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
from atlas.biobtree import map_all
from atlas.section import Section

_CHEMICAL_CHAIN = ">>chembl_molecule>>pubchem>>pharmgkb"
_GUIDELINE_CHAIN = ">>chembl_molecule>>pubchem>>pharmgkb>>pharmgkb_guideline"


def _int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def collect(a):
    chem = map_all(a.chembl_id, _CHEMICAL_CHAIN, cap=1)
    c0 = chem[0] if chem else {}
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
              "variant_annotation_count", "guidelines", "guideline_count"),
    datasets=("chembl_molecule", "pubchem", "pharmgkb", "pharmgkb_guideline"),
    chains=(_CHEMICAL_CHAIN, _GUIDELINE_CHAIN),
    collect_fn=collect,
)
