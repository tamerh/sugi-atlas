"""§2 — targets. Primary = GtoPdb curated mechanism targets (gene + action +
pAffinity; covers antibodies). Secondary = the raw ChEMBL bioactivity target
set (count + sample). Each primary target is annotated with its DepMap
cancer-dependency signal (>>hgnc>>depmap) — the target-tractability view:
"is this drug's target a cancer dependency?". Gene symbols become /atlas/gene/
links once links.py exists (post-build pass)."""
from atlas.biobtree import map_all, entry
from atlas.section import Section


def _depmap(hgnc_id):
    """DepMap per-gene CRISPR fitness aggregate; {} on miss (non-cancer-screened
    genes / no hgnc)."""
    if not hgnc_id:
        return {}
    r = map_all(hgnc_id, ">>hgnc>>depmap", cap=1)
    return r[0] if r else {}


def _hgnc_symbol(hgnc_id):
    """HGNC id -> primary symbol (the chembl_mechanism>>hgnc projection carries no
    symbol). None on miss."""
    try:
        h = (entry(hgnc_id, "hgnc").get("Attributes") or {}).get("Hgnc") or {}
        syms = h.get("symbols") or []
        return syms[0] if syms else None
    except Exception:
        return None


def _mechanisms(chembl_id):
    """ChEMBL curated mechanism-of-action (drug_mechanism table). The ONLY target
    edge for RNA therapeutics (siRNA/ASO), which carry no bioactivity target —
    e.g. inclisiran -> 'PCSK9 mRNA RNAi inhibitor' -> PCSK9. Returns (moa rows,
    resolved target genes)."""
    if not chembl_id:
        return [], []
    moa = [{"mechanism_of_action": r.get("mechanism_of_action"),
            "action_type": r.get("action_type"),
            "target_name": r.get("target_name"),
            "target_type": r.get("target_type")}
           for r in map_all(chembl_id, ">>chembl_molecule>>chembl_mechanism", cap=2)]
    genes = []
    seen = set()
    for r in map_all(chembl_id, ">>chembl_molecule>>chembl_mechanism>>hgnc", cap=2):
        hid = r.get("id")
        if hid and hid not in seen:
            seen.add(hid)
            genes.append({"hgnc_id": hid, "gene_symbol": _hgnc_symbol(hid)})
    return moa, genes


def collect(a):
    primary = []
    for t in a.targets:
        dm = _depmap(t.hgnc_id)
        primary.append({
            "gene_symbol": t.gene_symbol, "hgnc_id": t.hgnc_id, "uniprot": t.uniprot,
            "target_name": t.target_name, "target_type": t.target_type,
            "source": t.source, "action": t.action, "affinity": t.affinity,
            "dep_pct": dm.get("pct_dependent"),
            "dep_selective": dm.get("strongly_selective") == "true",
            "dep_common_essential": dm.get("common_essential") == "true",
        })
    bioactivity = [{"chembl_target_id": t.get("chembl_target_id"),
                    "name": t.get("name"), "type": t.get("type")}
                   for t in a.bioactivity_targets]
    mechanisms, mechanism_genes = _mechanisms(a.chembl_id)
    return {
        "section": "02_targets",
        "primary_targets": primary,
        "primary_source": (a.targets[0].source if a.targets else None),
        "bioactivity_target_count": len(bioactivity),
        "bioactivity_targets": bioactivity[:30],
        "mechanisms": mechanisms,
        "mechanism_genes": mechanism_genes,
    }


SECTION = Section(
    id="2", name="targets",
    description=("Primary mechanism targets (GtoPdb-curated: gene + action + "
                 "pAffinity; covers antibodies) annotated with DepMap cancer-"
                 "dependency + secondary ChEMBL bioactivity target set"),
    needs=("targets", "bioactivity_targets"),
    produces=("primary_targets", "primary_source", "bioactivity_target_count",
              "bioactivity_targets", "mechanisms", "mechanism_genes"),
    datasets=("gtopdb_ligand", "gtopdb_interaction", "gtopdb", "uniprot",
              "hgnc", "chembl_target", "chembl_mechanism", "depmap"),
    chains=(">>gtopdb_ligand>>gtopdb_interaction>>gtopdb>>uniprot>>hgnc",
            ">>chembl_molecule>>chembl_target", ">>hgnc>>depmap",
            ">>chembl_molecule>>chembl_mechanism",
            ">>chembl_molecule>>chembl_mechanism>>hgnc"),
    collect_fn=collect,
)
