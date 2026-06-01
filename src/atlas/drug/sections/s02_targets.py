"""§2 — targets. Primary = GtoPdb curated mechanism targets (gene + action +
pAffinity; covers antibodies). Secondary = the raw ChEMBL bioactivity target
set (count + sample). Each primary target is annotated with its DepMap
cancer-dependency signal (>>hgnc>>depmap) — the target-tractability view:
"is this drug's target a cancer dependency?". Gene symbols become /atlas/gene/
links once links.py exists (post-build pass)."""
from atlas.biobtree import map_all
from atlas.section import Section


def _depmap(hgnc_id):
    """DepMap per-gene CRISPR fitness aggregate; {} on miss (non-cancer-screened
    genes / no hgnc)."""
    if not hgnc_id:
        return {}
    r = map_all(hgnc_id, ">>hgnc>>depmap", cap=1)
    return r[0] if r else {}


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
    return {
        "section": "02_targets",
        "primary_targets": primary,
        "primary_source": (a.targets[0].source if a.targets else None),
        "bioactivity_target_count": len(bioactivity),
        "bioactivity_targets": bioactivity[:30],
    }


SECTION = Section(
    id="2", name="targets",
    description=("Primary mechanism targets (GtoPdb-curated: gene + action + "
                 "pAffinity; covers antibodies) annotated with DepMap cancer-"
                 "dependency + secondary ChEMBL bioactivity target set"),
    needs=("targets", "bioactivity_targets"),
    produces=("primary_targets", "primary_source", "bioactivity_target_count",
              "bioactivity_targets"),
    datasets=("gtopdb_ligand", "gtopdb_interaction", "gtopdb", "uniprot",
              "hgnc", "chembl_target", "depmap"),
    chains=(">>gtopdb_ligand>>gtopdb_interaction>>gtopdb>>uniprot>>hgnc",
            ">>chembl_molecule>>chembl_target", ">>hgnc>>depmap"),
    collect_fn=collect,
)
