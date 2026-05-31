"""§2 — targets. Primary = GtoPdb curated mechanism targets (gene + action +
pAffinity; covers antibodies). Secondary = the raw ChEMBL bioactivity target
set (count + sample). Both pre-resolved on the anchor; no extra calls. Gene
symbols become /atlas/gene/ links once links.py exists (post-build pass)."""
from atlas.section import Section


def collect(a):
    primary = [{
        "gene_symbol": t.gene_symbol, "hgnc_id": t.hgnc_id, "uniprot": t.uniprot,
        "target_name": t.target_name, "target_type": t.target_type,
        "source": t.source, "action": t.action, "affinity": t.affinity,
    } for t in a.targets]
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
                 "pAffinity; covers antibodies) + secondary ChEMBL bioactivity "
                 "target set"),
    needs=("targets", "bioactivity_targets"),
    produces=("primary_targets", "primary_source", "bioactivity_target_count",
              "bioactivity_targets"),
    datasets=("gtopdb_ligand", "gtopdb_interaction", "gtopdb", "uniprot",
              "hgnc", "chembl_target"),
    chains=(">>gtopdb_ligand>>gtopdb_interaction>>gtopdb>>uniprot>>hgnc",
            ">>chembl_molecule>>chembl_target"),
    collect_fn=collect,
)
