"""§4 — indications. The drug's labelled/clinical disease uses, from the
chembl_molecule `indications` field (each {efo, mesh, highestDevelopmentPhase}),
deduped + cross-walked to MONDO on the anchor. Grouped by max phase; the MONDO
id + slug become /atlas/disease/ links once links.py exists."""
from atlas.section import Section


def collect(a):
    inds = [{
        "name": i.name, "efo_id": i.efo_id, "mesh_id": i.mesh_id,
        "mondo_id": i.mondo_id, "slug": i.slug, "max_phase": i.max_phase,
    } for i in a.indications]
    approved = [i for i in inds if i["max_phase"] >= 4]
    return {
        "section": "04_indications",
        "indication_count": len(inds),
        "approved_count": len(approved),
        "indications": inds,
    }


SECTION = Section(
    id="4", name="indications",
    description=("Drug indications (chembl_molecule.indications → efo/mesh, "
                 "cross-walked to MONDO), grouped by max development phase"),
    needs=("indications",),
    produces=("indication_count", "approved_count", "indications"),
    datasets=("chembl_molecule", "efo", "mesh", "mondo"),
    chains=(">>efo>>mondo", ">>mesh>>mondo"),
    collect_fn=collect,
)
