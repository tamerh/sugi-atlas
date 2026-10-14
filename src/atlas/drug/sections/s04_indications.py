"""§4 — indications. The drug's labelled/clinical disease uses, from the
chembl_molecule `indications` field (each {efo, mesh, highestDevelopmentPhase}),
deduped + cross-walked to MONDO on the anchor. Grouped by max phase; the MONDO
id + slug become /atlas/disease/ links once links.py exists."""
from atlas.section import Section
from atlas.indication import approved_indication, late_stage_use


def collect(a):
    # `approved` is computed ONCE here (the anchor has both ATC + the indication
    # list) and read downstream by r_indications, the mesh, and the disease-side
    # indication index — phase 4, or an anticancer drug at phase 3 vs a cancer
    # (ChEMBL logs imatinib→CML at phase 3 though FDA-approved). See atlas.indication.
    atc = list(a.atc_codes)
    # Molecule-level FDA approval (PubChem is_fda_approved, or ChEMBL max_phase 4)
    # gates the anticancer phase-3 upgrade — see atlas.indication.
    mol_approved = bool(a.is_fda_approved) or (a.max_phase or 0) >= 4
    inds = [{
        "name": i.name, "efo_id": i.efo_id, "mesh_id": i.mesh_id,
        "mondo_id": i.mondo_id, "slug": i.slug, "max_phase": i.max_phase,
        # Approved = ChEMBL phase 4 OR a DrugCentral on-label use (biobtree v2.17).
        "approved": approved_indication(atc, i.max_phase, i.name, mol_approved) or bool(i.label),
        "label": bool(i.label),
        "late": late_stage_use(i.max_phase, mol_approved) and not i.label,
    } for i in a.indications]
    approved = [i for i in inds if i["approved"]]
    return {
        "section": "04_indications",
        "indication_count": len(inds),
        "approved_count": len(approved),
        "indications": inds,
        # DrugCentral's on-label indications as labelled (free text; those with an
        # exact Mondo term are also rows above, flagged `label`).
        "label_indications": list(getattr(a, "label_indications", ()) or ()),
    }


SECTION = Section(
    id="4", name="indications",
    description=("Drug indications (chembl_molecule.indications → efo/mesh, "
                 "cross-walked to MONDO), grouped by max development phase"),
    needs=("indications",),
    produces=("indication_count", "approved_count", "indications", "label_indications"),  # indications carry `approved`
    datasets=("chembl_molecule", "efo", "mesh", "mondo", "drugcentral"),
    chains=(">>efo>>mondo", ">>mesh>>mondo", ">>chembl_molecule>>drugcentral>>mondo"),
    collect_fn=collect,
)
