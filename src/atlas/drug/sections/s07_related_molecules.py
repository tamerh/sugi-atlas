"""§7 — related molecules. Competitor drugs that share ≥1 of this drug's
PRIMARY targets (the GtoPdb-curated set, not the 85-target bioactivity cloud —
fanning over the full set would flood the list with off-target overlap). For
each primary target's UniProt, pull phase-≥2 ChEMBL molecules hitting it,
aggregate by how many of the drug's targets each shares. The drug itself + its
salt forms are excluded. Competitor names become /atlas/drug/ links once
links.py exists."""
from atlas.biobtree import map_all
from atlas.section import Section


def _phase(v):
    try:
        return int(v or 0)
    except (TypeError, ValueError):
        return 0


def collect(a):
    self_ids = {a.chembl_id} | set(a.child_chembls)
    if a.parent_chembl:
        self_ids.add(a.parent_chembl)

    comp = {}  # chembl_id -> {id, name, phase, shared: set}
    for t in a.targets:
        if not t.uniprot:
            continue
        label = t.gene_symbol or t.target_name
        for m in map_all(t.uniprot,
                         ">>uniprot>>chembl_target>>chembl_molecule[highestDevelopmentPhase>=2]",
                         cap=5):
            mid = m.get("id")
            if not mid or mid in self_ids:
                continue
            e = comp.setdefault(mid, {"id": mid, "name": m.get("name"),
                                      "phase": _phase(m.get("highestDevelopmentPhase")),
                                      "shared": set()})
            e["shared"].add(label)

    ranked = sorted(comp.values(),
                    key=lambda d: (-len(d["shared"]), -d["phase"], d["name"] or ""))
    return {
        "section": "07_related_molecules",
        "competitor_count": len(comp),
        "related_molecules": [{"id": d["id"], "name": d["name"], "phase": d["phase"],
                               "shared_targets": sorted(d["shared"]),
                               "shared_count": len(d["shared"])}
                              for d in ranked[:30]],
    }


SECTION = Section(
    id="7", name="related_molecules",
    description=("Competitor drugs sharing ≥1 primary target (phase ≥2), "
                 "aggregated by shared-target count — fans over the curated "
                 "GtoPdb targets, not the full bioactivity set"),
    needs=("targets", "chembl_id", "parent_chembl", "child_chembls"),
    produces=("competitor_count", "related_molecules"),
    datasets=("uniprot", "chembl_target", "chembl_molecule"),
    chains=(">>uniprot>>chembl_target>>chembl_molecule[highestDevelopmentPhase>=2]",),
    collect_fn=collect,
)
