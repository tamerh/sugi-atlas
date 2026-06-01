"""§7 — related molecules. Competitor / same-mechanism molecules that share
≥1 of this drug's PRIMARY targets (the GtoPdb-curated set, not the full
bioactivity cloud — fanning over that floods the list with off-target overlap).

Two complementary biobtree sources, merged + deduplicated by drug name:
  - ChEMBL   — `>>uniprot>>chembl_target>>chembl_molecule[highestDevelopmentPhase>=2]`
               (clinical-stage candidates; carries development phase)
  - PubChem  — `>>uniprot>>pubchem_activity>>pubchem`, filtered to
               `compound_type == "drug"` (approved / known drugs acting on the
               target; broadens coverage beyond ChEMBL, carries FDA flag)

Each molecule is aggregated by how many of the drug's primary targets it shares.
The drug itself + its salt forms are excluded. Competitor names become
/atlas/drug/ links once links.py exists."""
from atlas.biobtree import map_all
from atlas.section import Section

# Trailing salt / hydrate words to strip when normalising a name for dedup
# (so "regorafenib anhydrous" == "regorafenib", "imatinib mesylate" == self).
_SALT_SUFFIXES = (
    " anhydrous", " mesylate", " hydrochloride", " dihydrochloride",
    " sulfate", " sulphate", " maleate", " citrate", " tartrate", " fumarate",
    " sodium", " potassium", " calcium", " phosphate", " acetate", " besylate",
    " hydrobromide", " succinate", " hemifumarate",
)


def _phase(v):
    try:
        return int(v or 0)
    except (TypeError, ValueError):
        return 0


def _fda(v):
    return v in (True, "true", "True", 1, "1")


def _strip_salt(name, lower):
    """Strip trailing salt/hydrate words. `lower` controls case folding —
    True for the dedup key, False to keep display casing."""
    n = (name or "").strip()
    if lower:
        n = n.lower()
    changed = True
    while changed:
        changed = False
        for suf in _SALT_SUFFIXES:
            if n.lower().endswith(suf):
                n = n[: -len(suf)].rstrip()
                changed = True
    return n


def _norm(name):
    return _strip_salt(name, lower=True)


def _clean_display(name):
    return _strip_salt(name, lower=False)


def collect(a):
    self_ids = {a.chembl_id} | set(a.child_chembls)
    if a.parent_chembl:
        self_ids.add(a.parent_chembl)
    self_key = _norm(a.canonical_name)

    comp = {}  # key (normalised name | id) -> record

    def _add(key, name, phase, source, fda, label):
        e = comp.setdefault(key, {"name": name, "phase": phase,
                                  "sources": set(), "shared": set(), "fda": fda})
        e["sources"].add(source)
        e["shared"].add(label)
        # Prefer the cleaner display name (shorter wins → drops salt/hydrate
        # suffixes like "… ANHYDROUS" when a bare form is also seen).
        if name and (not e["name"] or len(name) < len(e["name"])):
            e["name"] = name
        if phase > e["phase"]:
            e["phase"] = phase
        e["fda"] = e["fda"] or fda

    for t in a.targets:
        if not t.uniprot:
            continue
        label = t.gene_symbol or t.target_name

        # ChEMBL clinical-stage molecules (phase ≥2).
        for m in map_all(t.uniprot,
                         ">>uniprot>>chembl_target>>chembl_molecule[highestDevelopmentPhase>=2]",
                         cap=5):
            mid = m.get("id")
            if not mid or mid in self_ids:
                continue
            name = m.get("name")
            key = _norm(name) or mid
            if key == self_key:
                continue
            _add(key, name, _phase(m.get("highestDevelopmentPhase")), "ChEMBL", False, label)

        # PubChem drug-class molecules acting on the same target. The route
        # returns ~600 mixed compounds per target; we keep only compound_type
        # == "drug" (named, mostly approved) to avoid flooding with the
        # literature/bioassay cloud.
        for m in map_all(t.uniprot, ">>uniprot>>pubchem_activity>>pubchem", cap=5):
            if m.get("compound_type") != "drug":
                continue
            name = (m.get("title") or "").strip()
            if not name:
                continue
            key = _norm(name)
            if key == self_key:
                continue
            _add(key, name, 0, "PubChem", _fda(m.get("is_fda_approved")), label)

    ranked = sorted(comp.values(),
                    key=lambda d: (-len(d["shared"]), -d["phase"],
                                   0 if d["fda"] else 1, (d["name"] or "").lower()))
    return {
        "section": "07_related_molecules",
        "competitor_count": len(comp),
        "related_molecules": [{"name": _clean_display(d["name"]), "phase": d["phase"],
                               "fda": d["fda"],
                               "sources": sorted(d["sources"]),
                               "shared_targets": sorted(d["shared"]),
                               "shared_count": len(d["shared"])}
                              for d in ranked[:60]],
    }


SECTION = Section(
    id="7", name="related_molecules",
    description=("Competitor / same-mechanism molecules sharing ≥1 curated "
                 "primary target, merged from ChEMBL (phase ≥2) and PubChem "
                 "(drug-class bioactivity), aggregated by shared-target count"),
    needs=("targets", "chembl_id", "canonical_name", "parent_chembl", "child_chembls"),
    produces=("competitor_count", "related_molecules"),
    datasets=("uniprot", "chembl_target", "chembl_molecule",
              "pubchem_activity", "pubchem"),
    chains=(">>uniprot>>chembl_target>>chembl_molecule[highestDevelopmentPhase>=2]",
            ">>uniprot>>pubchem_activity>>pubchem"),
    collect_fn=collect,
)
