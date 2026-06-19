"""§14 — DrugCentral regulatory approval + curated mechanism.

DrugCentral's authoritative signal that ChEMBL/GtoPdb don't carry cleanly:
multi-agency regulatory approval (FDA / EMA / PMDA) and the curated
mechanism-of-action target set + action types.

Reached by the direct chembl_molecule→drugcentral edge (restored in the biobtree
reindex that fixed BIOBTREE_ISSUES #57; we previously used an InChIKey-search
workaround). The edge row carries the approval flags + target_count directly
(pmda now explicit true/false, not absent); one entry fetch adds moa_targets +
action_types. Small molecules resolve 1:1; biologics return nothing (DrugCentral
is small-molecule — correct, not a miss).

We deliberately do NOT re-list DrugCentral's full target set — the drug's targets
are already resolved + linked in §2/§3/§8. This section is the regulatory +
MOA layer those sections lack."""
from atlas.biobtree import map_all, entry
from atlas.section import Section

_AGENCIES = (("fda_approved", "FDA"), ("ema_approved", "EMA"), ("pmda_approved", "PMDA"))


def _int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def _attrs(e):
    a = (e.get("Attributes") or {})
    return a.get("Drugcentral") or a.get("DrugCentral") or {}


def _approved(v):
    return v in (True, "true")


def collect(a):
    hits = map_all(a.chembl_id, ">>chembl_molecule>>drugcentral")
    hit = next((r for r in hits if (r.get("id") or "").strip()), None)
    if not hit:
        return {"section": "14_drugcentral", "found": False}

    # Approval flags + target_count are on the edge row (pmda now explicit, so an
    # un-flagged agency is a real "not approved", not missing data). moa_targets /
    # action_types still come from the entry payload.
    dc = _attrs(entry(hit["id"], "drugcentral"))
    approvals = [label for key, label in _AGENCIES
                 if _approved(hit.get(key)) or _approved(dc.get(key))]
    return {
        "section": "14_drugcentral",
        "found": True,
        "drugcentral_id": hit["id"],
        "drugcentral_name": (hit.get("name") or dc.get("name") or "").strip(),
        "approvals": approvals,
        "action_types": [t for t in (dc.get("action_types") or []) if t],
        "moa_targets": [t for t in (dc.get("moa_targets") or []) if t],
        "target_count": _int(hit.get("target_count") or dc.get("target_count")),
    }


SECTION = Section(
    id="14", name="drugcentral",
    description=("DrugCentral regulatory approval (FDA / EMA / PMDA) + curated "
                 "mechanism (action types, MOA targets) via the direct "
                 "chembl_molecule→drugcentral edge."),
    needs=("chembl_id",),
    produces=("found", "drugcentral_id", "drugcentral_name", "approvals",
              "action_types", "moa_targets", "target_count"),
    datasets=("chembl_molecule", "drugcentral"),
    chains=(">>chembl_molecule>>drugcentral",),
    collect_fn=collect,
)
