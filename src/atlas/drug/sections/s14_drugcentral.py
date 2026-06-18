"""§14 — DrugCentral regulatory approval + curated mechanism.

DrugCentral's authoritative signal that ChEMBL/GtoPdb don't carry cleanly:
multi-agency regulatory approval (FDA / EMA / PMDA) and the curated
mechanism-of-action target set + action types.

There is no chembl_molecule→drugcentral edge in this biobtree build (db_v38 was
built before the InChIKey-lookup edge fix; the clean edge lands in db_v39 — see
BIOBTREE_ISSUES #57). Until then we reach DrugCentral by its INCHIKEY — a
structural key, NOT fragile name matching — in exactly two calls: a
source-scoped search on the drug's InChIKey, then one entry fetch. Small
molecules resolve 1:1; biologics (no InChIKey) cleanly return nothing (DrugCentral
is small-molecule, so that's correct, not a miss).

We deliberately do NOT re-list DrugCentral's full target set — the drug's targets
are already resolved + linked in §2/§3/§8. This section is the regulatory +
MOA layer those sections lack."""
from atlas.biobtree import search, rows, entry
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


def collect(a):
    ik = getattr(a, "inchi_key", None)
    if not ik:
        return {"section": "14_drugcentral", "found": False}
    hits = rows(search(ik, source="drugcentral"))
    hit = next((r for r in hits if (r.get("id") or "").strip()), None)
    if not hit:
        return {"section": "14_drugcentral", "found": False}

    dc = _attrs(entry(hit["id"], "drugcentral"))
    if not dc:
        return {"section": "14_drugcentral", "found": False}

    # DrugCentral booleans are present-and-true or absent; absent ≠ "rejected",
    # so we only assert the agencies it explicitly flags approved.
    approvals = [label for key, label in _AGENCIES if dc.get(key) in (True, "true")]
    return {
        "section": "14_drugcentral",
        "found": True,
        "drugcentral_id": hit["id"],
        "drugcentral_name": (dc.get("name") or "").strip(),
        "approvals": approvals,
        "action_types": [t for t in (dc.get("action_types") or []) if t],
        "moa_targets": [t for t in (dc.get("moa_targets") or []) if t],
        "target_count": _int(dc.get("target_count")),
    }


SECTION = Section(
    id="14", name="drugcentral",
    description=("DrugCentral regulatory approval (FDA / EMA / PMDA) + curated "
                 "mechanism (action types, MOA targets), reached by InChIKey "
                 "lookup (db_v38 lacks the chembl_molecule→drugcentral edge)."),
    needs=("chembl_id",),
    produces=("found", "drugcentral_id", "drugcentral_name", "approvals",
              "action_types", "moa_targets", "target_count"),
    datasets=("drugcentral",),
    chains=(),
    collect_fn=collect,
)
