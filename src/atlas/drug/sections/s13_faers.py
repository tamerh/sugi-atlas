"""§13 — post-marketing adverse events (openFDA FAERS, CC0).

The FDA Adverse Event Reporting System: spontaneous post-marketing reports of
adverse events co-occurring with a drug. biobtree exposes it as master (per
drug name: total_reports, distinct_reactions) + faers_reaction child (per
drug×reaction: MedDRA preferred term, report_count, PRR disproportionality):

    chembl_molecule >> faers >> faers_reaction

A drug resolves to SEVERAL FAERS name records (generic + brand + salt, e.g.
imatinib / glivec / imatinib mesylate), each with its own counts — so reactions
are aggregated by preferred term across those records (report_count summed; PRR
pooled report-count-weighted).

Two read-outs, because sort order changes the story:
  - most-reported  — by volume (the events clinicians see; dominated by common,
    non-specific terms like nausea / drug ineffective);
  - disproportionate — by PRR (the safety SIGNAL: events reported more for this
    drug than the FAERS background), gated by a report-count floor so a PRR
    spike off a handful of reports isn't mistaken for signal.

CRITICAL framing: FAERS is co-occurrence, NOT causation — reports are
unverified, voluntary, and confounded by indication/notoriety. PRR is a
disproportionality heuristic, not a risk estimate.

`serious_count` is now in the lite projection (biobtree restored it) and is
summed per PT + rendered as a "Serious" column. `outcome` is also projected but
deliberately NOT surfaced: it's a single ICH reaction-outcome code per PT
(1=recovered … 5=fatal, 6=unknown, dominated by 6) and biobtree's per-PT
aggregation (modal vs most-severe) is unconfirmed — see BIOBTREE_ISSUES."""
from atlas.biobtree import map_all
from atlas.section import Section

_MASTER_CHAIN = ">>chembl_molecule>>faers"
_REACTION_CHAIN = ">>chembl_molecule>>faers>>faers_reaction"

# Report-count floor for the disproportionality (PRR) view — a high PRR off a
# few reports is noise, not signal (same lesson as under-powered ORA).
_PRR_MIN_REPORTS = 10
# How many rows each view keeps in the bundle (render caps again at ROW_CAP).
_TOP_N = 50
# map_all `cap` bounds PAGES (~100 rows each). FAERS reaction lists are large
# (paracetamol ~25.8k rows) and are NOT globally sorted by report_count, so the
# default cap=60 (~6.1k) silently truncated mid-way and corrupted the per-PT
# aggregation across name records. biobtree already aggregates (drops count=1
# singletons), so the full set is bounded + cheap to pull (sub-second on
# localhost); 500 pages (~50k rows) clears the biggest drugs with headroom.
_FETCH_CAP = 500


def _int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def _float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _aggregate_reactions(rows, prr_min_reports=_PRR_MIN_REPORTS, top_n=_TOP_N):
    """Group faers_reaction rows by MedDRA preferred term across the drug's FAERS
    name records: report_count summed, PRR pooled report-count-weighted (a single
    representative disproportionality dominated by the largest record). Returns
    (most_reported, disproportionate, distinct_count) — two views (by volume; by
    PRR among reactions with >= prr_min_reports reports) + the distinct PT count.
    Pure + deterministic (ties broken by reaction name)."""
    agg: dict = {}
    for r in rows:
        pt = (r.get("reaction") or "").strip()
        rc = _int(r.get("report_count")) or 0
        prr = _float(r.get("prr"))
        if not pt or rc <= 0:
            continue
        g = agg.setdefault(pt, {"reaction": pt, "report_count": 0, "serious_count": 0,
                                "_prr_num": 0.0, "_prr_den": 0, "records": 0})
        g["report_count"] += rc
        g["serious_count"] += _int(r.get("serious_count")) or 0
        g["records"] += 1
        if prr is not None:
            g["_prr_num"] += prr * rc
            g["_prr_den"] += rc

    reactions = []
    for g in agg.values():
        prr = round(g["_prr_num"] / g["_prr_den"], 2) if g["_prr_den"] else None
        reactions.append({"reaction": g["reaction"], "report_count": g["report_count"],
                          "serious_count": g["serious_count"], "prr": prr,
                          "records": g["records"]})

    most_reported = sorted(reactions, key=lambda x: (-x["report_count"], x["reaction"]))[:top_n]
    disproportionate = sorted(
        [x for x in reactions if x["prr"] is not None and x["report_count"] >= prr_min_reports],
        key=lambda x: (-(x["prr"] or 0.0), -x["report_count"], x["reaction"]))[:top_n]
    return most_reported, disproportionate, len(reactions)


def collect(a):
    master = map_all(a.chembl_id, _MASTER_CHAIN)
    name_records = [{"drug_name": m.get("drug_name"),
                     "total_reports": _int(m.get("total_reports")),
                     "distinct_reactions": _int(m.get("distinct_reactions"))}
                    for m in master if m.get("drug_name")]
    total_reports = sum(r["total_reports"] or 0 for r in name_records)

    most_reported, disproportionate, distinct_reactions = _aggregate_reactions(
        map_all(a.chembl_id, _REACTION_CHAIN, cap=_FETCH_CAP))

    return {
        "section": "13_faers",
        "chembl_id": a.chembl_id,
        "total_reports": total_reports,
        "distinct_reactions": distinct_reactions,
        "name_records": name_records,
        "name_record_count": len(name_records),
        "most_reported": most_reported,
        "disproportionate": disproportionate,
        "prr_min_reports": _PRR_MIN_REPORTS,
    }


SECTION = Section(
    id="13", name="faers",
    description=("Post-marketing adverse events (openFDA FAERS): MedDRA reactions "
                 "with report counts + PRR disproportionality, aggregated across "
                 "the drug's FAERS name records via chembl_molecule→faers→"
                 "faers_reaction. Co-occurrence, not causation."),
    needs=("chembl_id",),
    produces=("total_reports", "distinct_reactions", "name_records",
              "name_record_count", "most_reported", "disproportionate",
              "prr_min_reports"),
    datasets=("chembl_molecule", "faers", "faers_reaction"),
    chains=(_MASTER_CHAIN, _REACTION_CHAIN),
    collect_fn=collect,
)
