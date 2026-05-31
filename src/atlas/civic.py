"""Shared CIViC clinical-evidence aggregation.

`civic_evidence` rows (reached via `>>hgnc>>civic_evidence` on the gene side or
`>>mondo>>civic_evidence` on the disease side) carry the drug × variant ×
indication triple:

    id | molecular_profile | disease | therapies | evidence_type | evidence_level | significance

This module collapses the PREDICTIVE subset — the drug-actionable
precision-medicine view ("therapy T for variant/subtype V in disease D") — to
unique (profile, therapy, indication, effect) associations, ranked by CIViC
evidence level (A validated → E inferential), keeping a representative evidence
id and a supporting-item count per association. Used by gene §10 (anchored via
hgnc) and disease §13 (anchored via mondo) so the two entity types share one
implementation of the join.
"""
from collections import Counter

# CIViC evidence levels (validated → inferential) and clinical-significance
# ordering — rank so the strongest, most-actionable rows surface first.
# Source: CIViC evidence-level / significance documentation.
LEVEL_RANK = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4}
SIG_RANK = {"Sensitivity/Response": 0, "Reduced Sensitivity": 1,
            "Resistance": 2, "Adverse Response": 3}


def clean_profile(p):
    """biobtree prefixes some fusion/expression molecular profiles with
    'v::' (e.g. 'v::ALK Fusion') — strip it for display."""
    p = (p or "").strip()
    return p[3:] if p.startswith("v::") else p


def therapy_label(t):
    """Combination therapies arrive comma-joined ('Cetuximab,Adagrasib');
    render as 'Cetuximab + Adagrasib'."""
    return " + ".join(s.strip() for s in (t or "").split(",") if s.strip())


def aggregate_predictive(rows, top=30):
    """civic_evidence target dicts -> (ranked associations[:top], stats).

    Dedupes the Predictive subset by (profile, therapy, indication, effect),
    keeping the best (lowest-rank) evidence level + a representative evidence
    id + a supporting-item count. `stats` carries evidence_type_counts,
    predictive_total, association_total, evidence_total for the header note.
    """
    type_counts = dict(Counter(r.get("evidence_type") for r in rows
                               if r.get("evidence_type")).most_common())
    predictive = [r for r in rows if r.get("evidence_type") == "Predictive"]
    assoc = {}
    for r in predictive:
        therapy = (r.get("therapies") or "").strip()
        profile = clean_profile(r.get("molecular_profile"))
        key = (profile, therapy, r.get("disease"), r.get("significance"))
        lvl = (r.get("evidence_level") or "").strip()
        cur = assoc.get(key)
        if cur is None:
            assoc[key] = {"profile": profile, "therapy": therapy,
                          "disease": r.get("disease"),
                          "significance": r.get("significance"),
                          "level": lvl, "evidence_id": r.get("id"), "n": 1}
        else:
            cur["n"] += 1
            if LEVEL_RANK.get(lvl, 9) < LEVEL_RANK.get(cur["level"], 9):
                cur["level"], cur["evidence_id"] = lvl, r.get("id")
    ranked = sorted(assoc.values(),
                    key=lambda d: (LEVEL_RANK.get(d["level"], 9),
                                   SIG_RANK.get(d["significance"], 4),
                                   -d["n"], d["profile"] or ""))
    stats = {"evidence_type_counts": type_counts,
             "predictive_total": len(predictive),
             "association_total": len(ranked),
             "evidence_total": len(rows)}
    return ranked[:top], stats
