"""Per-entity evidence_score score for search ranking + featured-card curation.

The atlas already computes, during mining, every signal that indicates how
well-studied an entity is (curated clinical evidence, drugs, variants, trials,
…) — they're just buried in body text. This module surfaces them:

  - `evidence_components` — the raw integer counts, with STABLE keys per
    entity type (every gene has `variant_count`, etc.; 0 when absent), so the
    web layer can re-weight client-side without a backend change.
  - `evidence_score` — a single 0-100 score, the entity's PERCENTILE RANK within
    its own type. Per-type normalization keeps a cross-type result list fair
    (a gene's CIViC count and a disease's GWAS count aren't the same scale);
    percentile is robust to the heavy skew of biomedical counts. It is NOT
    pre-damped — the consumer applies `log(1+evidence_score)` as a relevance
    multiplier so evidence_score boosts but never drowns an exact-name match.

The percentile needs the corpus-wide distribution, so it's computed in the
batch MERGE phase and written to `<dist>/atlas/evidence.json`; a single-entity
rebuild reuses that frozen distribution (the `_dist` block) so its rank doesn't
drift. The raw composite that drives the rank weights curated/clinical evidence
above raw counts and damps inflated ones (trials, patents), mirroring the
disease corpus `signal_score`.
"""
import bisect
import json
import math
import os

# (component_key, section_id, bundle_field, weight) per entity type.
# Component keys are the STABLE frontmatter schema; fields are real bundle
# counts (verified against built bundles). Weights drive the composite raw
# signal only — curated/clinical high, broad/inflated counts damped.
_SPEC = {
    "gene": [
        ("civic_count",       "10", "civic_evidence_total", 2.0),
        ("drug_count",        "10", "molecule_count",        1.5),
        ("variant_count",     "6",  "clinvar_total",         1.0),
        ("trial_count",       "10", "disease_trial_count",   1.0),
        ("gwas_count",        "12", "gwas_total",            0.7),
        ("interaction_count", "8",  "string_count",          0.5),
    ],
    "disease": [
        ("civic_count",   "13", "civic_evidence_total", 2.0),
        ("gene_count",    "5",  "gene_count",           1.0),
        ("variant_count", "3",  "clinvar_total",        1.0),
        ("drug_count",    "10", "phased_count",         1.0),
        ("gwas_count",    "2",  "assoc_total",          0.7),
        ("trial_count",   "13", "trial_count",          0.4),
        # Orphanet/HPO phenotype count. Weight 0.0 → reported in the frontmatter
        # (so the corpus-depth stats can tell a symptom-characterized disease from
        # a truly empty one) but it does NOT enter the evidence_score, leaving the
        # ranking unchanged.
        ("phenotype_count", "1", "phenotype_count",     0.0),
    ],
    "drug": [
        ("civic_count",         "10", "civic_evidence_total",   2.0),
        ("indication_count",    "4",  "indication_count",       1.5),
        ("trial_count",         "5",  "trial_count",            1.0),
        ("target_count",        "2",  "bioactivity_target_count", 0.7),
        ("activity_count",      "3",  "activity_total",         0.5),
        ("patent_family_count", "11", "patent_family_total",    0.3),
    ],
}


def _int(v):
    try:
        return int(v) if v is not None else 0
    except (TypeError, ValueError):
        return 0


def components(entity_type, bundle):
    """Raw integer component counts for an entity, with the fixed per-type keys
    always present (0 when the section/field is absent)."""
    out = {}
    for key, sid, field, _w in _SPEC.get(entity_type, []):
        sec = bundle.get(sid) if isinstance(bundle.get(sid), dict) else {}
        out[key] = _int((sec or {}).get(field))
    return out


def raw_signal(entity_type, comps):
    """Weighted sum of log1p(count) — the composite that drives the rank. log1p
    keeps a single huge count (e.g. STRING partners) from dominating; the
    percentile makes the absolute scale irrelevant, only the ordering matters."""
    by_key = {k: w for k, _s, _f, w in _SPEC.get(entity_type, [])}
    return sum(w * math.log1p(comps.get(k, 0)) for k, w in by_key.items())


def percentiles(raw_by_slug):
    """slug -> 0-100 percentile rank of its raw signal (within this set). Ties
    take the lower rank; the top entity is 100, the bottom 0."""
    vals = sorted(raw_by_slug.values())
    n = len(vals)
    if n <= 1:
        return {s: 100 for s in raw_by_slug}
    return {s: round(100 * bisect.bisect_left(vals, v) / (n - 1))
            for s, v in raw_by_slug.items()}


# ---- render-side: load the frozen distribution, look an entity up -----------
_PROM = {}    # {entity_type: {slug: score}}
_DIST = {}    # {entity_type: [sorted raw signals]} — for single-entity reuse
_ROOT = None


def load(dist_root):
    """Load <dist>/atlas/evidence.json (idempotent per root)."""
    global _PROM, _DIST, _ROOT
    if dist_root == _ROOT:
        return
    _PROM, _DIST, _ROOT = {}, {}, dist_root
    path = os.path.join(dist_root, "atlas", "evidence.json")
    if os.path.exists(path):
        d = json.load(open(path))
        _DIST = d.pop("_dist", {})
        _PROM = d


def reset():
    global _PROM, _DIST, _ROOT
    _PROM, _DIST, _ROOT = {}, {}, None


def lookup(entity_type, slug, raw):
    """The entity's evidence_score: its batch-computed score if known, else its
    percentile against the last batch's frozen distribution (single-entity
    rebuild), else None (no batch has run — emit components only)."""
    if entity_type in _PROM and slug in _PROM[entity_type]:
        return _PROM[entity_type][slug]
    dist = _DIST.get(entity_type)
    if dist:
        n = len(dist)
        return round(100 * bisect.bisect_left(dist, raw) / (n - 1)) if n > 1 else 100
    return None
