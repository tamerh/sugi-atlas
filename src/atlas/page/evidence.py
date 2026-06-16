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
    # Pathway has a FLAT bundle (no section ids) — components() reads it directly;
    # the sid/field slots are unused, only the (key, weight) pair matters.
    "pathway": [
        ("member_count", "", "member_count", 1.0),
        ("child_count",  "", "child_count",  0.5),
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
    if entity_type == "pathway":
        # Flat bundle: read the fields directly off it (no section nesting).
        return {key: _int(bundle.get(field))
                for key, _sid, field, _w in _SPEC.get("pathway", [])}
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


# Per-component corpus-relative ranking ("top 1% of genes by …"). Only the
# components worth a headline rank get an entry → (plural noun, min absolute
# count). The floor matters: most entities score 0 on most components, so a
# "top 2%" on a count of 3 would be noise — require a real count first.
_RANK_SPEC = {
    ("gene", "drug_count"):        ("genes", 10),
    ("gene", "gwas_count"):        ("genes", 25),
    ("gene", "variant_count"):     ("genes", 50),
    ("gene", "interaction_count"): ("genes", 50),
    ("disease", "trial_count"):    ("diseases", 25),
    ("disease", "gwas_count"):     ("diseases", 25),
    ("disease", "gene_count"):     ("diseases", 25),
    ("drug", "target_count"):      ("drugs", 10),
    ("drug", "indication_count"):  ("drugs", 5),
}
_RANK_FLOOR_PCT = 90              # only surface a rank in the top decile


# ---- render-side: load the frozen distribution, look an entity up -----------
_PROM = {}    # {entity_type: {slug: score}}
_DIST = {}    # {entity_type: [sorted raw signals]} — for single-entity reuse
_DIST_COMP = {}   # {entity_type: {component_key: [sorted counts]}} — per-metric rank
_COMP = {}    # {entity_type: {slug: {component_key: count}}} — per-slug lookup
_ROOT = None


def load(dist_root):
    """Load <dist>/atlas/evidence.json (idempotent per root)."""
    global _PROM, _DIST, _DIST_COMP, _COMP, _ROOT
    if dist_root == _ROOT:
        return
    _PROM, _DIST, _DIST_COMP, _COMP, _ROOT = {}, {}, {}, {}, dist_root
    path = os.path.join(dist_root, "atlas", "evidence.json")
    if os.path.exists(path):
        d = json.load(open(path))
        _DIST = d.pop("_dist", {})
        _DIST_COMP = d.pop("_dist_components", {})
        _COMP = d.pop("_components", {})
        _PROM = d


def reset():
    global _PROM, _DIST, _DIST_COMP, _COMP, _ROOT
    _PROM, _DIST, _DIST_COMP, _COMP, _ROOT = {}, {}, {}, {}, None


def components_for(entity_type, slug):
    """The frozen per-component counts for another built entity (by slug) — lets
    a render read a *related* page's headline numbers (a thin disease quoting its
    parent's evidence, a lncRNA quoting its sense gene). {} when unknown."""
    return _COMP.get(entity_type, {}).get(slug) or {}


def illumination(comps):
    """An Atlas-derived 'illumination level' for a gene from its evidence
    components: drug-illuminated > chemically-probed > functionally-characterized
    > understudied. Deterministic, self-contained, with the rule inline. NOT a
    Pharos TDL clone — the corpus has effectively no truly-dark genes, so we frame
    it as illumination and never borrow the "Tdark" label."""
    def g(k):
        return int((comps or {}).get(k) or 0)
    if g("drug_count") and g("civic_count"):
        return "drug-illuminated"          # has a drug AND curated precision evidence
    if g("drug_count"):
        return "chemically-probed"         # ChEMBL molecules but no curated actionability
    if g("gwas_count") or g("variant_count") or g("interaction_count"):
        return "functionally-characterized"
    return "understudied"


def component_percentile(entity_type, key, count):
    """Percentile rank (0-100) of `count` within the frozen per-component
    distribution for this entity type. None when no distribution is loaded
    (no batch has run — e.g. a unit test calling a renderer directly)."""
    dist = _DIST_COMP.get(entity_type, {}).get(key)
    if not dist:
        return None
    n = len(dist)
    if n <= 1:
        return 100
    return round(100 * bisect.bisect_left(dist, count) / (n - 1))


def rank_clause(entity_type, key, count):
    """A corpus-relative ranking parenthetical for a headline count —
    ' (top 1% of genes corpus-wide)' — or '' when the metric isn't rank-worthy,
    the count is below its floor, the rank isn't in the top decile, or no
    distribution is loaded. Honest by construction: it always names the
    reference set (the Atlas corpus), never implies completeness of biology."""
    spec = _RANK_SPEC.get((entity_type, key))
    if not spec:
        return ""
    noun, floor = spec
    if (count or 0) < floor:
        return ""
    pct = component_percentile(entity_type, key, count)
    if pct is None or pct < _RANK_FLOOR_PCT:
        return ""
    gap = 100 - pct
    band = "top 1%" if gap < 1 else f"top {round(gap)}%"
    return f" ({band} of {noun} corpus-wide)"


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
