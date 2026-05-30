"""§2 — GWAS landscape: mondo → GWAS assocs + studies. Counts, top hits,
study-level metadata.

Deterministic collector. Strategy:
  - assoc/study totals come straight off DiseaseAnchors.xref_counts (free).
  - top assocs: map_all mondo→gwas (capped ~1000 rows), sort by p_value asc,
    slice top 50, then one entry() per top assoc to pull rsid / risk_allele /
    odds_ratio / mapped hgnc ids that aren't in the list view.
  - studies: map_all mondo→gwas_study (capped), keep top 10 by case count
    (derived from initial_sample_size on entry()).
"""
from collections import OrderedDict
import re

from atlas.section import Section
from atlas.biobtree import map_all, entry

CHAINS   = (">>mondo>>gwas", ">>mondo>>gwas_study")
DATASETS = ("mondo", "gwas", "gwas_study")

# Bounds — keep biobtree calls predictable.
_GWAS_PAGE_CAP   = 10   # × 100 rows = up to ~1000 assocs scanned
_STUDY_PAGE_CAP  = 5    # × 100 rows = up to ~500 studies scanned
_TOP_ASSOCS_N    = 50
_TOP_STUDIES_N   = 10
_ENRICH_STUDIES_N = 50   # how many list-view studies to entry()-enrich before ranking

_CASES_RE = re.compile(r"([\d,]+)\s+[^,]*?cases", re.IGNORECASE)


def _parse_pvalue(s):
    """GWAS p_value comes as a string like '4.000000e-08'; missing/garbage → inf."""
    if s is None:
        return float("inf")
    try:
        v = float(s)
        return v if v > 0 else float("inf")
    except (TypeError, ValueError):
        return float("inf")


def _parse_cases(initial_sample_size):
    """Best-effort cases count for ranking studies; 0 if no match."""
    if not initial_sample_size:
        return 0
    m = _CASES_RE.search(initial_sample_size)
    if not m:
        return 0
    try:
        return int(m.group(1).replace(",", ""))
    except ValueError:
        return 0


def _gwas_attrs(eid):
    e = entry(eid, "gwas") or {}
    return ((e.get("Attributes") or {}).get("Gwas")) or {}


def _study_attrs(sid):
    e = entry(sid, "gwas_study") or {}
    return ((e.get("Attributes") or {}).get("GwasStudy")) or {}


def collect(a):
    xc = dict(a.xref_counts or {})
    assoc_total = int(xc.get("gwas", 0))
    study_total = int(xc.get("gwas_study", 0))

    # ---- top assocs ------------------------------------------------------
    rows = map_all(a.mondo_id, ">>mondo>>gwas", cap=_GWAS_PAGE_CAP)
    rows_sorted = sorted(rows, key=lambda r: _parse_pvalue(r.get("p_value")))
    # Walk the sorted list, fetch each entry, and dedupe by rsID as we go
    # (best p-value wins because the list is already p-sorted). Stops once
    # we have _TOP_ASSOCS_N unique rsIDs — bounds entry cost the same as the
    # old slice-then-fetch path but produces 50 unique variants instead of
    # the prior 50 repeated rows for the same handful of lead SNPs.
    seen_rsids = set()
    top_rows = []
    top_attrs = []  # parallel array: cached gwas entries to avoid re-fetching
    for r in rows_sorted:
        gid = r.get("id")
        if not gid:
            continue
        attrs = _gwas_attrs(gid)
        rsid = attrs.get("snp_id") or r.get("snp_id")
        if not rsid:
            # No rsID — keep at most one such row to retain coverage
            rsid = f"_no_rsid:{gid}"
        if rsid in seen_rsids:
            continue
        seen_rsids.add(rsid)
        top_rows.append(r)
        top_attrs.append(attrs)
        if len(top_rows) >= _TOP_ASSOCS_N:
            break

    top_assocs = []
    gene_ids = set()
    for r, attrs in zip(top_rows, top_attrs):
        gid = r.get("id")
        rsid = attrs.get("snp_id") or None
        # strongest_snp_risk_allele looks like 'rs380390-C'; split off the allele
        risk_allele = None
        srsa = attrs.get("strongest_snp_risk_allele")
        if srsa and "-" in srsa:
            risk_allele = srsa.rsplit("-", 1)[1] or None
        or_beta = attrs.get("or_beta")
        try:
            odds_ratio = float(or_beta) if or_beta not in (None, "", "NR") else None
        except (TypeError, ValueError):
            odds_ratio = None
        gene_symbol = attrs.get("mapped_gene") or r.get("mapped_gene") or None
        for hg in attrs.get("snp_gene_ids") or []:
            gene_ids.add(hg)
        top_assocs.append({
            "id": gid,
            "rsid": rsid,
            "pvalue": _parse_pvalue(attrs.get("p_value") or r.get("p_value")),
            "gene_symbol": gene_symbol,
            "risk_allele": risk_allele,
            "odds_ratio": odds_ratio,
        })

    # ---- studies ---------------------------------------------------------
    study_rows = map_all(a.mondo_id, ">>mondo>>gwas_study", cap=_STUDY_PAGE_CAP)
    # de-dup by id, preserve order
    seen = OrderedDict()
    for s in study_rows:
        sid = s.get("id")
        if sid and sid not in seen:
            seen[sid] = s

    enriched = []
    for sid, s in list(seen.items())[:_ENRICH_STUDIES_N]:
        attrs = _study_attrs(sid)
        iss = attrs.get("initial_sample_size") or ""
        cases = _parse_cases(iss)
        # controls = anything between cases and end that mentions 'control'
        controls = 0
        m = re.search(r"([\d,]+)\s+[^,]*?controls", iss, re.IGNORECASE)
        if m:
            try:
                controls = int(m.group(1).replace(",", ""))
            except ValueError:
                controls = 0
        year = None
        pubdate = attrs.get("publication_date") or s.get("publication_date")
        if pubdate and len(pubdate) >= 4 and pubdate[:4].isdigit():
            year = int(pubdate[:4])
        enriched.append({
            "id": sid,
            "lead_author": attrs.get("first_author") or s.get("first_author"),
            "year": year,
            "sample_size_cases": cases,
            "sample_size_controls": controls,
            "title": attrs.get("study"),
        })
    # rank by case count desc, then year desc; fall back to anchor's ordering
    enriched.sort(key=lambda d: (d["sample_size_cases"], d["year"] or 0), reverse=True)
    studies = enriched[:_TOP_STUDIES_N]

    return {
        "section": "02_gwas_landscape",
        "mondo_id": a.mondo_id,
        "assoc_total": assoc_total,
        "study_total": study_total,
        "top_assocs": top_assocs,
        "studies": studies,
        "unique_gene_count": len(gene_ids),
    }


SECTION = Section(
    id="2", name="gwas_landscape",
    description=("GWAS associations + studies for the disease — total counts, "
                 "top assocs by p-value, study-level meta (lead author, year, "
                 "case/control counts)."),
    needs=("mondo_id", "xref_counts"),
    produces=("assoc_total", "study_total", "top_assocs", "studies",
              "unique_gene_count"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
