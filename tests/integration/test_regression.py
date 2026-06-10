"""Regression canaries the structural suite can't see.

The rest of the integration suite proves every page is *well-formed*; it has no
notion of how many pages should exist, nor that specific features are present.
Those are the two gaps a silent regression slips through — a collector that
drops a whole entity type still yields valid pages, and a refactor that quietly
removes the CIViC legend or the variant links leaves the markdown valid. These
checks close that:

  1. coverage floors  — catch catastrophic / large page loss
  2. golden pages     — pin feature *presence/shape* on a few canonical pages
  3. aggregate ratios — assert a feature holds across the whole corpus

All assertions are SHAPE-based, never exact counts — a biobtree data refresh
shifts values under the same code, so "CYP2D6 has a CPIC guideline" is durable
where "CYP2D6 has 120 annotations" is not.
"""
import re
from collections import Counter

import pytest

from ._harness import report
from atlas.civic import LEGEND as CIVIC_LEGEND

pytestmark = pytest.mark.integration


def test_corpus_coverage_floors(pages):
    """Per-entity page-count floors. Adapts to the two legitimate build sizes:
    the dense gate set (~1k) and the full corpus (~52k). Guards against a whole
    collector silently producing little/nothing — every surviving page can pass
    the structural suite while half the corpus is missing, and only a count
    floor notices."""
    counts = Counter(p.entity for p in pages)
    total = sum(counts.values())
    full = total > 10_000                       # cleanly separates dense (~1k) from full (~52k)
    floors = ({"gene": 27_000, "disease": 17_000, "drug": 4_400} if full
              else {"gene": 400, "disease": 250, "drug": 250})
    bad = [f"{e}: {counts.get(e, 0)} < floor {f}" for e, f in floors.items()
           if counts.get(e, 0) < f]
    assert not bad, (report(bad) +
                     f"\n(total={total}, mode={'full' if full else 'dense'})")


# Must-have genes — a transient biobtree timeout under full-regen load dropped
# TTN and BRCA2 from v1.2.0 (both among the heaviest genes). The count floor
# above can't see a 2-gene loss; this names specific flagships so a silent drop
# fails the full-corpus integration sweep (before the archive is cut).
_FLAGSHIP_GENES = ["TP53", "BRCA1", "BRCA2", "EGFR", "KRAS", "TTN", "CFTR",
                   "CYP2D6", "APOE", "MTHFR"]


def test_flagship_genes_present(pages):
    """On a full-corpus build, a handful of must-have genes must exist. Skips on
    the dense gate set (which doesn't include all of them)."""
    if sum(1 for _ in pages) <= 10_000:
        pytest.skip("dense build — flagship presence is checked on the full corpus")
    have = {p.slug for p in pages if p.entity == "gene"}
    missing = [g for g in _FLAGSHIP_GENES if g not in have]
    assert not missing, f"flagship genes missing from the corpus: {missing}"


# (entity, slug, [(kind, needle)]) — kind is 'contains' or 'absent'. Each needle
# is a SHAPE marker for a feature we shipped, chosen to survive data refreshes.
# Pages absent from the current build (a dense subset) are skipped, not failed.
GOLDEN = [
    ("gene", "CYP2D6", [
        ("contains", "CPIC guideline:"),     # the real PharmGKB signal we keep
        ("absent",   "VIP="),                # the broken is_vip flag we dropped
    ]),
    ("gene", "EGFR", [
        ("contains", CIVIC_LEGEND),          # CIViC A→E legend (EGFR has CIViC)
        ("contains", "/atlas/"),             # cross-entity links render
    ]),
    ("drug", "imatinib", [
        ("contains", "approved indication"), # indications are tiered
        ("absent",   "phase 1–3"),      # investigational floor raised to phase 2–3
        ("contains", "uniprotkb/"),          # bioactivity UniProt links
    ]),
    ("disease", "glioblastoma", [
        ("contains", "ncbi.nlm.nih.gov/snp/"),  # GWAS rsIDs link to dbSNP
        ("contains", "/atlas/gene/"),           # cohort gene symbols link out
    ]),
]


def test_golden_pages(pages):
    """Feature presence/shape on canonical pages — guards the specific things
    this session fixed against silent removal."""
    by_key = {(p.entity, p.slug): p for p in pages}
    found, bad = 0, []
    for entity, slug, checks in GOLDEN:
        p = by_key.get((entity, slug))
        if not p:                            # not in this build (dense subset)
            continue
        found += 1
        for kind, needle in checks:
            present = needle in p.body
            if kind == "contains" and not present:
                bad.append(f"{entity}/{slug}: missing {needle!r}")
            elif kind == "absent" and present:
                bad.append(f"{entity}/{slug}: should not contain {needle!r}")
    if not found:
        pytest.skip("no golden pages in this build")
    assert not bad, report(bad)


def test_civic_legend_present_wherever_rendered(pages):
    """Every page that renders a CIViC table must carry the A→E legend — the
    column is the sort key, so it can't go unexplained. Keyed on the table's
    caption ('predictive associations'), not the {#civic} anchor: the canonical
    section emits that anchor even for the 'No CIViC evidence' placeholder."""
    bad = [f"{p.entity}/{p.slug}" for p in pages
           if "predictive associations" in p.body and CIVIC_LEGEND not in p.body]
    assert not bad, report(bad)


_GWAS_HEAD = "### Top associations by p-value"
_RSID = re.compile(r"\brs\d{2,}\b")


def test_gwas_disease_pages_link_rsids(pages):
    """Among disease pages with a GWAS top-associations table that mentions an
    rsID, the overwhelming majority must link it to dbSNP. Not 100% — a table
    of only non-rs variant ids (chr:pos) legitimately has no dbSNP link — so a
    0.9 floor guards the linking without false alarms on the long tail."""
    rel = [p for p in pages
           if p.entity == "disease" and _GWAS_HEAD in p.body and _RSID.search(p.body)]
    if not rel:
        pytest.skip("no GWAS disease pages with rsIDs in this build")
    linked = [p for p in rel if "ncbi.nlm.nih.gov/snp/" in p.body]
    ratio = len(linked) / len(rel)
    assert ratio >= 0.9, (f"only {len(linked)}/{len(rel)} ({ratio:.0%}) GWAS "
                          "disease pages link rsIDs to dbSNP")


_SHOWING_OF = re.compile(r"showing ([\d,]+) of ([\d,]+)")


def test_disclosure_matches_table_rows(pages):
    """capped_table() emits a SINGLE top caption "showing N of T …:" directly
    above the table — N must equal the table's real row count. Guards the bug
    where the render cap was used as N while the table actually had fewer rows
    (collector cap below the render cap, or table()'s identical-row dedup). A
    small tolerance absorbs residual dedup."""
    bad = []
    for p in pages:
        pending, cur = None, 0          # pending = N claimed by a caption awaiting its table
        for ln in p.body.splitlines():
            s = ln.strip()
            if s.startswith("|"):
                cur += 1
                continue
            if cur:                                  # a table just ended
                if pending is not None:
                    rows = max(0, cur - 2)           # minus header + separator
                    if rows < pending - max(2, pending // 10):
                        bad.append(f"{p.entity}/{p.slug}: 'showing {pending}' "
                                   f"over a {rows}-row table")
                    pending = None
                cur = 0
            if not s:
                continue                             # blank lines keep the link
            m = _SHOWING_OF.search(s)
            if m:
                pending = int(m.group(1).replace(",", ""))
            elif pending is not None:
                pending = None                       # caption not followed by a table
    assert not bad, report(bad)


def test_no_degenerate_truncation_disclosure(pages):
    """A '+N more' / '(showing top K)' disclosure must never be degenerate
    (+0 / top 0 / negative) — that would mean the cap math produced nonsense."""
    bad = [f"{p.entity}/{p.slug}" for p in pages
           if "+0 more" in p.body or "showing top 0)" in p.body or "(+0 more" in p.body]
    assert not bad, report(bad)
