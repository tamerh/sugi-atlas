"""Frozen page contract (docs/PAGE_CONTRACT.md) over the built corpus."""
import re

import pytest

from ._harness import H2_IDS, H3_IDS, H4_IDS, report

pytestmark = pytest.mark.integration

_HEADING = re.compile(r"^#{2,6} .+?\{#([a-z0-9-]+)\}\s*$")


def test_canonical_h2_set_and_order(pages):
    """Every page emits exactly its entity's canonical H2 ids, in frozen order."""
    bad = []
    for p in pages:
        ids = [i for _label, i in p.h2]
        if ids != H2_IDS[p.entity]:
            bad.append(f"{p.entity}/{p.slug}: {ids}")
    assert not bad, report(bad)


def test_every_h2_has_explicit_kebab_id(pages):
    """No auto-derived heading ids — every '## ' carries an explicit {#id}."""
    bad = [f"{p.entity}/{p.slug}: '## {label}'"
           for p in pages for label, i in p.h2 if not i]
    assert not bad, report(bad)


def test_summary_first_related_last(pages):
    bad = []
    for p in pages:
        ids = [i for _l, i in p.h2 if i]
        if not ids or ids[0] != "summary" or ids[-1] != "related":
            bad.append(f"{p.entity}/{p.slug}: {ids[:1]}…{ids[-1:]}")
    assert not bad, report(bad)


def test_section_h3_have_explicit_ids(pages):
    """Section H3 and H4 headings carry backend-owned {#id}s — never Hugo's
    prose-derived autoHeadingID (…generif-showing-40 breaks on prose change).
    H4s are now table-block deep-link targets (gene §drug-data), so they need
    explicit ids too."""
    bad = []
    for p in pages:
        for line in p.body.splitlines():
            if (line.startswith("### ") or line.startswith("#### ")) and "{#" not in line:
                bad.append(f"{p.entity}/{p.slug}: '{line.strip()[:55]}'")
    assert not bad, report(bad)


def test_section_h3_ids_match_contract(pages):
    """Section H3 ids stay within the frozen contract set — a new/renamed id
    (anchor-API drift) fails here, not silently in the wild."""
    bad = []
    for p in pages:
        ids = {i for _l, i in p.h3 if i}
        extra = ids - H3_IDS[p.entity]
        if extra:
            bad.append(f"{p.entity}/{p.slug}: not in contract → {sorted(extra)}")
    assert not bad, report(bad)


def test_section_h4_ids_match_contract(pages):
    """H4 table-block ids (e.g. gene §drug-data: #bindingdb, #civic) stay within
    the frozen set — a new/renamed anchor fails here, not silently in the wild."""
    bad = []
    for p in pages:
        ids = {i for _l, i in p.h4 if i}
        extra = ids - H4_IDS[p.entity]
        if extra:
            bad.append(f"{p.entity}/{p.slug}: not in contract → {sorted(extra)}")
    assert not bad, report(bad)


def test_no_duplicate_anchor_ids(pages):
    """Anchor ids are a stable API — duplicates within a page break deep links."""
    bad = []
    for p in pages:
        ids = [m.group(1) for line in p.body.splitlines()
               if (m := _HEADING.match(line))]
        dups = sorted({x for x in ids if ids.count(x) > 1})
        if dups:
            bad.append(f"{p.entity}/{p.slug}: {dups}")
    assert not bad, report(bad)
