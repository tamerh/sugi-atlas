"""Frozen page contract (docs/PAGE_CONTRACT.md) over the built corpus."""
import re

import pytest

from ._harness import H2_IDS, report

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
