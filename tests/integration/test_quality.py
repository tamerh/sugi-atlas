"""Data-quality guards over the built corpus — each maps to a fix we shipped,
so a regression reappears as a failing gate here."""
import re

import pytest

from ._harness import report

pytestmark = pytest.mark.integration

_NAN_CELL = re.compile(r"\|\s*nan\s*\|", re.I)             # "| nan |" / "| NaN |"
# The UNAMBIGUOUS float32-repr signature — a long run of 9s or 0s in the decimal
# tail (0.20799999999996, 6.170000076293945). Legit precise/small values
# (0.00133 nM, MW 493.6038) never have this. A blunt "long decimal" regex can't
# tell float32 noise from a genuine multi-sig-fig measurement, so the targeted
# fnum rounding (odds_ratio/maf/resolution, unit-tested) is the real guard; this
# is the corpus net for the egregious storage-noise case.
_FLOAT_ARTIFACT = re.compile(r"\.\d*(?:9{6,}|0{5,}[1-9])")
_DUP_INCLUDING = re.compile(r"including (\w[\w '-]+?) and \1\b", re.I)
_NONE_CELL = re.compile(r"\|\s*None\s*\|")
_BAD_ATTR = re.compile(r"\{#")                             # {#id} only belongs on headings


def test_no_leaked_nan_in_tables(pages):
    bad = [f"{p.entity}/{p.slug}" for p in pages if _NAN_CELL.search(p.body)]
    assert not bad, report(bad)


def test_no_float32_artifacts(pages):
    bad = []
    for p in pages:
        m = _FLOAT_ARTIFACT.search(p.body)
        if m:
            bad.append(f"{p.entity}/{p.slug}: …{m.group(0)}…")
    assert not bad, report(bad)


def test_no_x_and_x_dup_in_prose(pages):
    bad = []
    for p in pages:
        m = _DUP_INCLUDING.search(p.body)
        if m:
            bad.append(f"{p.entity}/{p.slug}: '{m.group(0)}'")
    assert not bad, report(bad)


def test_no_none_cells(pages):
    bad = [f"{p.entity}/{p.slug}" for p in pages if _NONE_CELL.search(p.body)]
    assert not bad, report(bad)


def test_anchor_attrs_only_on_headings(pages):
    """`{#id}` must render only as a heading attribute — a leak into prose means
    the goldmark-attribute syntax landed somewhere it won't be parsed."""
    bad = []
    for p in pages:
        for line in p.body.splitlines():
            if "{#" in line and not line.lstrip().startswith("#"):
                bad.append(f"{p.entity}/{p.slug}: {line.strip()[:60]}")
                break
    assert not bad, report(bad)


def test_pages_non_trivial(pages):
    """No empty/truncated pages."""
    bad = [f"{p.entity}/{p.slug}: {len(p.raw)}b" for p in pages if len(p.raw) < 2000]
    assert not bad, report(bad)
