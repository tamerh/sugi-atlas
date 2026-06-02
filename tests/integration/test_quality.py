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


def test_tables_have_consistent_columns(pages):
    """Every table row has the header's column count — guards the unescaped-pipe
    column-shift corruption class (the MeSH-row bug)."""
    bad = []
    for p in pages:
        for ti, tbl in enumerate(p.tables()):
            if len(tbl) < 2:
                continue
            ncol = len(tbl[0])
            for ri, row in enumerate(tbl[1:], 1):
                if len(row) != ncol:
                    bad.append(f"{p.entity}/{p.slug} table#{ti} row{ri}: "
                               f"{len(row)}≠{ncol} cols {row}")
                    break
    assert not bad, report(bad)


_ONTOLOGY_LABEL = re.compile(
    r"\[((?:MONDO|EFO|MESH|MP|HP|HPO|DOID|ORPHANET|ORPHA|NCIT|GO|CHEBI|OMIM)[:_]\d+)\]\(",
    re.I)


def test_no_duplicate_table_rows(pages):
    """No identical data row repeated in a table (the GenCC ×19 dedup class)."""
    bad = []
    for p in pages:
        for ti, tbl in enumerate(p.tables()):
            if len(tbl) < 3:
                continue
            seen = set()
            for row in tbl[2:]:                       # skip header + separator
                key = tuple(c.strip() for c in row)
                if any(key) and key in seen:
                    bad.append(f"{p.entity}/{p.slug} table#{ti}: dup row {list(key)[:3]}")
                    break
                seen.add(key)
    assert not bad, report(bad)


_ENTITY = re.compile(r"&(?:[a-zA-Z]{2,}|#\d+);")


def test_no_html_entity_leaks(pages):
    """HTML entities (&alpha;, &amp;) are unescaped at render — a leak means a
    table cell or label bypassed the unescape."""
    bad = []
    for p in pages:
        text = re.sub(r"<script.*?</script>", "", p.body, flags=re.S)  # JSON-LD ok
        m = _ENTITY.search(text)
        if m:
            i = m.start()
            bad.append(f"{p.entity}/{p.slug}: …{text[max(0, i-15):i+12]}…")
    assert not bad, report(bad)


def test_no_ontology_id_as_link_label(pages):
    """A linked name must be a label, never a raw ontology accession
    (the #11 'MP:0001914 as a disease name' class)."""
    bad = []
    for p in pages:
        m = _ONTOLOGY_LABEL.search(p.body)
        if m:
            bad.append(f"{p.entity}/{p.slug}: [{m.group(1)}](…)")
    assert not bad, report(bad)
