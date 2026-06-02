"""Frontmatter schema (P2/P3) over the built corpus."""
import pytest

from ._harness import ID_RE, report

pytestmark = pytest.mark.integration

_REQUIRED = ("title", "identifier", "entity_type", "generated_at",
             "atlas_version", "generated_by")


def test_frontmatter_parses(pages):
    bad = [f"{p.entity}/{p.slug}" for p in pages
           if not p.fm or p.fm.get("__yaml_error__")]
    assert not bad, report(bad)


def test_required_fields_present(pages):
    bad = [f"{p.entity}/{p.slug}: missing {k}"
           for p in pages for k in _REQUIRED if not p.fm.get(k)]
    assert not bad, report(bad)


def test_identifier_is_typed(pages):
    """gene → HGNC symbol, disease → MONDO:n, drug → CHEMBLn (the trustworthy
    cross-type key templates rely on)."""
    bad = []
    for p in pages:
        ident = str(p.fm.get("identifier") or "")
        if not ID_RE[p.entity].match(ident):
            bad.append(f"{p.entity}/{p.slug}: identifier={ident!r}")
    assert not bad, report(bad)


def test_no_hugo_reserved_aliases_key(pages):
    """`aliases:` would make Hugo emit 301 redirects — search aliases live in
    `alt_names:`."""
    bad = [f"{p.entity}/{p.slug}" for p in pages if "aliases" in (p.fm or {})]
    assert not bad, report(bad)


def test_title_not_shouting(pages):
    """Disease/drug titles are de-SHOUTed (#12). Genes are exempt — HGNC symbols
    are upper by convention."""
    bad = []
    for p in pages:
        t = str(p.fm.get("title") or "")
        if p.entity != "gene" and len(t) > 3 and t == t.upper():
            bad.append(f"{p.entity}/{p.slug}: {t!r}")
    assert not bad, report(bad)


def test_list_and_map_fields_well_typed(pages):
    """When present: alt_names/tldr are lists, section_defaults is a map with
    summary open."""
    bad = []
    for p in pages:
        for k in ("alt_names", "tldr"):
            if k in p.fm and not isinstance(p.fm[k], list):
                bad.append(f"{p.entity}/{p.slug}: {k} not a list")
        sd = p.fm.get("section_defaults")
        if sd is not None and (not isinstance(sd, dict) or sd.get("summary") != "open"):
            bad.append(f"{p.entity}/{p.slug}: section_defaults={sd}")
    assert not bad, report(bad)


def test_tldr_coverage(pages):
    """Soft gate: the vast majority of pages should carry a TL;DR."""
    with_tldr = sum(1 for p in pages if p.fm.get("tldr"))
    frac = with_tldr / len(pages)
    assert frac >= 0.9, f"only {frac:.0%} of pages have a tldr ({with_tldr}/{len(pages)})"
