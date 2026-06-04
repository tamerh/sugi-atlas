"""Frontmatter schema (P2/P3) over the built corpus."""
import pytest

from collections import defaultdict

from ._harness import ID_RE, H2_IDS, report

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
        # All-caps WORDS shout (GLEEVEC); all-caps codes with a digit (N6022,
        # K-877) are identifiers and stay as-is.
        if (p.entity != "gene" and len(t) > 3 and t == t.upper()
                and not any(c.isdigit() for c in t)):
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
    """Soft gate: the vast majority of pages should carry a TL;DR. The full
    corpus has a long tail of legitimately-thin entities (ncRNA genes, sparse
    rare diseases) with no key facts to summarize, so the floor is 0.85 — the
    dense set sat near 1.0; this only catches a corpus-wide tldr regression."""
    with_tldr = sum(1 for p in pages if p.fm.get("tldr"))
    frac = with_tldr / len(pages)
    assert frac >= 0.85, f"only {frac:.0%} of pages have a tldr ({with_tldr}/{len(pages)})"


def test_evidence_score_present_and_typed(pages):
    """Every page carries a 0-100 integer `evidence_score` and a
    `evidence_components` map of integer counts whose keys are STABLE across
    all pages of a type (the contract the web search ranking relies on)."""
    bad = []
    type_keys = {}
    for p in pages:
        prom = p.fm.get("evidence_score")
        if not isinstance(prom, int) or not (0 <= prom <= 100):
            bad.append(f"{p.entity}/{p.slug}: evidence_score={prom!r}")
        comps = p.fm.get("evidence_components")
        if not isinstance(comps, dict) or not comps:
            bad.append(f"{p.entity}/{p.slug}: evidence_components missing/empty")
            continue
        if any(not isinstance(v, int) for v in comps.values()):
            bad.append(f"{p.entity}/{p.slug}: non-integer component value")
        keys = frozenset(comps)
        type_keys.setdefault(p.entity, keys)
        if keys != type_keys[p.entity]:
            bad.append(f"{p.entity}/{p.slug}: component keys differ for type")
    assert not bad, report(bad)


def test_section_defaults_keys_are_valid_anchors(pages):
    """section_defaults hints reference real canonical anchor ids."""
    bad = []
    for p in pages:
        for k in (p.fm.get("section_defaults") or {}):
            if k not in H2_IDS[p.entity]:
                bad.append(f"{p.entity}/{p.slug}: section_defaults key '{k}'")
    assert not bad, report(bad)


def test_identifier_unique_per_entity(pages):
    """No two pages of a type share an identifier (slug/id collision)."""
    seen = defaultdict(list)
    for p in pages:
        seen[(p.entity, str(p.fm.get("identifier")))].append(p.slug)
    bad = [f"{et} {ident}: {slugs}" for (et, ident), slugs in seen.items()
           if ident not in ("None", "") and len(slugs) > 1]
    assert not bad, report(bad)
