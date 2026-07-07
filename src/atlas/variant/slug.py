"""Deterministic slugs for variant pages, from GENE + HGVS.

biobtree can't resolve HGVS free-text (`search("p.Pro334Arg")` → 0), so variant
pages are enumerated gene-first and the URL must be derivable from `GENE` + the
normalized HGVS — in BOTH the protein (p.) and coding (c.) forms, since the
search demand comes as both `acta1 "pro309ala"` and `"c.925c>g" acta1`.
Canonical = the p. form when present (the more-searched shape); the c. form is a
same-page alias.
"""
import re

# ClinVar `name` looks like: "NM_001100.4(ACTA1):c.983_985del (p.Lys328del)".
_C_RE = re.compile(r"(c\.[^\s)]+)")
_P_RE = re.compile(r"\(p\.([^)]+)\)")


def parse_hgvs(name):
    """(c_form, p_form) from a ClinVar `name`; either may be None."""
    c = _C_RE.search(name or "")
    p = _P_RE.search(name or "")
    return (c.group(1) if c else None), (f"p.{p.group(1)}" if p else None)


def _norm(s):
    """Slug-normalize an HGVS/gene token: lowercase, non-alnum → '-', collapse."""
    s = (s or "").lower()
    s = re.sub(r"[^a-z0-9]+", "-", s)
    return s.strip("-")


def variant_slugs(gene, c_form, p_form):
    """(canonical_slug, [all_slugs]) — canonical is the p. form when present,
    else c.; every distinct form is emitted as an alias so either query shape
    lands the same page. None,[] if neither HGVS parses."""
    g = _norm(gene)
    if not g:
        return None, []
    out = []
    for form in (p_form, c_form):        # p. first → canonical
        if form:
            slug = f"{g}-{_norm(form)}"
            if slug not in out:
                out.append(slug)
    return (out[0] if out else None), out
