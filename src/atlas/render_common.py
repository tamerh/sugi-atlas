"""Shared markdown rendering helpers — used by both atlas.gene.render and
atlas.disease.render so we don't duplicate the table primitive.

Keep this module tiny: zero imports beyond stdlib. Anything section-aware
belongs in the entity-specific renderer."""
import html


def fnum(v, nd=2):
    """Round float-ish display values to `nd` decimals so float32 artifacts
    (e.g. 6.170000076293945, 0.19679999999999997) don't leak into pages.
    Ints, None, and non-numeric strings pass through unchanged; an integral
    result renders without a trailing '.0'."""
    if isinstance(v, bool) or v is None:
        return v
    f = None
    if isinstance(v, (int, float)):
        f = float(v)
    elif isinstance(v, str):
        try:
            f = float(v)
        except ValueError:
            return v
    if f is None:
        return v
    r = round(f, nd)
    return int(r) if r == int(r) else r


def phase_label(p):
    """Normalize a clinical-trial phase for display. biobtree emits 'NaN' for
    trials with no interventional phase (observational / natural-history), which
    naively uppercases to a confusing 'NAN' — map those to 'Not specified'."""
    p = (p or "").strip().upper()
    return "Not specified" if p in ("", "NAN", "NA") else p


def table(headers, rows):
    """GitHub-flavored markdown table. Empty cells render as blank; HTML
    entities in cell values are unescaped (UniProt names often contain &alpha;
    etc. from the upstream record)."""
    out = ["| " + " | ".join(headers) + " |",
           "| " + " | ".join("---" for _ in headers) + " |"]
    for r in rows:
        out.append("| " + " | ".join("" if c is None else html.unescape(str(c)) for c in r) + " |")
    return "\n".join(out)
