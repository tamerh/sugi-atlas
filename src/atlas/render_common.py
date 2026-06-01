"""Shared markdown rendering helpers — used by both atlas.gene.render and
atlas.disease.render so we don't duplicate the table primitive.

Keep this module tiny: zero imports beyond stdlib. Anything section-aware
belongs in the entity-specific renderer."""
import html


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
