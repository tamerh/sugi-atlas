"""Shared markdown rendering helpers — used by both atlas.gene.render and
atlas.disease.render so we don't duplicate the table primitive.

Keep this module tiny: zero imports beyond stdlib. Anything section-aware
belongs in the entity-specific renderer."""
import html
import re

_HEADING = re.compile(r'^(#{2,5}) ', re.M)


def demote(md):
    """Bump every ATX heading one level deeper (## → ###) so a sub-section nests
    under its canonical H2. Non-heading lines pass through unchanged."""
    return _HEADING.sub(lambda m: "#" + m.group(0), md) if md else md


def emit_canonical(spec, anchors=None):
    """Emit the FROZEN canonical H2 sequence (docs/PAGE_CONTRACT.md) from a list
    of (label, id, body, placeholder): `## label {#id}` in the given order, body
    or an informative `*placeholder*` when empty — every section always emitted
    so the TOC is identical across every page of a type. `anchors` optionally
    maps an id → raw HTML prepended before its heading (e.g. the JSON-LD `@id`
    <a> for #protein). Sub-section headings are expected pre-demoted to H3."""
    anchors = anchors or {}
    out = []
    for label, anchor, body, placeholder in spec:
        body = (body or "").strip()
        content = body or (f"*{placeholder}*" if placeholder else "")
        if not content:
            continue
        out.append(f"{anchors.get(anchor, '')}## {label} {{#{anchor}}}\n\n{content}")
    return "\n\n".join(out)

# A raw ontology accession (MONDO:0004992, EFO:0010282, MP:0001914) leaking
# where a human-readable label belongs (audit #11). MP is a mouse-phenotype id,
# not even a human disease — such rows are unmapped noise, not data.
_ONTOLOGY_ID = re.compile(
    r'^(mondo|efo|mesh|hp|hpo|doid|umls|orphanet|orpha|ncit|snomedct|snomed'
    r'|meddra|mp|go|chebi|omim|gard|medgen|icd\d*)[:_]', re.I)


def is_ontology_id(s) -> bool:
    """True when `s` looks like a raw ontology accession rather than a label."""
    return bool(s) and bool(_ONTOLOGY_ID.match(str(s).strip()))


def display_name(s):
    """Title-case an all-caps (SHOUTING) label for display (audit #12: ChEMBL
    names like 'IMATINIB'/'WATER'); leave mixed-case strings untouched. NEVER
    apply to gene symbols — 'TP53'.isupper() is True but must stay upper."""
    return s.title() if (s and isinstance(s, str) and s.isupper()) else s


# GenCC classification strength, strongest first (audit #13 dedup ranking).
_GENCC_RANK = {"definitive": 6, "strong": 5, "moderate": 4, "supportive": 3,
               "limited": 2, "disputed evidence": 1, "refuted": 0,
               "animal model only": 0, "no known disease relationship": 0}


def gencc_rank(c):
    """Numeric strength of a GenCC classification label (higher = stronger)."""
    return _GENCC_RANK.get((c or "").strip().lower(), 0)


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
