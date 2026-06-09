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


def with_heading_id(md, anchor):
    """Add an explicit `{#anchor}` to the FIRST heading of a rendered section,
    so the H3 anchor is a stable backend-owned id (not Hugo's autoHeadingID,
    which is derived from the prose — `…generif-showing-40` breaks when the count
    changes). No-op if the section is empty or its first heading already carries
    an explicit id."""
    if not md:
        return md
    lines = md.split("\n")
    for i, ln in enumerate(lines):
        if re.match(r"^#{2,6} ", ln):
            if "{#" not in ln:
                lines[i] = ln.rstrip() + f" {{#{anchor}}}"
            break
    return "\n".join(lines)


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
    apply to gene symbols — 'TP53'.isupper() is True but must stay upper. Also
    leave drug development CODES alone — anything containing a digit (N6022,
    K-877, F19 131I, ABT199) is an identifier, not a shouting word, and must
    keep its exact case."""
    if not (s and isinstance(s, str) and s.isupper()):
        return s
    if any(c.isdigit() for c in s):           # code, not a word — preserve
        return s
    return s.title()


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


def _cell(c):
    """A markdown table cell: unescape HTML entities (UniProt names carry
    &alpha; etc.) and escape literal pipes — UniProt cleavage notation embeds
    '|' (e.g. "1838-Glu-|-Ala-1839"), which would otherwise shift the column."""
    if c is None:
        return ""
    s = html.unescape(str(c))
    if s.strip().lower() == "nan":          # biobtree returns the literal string
        return ""                           # "nan" for a missing value (e.g. an
                                            # ortholog row with no gene symbol)
    return s.replace("|", "\\|")


def more_line(total, shown, by=None):
    """The standard truncation-disclosure line for a capped table: when `total`
    exceeds `shown`, return a `*+N more (showing top X…)*` italic note, else "".
    Every capped table uses this so the corpus never prints a full count above a
    table that silently stops short. `by` optionally names the sort key
    ("by evidence level", "by phase") for the parenthetical."""
    try:
        total, shown = int(total), int(shown)
    except (TypeError, ValueError):
        return ""
    # shown <= 0 means the table rendered no rows (e.g. a collector that has a
    # count but an empty sample) — there's nothing to be "more than", so don't
    # emit a degenerate "(showing top 0)" line; the count already sits in the caption.
    if shown <= 0 or total <= shown:
        return ""
    tail = f" {by}" if by else ""
    return f"\n*+{total - shown} more (showing top {shown}{tail}).*"


def table(headers, rows):
    """GitHub-flavored markdown table. Empty cells blank; literal pipes in
    values escaped (no column shift); identical data rows collapsed (source
    sometimes repeats a row verbatim — Reactome pathway, Orphanet prevalence).
    All-blank rows are dropped, and a table with no surviving data rows renders
    as "" (not a dangling header-only table) — sub-blocks that turn out empty
    elide cleanly instead of printing `| Field | … |` with nothing under it."""
    out = ["| " + " | ".join(headers) + " |",
           "| " + " | ".join("---" for _ in headers) + " |"]
    seen = set()
    for r in rows:
        cells = tuple(_cell(c) for c in r)
        if not any(cells) or cells in seen:
            continue
        seen.add(cells)
        out.append("| " + " | ".join(cells) + " |")
    if len(out) == 2:                       # header only, no data rows
        return ""
    return "\n".join(out)
