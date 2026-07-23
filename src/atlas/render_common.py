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
    if " " not in s and any(c.isdigit() for c in s):   # single-token code — preserve
        return s                                        # (multi-word names w/ digits ARE names)
    return s.title()


# GenCC classification strength, strongest first (audit #13 dedup ranking).
_GENCC_RANK = {"definitive": 6, "strong": 5, "moderate": 4, "supportive": 3,
               "limited": 2, "disputed evidence": 1, "refuted": 0,
               "animal model only": 0, "no known disease relationship": 0}


def gencc_rank(c):
    """Numeric strength of a GenCC classification label (higher = stronger)."""
    return _GENCC_RANK.get((c or "").strip().lower(), 0)


def pval(s):
    """Tidy a GWAS-style p-value string for display: "8.000000e-11" → "8e-11",
    "1.500000e-08" → "1.5e-8". Passes non-numeric / empty through unchanged."""
    try:
        v = float(s)
    except (TypeError, ValueError):
        return s or ""
    if v <= 0:
        return s or ""
    mant, _, e = f"{v:.1e}".partition("e")
    return f"{mant.rstrip('0').rstrip('.')}e{int(e)}"


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
    """Normalize a clinical-trial phase for display: the raw CT.gov enum
    ('PHASE2', 'EARLY_PHASE1', 'PHASE1/PHASE2') → human form ('Phase 2', 'Early
    Phase 1', 'Phase 1/2'). Idempotent (an already-formatted 'Phase 3' passes
    through unchanged). biobtree emits 'NaN' for trials with no interventional
    phase (observational / natural-history) → 'Not specified'."""
    u = (p or "").strip().upper()
    if u in ("", "NAN", "NA", "N/A", "NONE"):
        return "Not specified"
    if u.replace(" ", "").replace("_", "").startswith("EARLYPHASE"):
        n = re.findall(r"\d", u)
        return "Early Phase " + n[0] if n else "Early Phase"
    nums = re.findall(r"PHASE\s*([0-9])", u)   # also matches already-formatted 'PHASE 3'
    if nums:
        return "Phase " + "/".join(nums)
    return (p or "").strip()                    # unknown label — leave as-is


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
    # Escape '*' too: chemical / patent-compound names embed it (stereochem
    # descriptors like '(R*)'), which markdown reads as emphasis and can leave a
    # literal, unbalanced '**' in the rendered HTML (test_no_unbalanced_bold — a
    # BindingDB ligand on CX3CR1 surfaced it in the full corpus). Table cells
    # carry data + links only, never intentional bold/italic (verified corpus-
    # wide: zero legitimate '**' cells), and no cell URL contains '*', so escaping
    # every '*' in a cell is safe and closes the whole class at the one choke point.
    return s.replace("|", "\\|").replace("*", "\\*")


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


def capped_table(headers, rows, cap, total=None, noun=""):
    """The one-stop capped table — the single source of truth for "this table is
    truncated" across the corpus. Slices `rows` to `cap`, renders, and prefixes
    ONE top caption: "showing N of T {noun}:" when truncated, else "T {noun}:".

    N is the ACTUAL number of rows rendered (computed AFTER table()'s identical-
    row dedup), so the count can never over-claim — and it lives in exactly one
    place (above the table), not split between a top description and a trailing
    "+N more". T is `total` (the full available count) or len(rows). Returns the
    markdown block, or "" when the table is empty so the sub-block elides.

    Callers pass the FULL row list + the cap; the helper owns the slice, the
    dedup-aware count, and the caption format."""
    t = len(rows) if total is None else (total or 0)
    md = table(headers, rows[:cap] if cap is not None else rows)
    if not md:
        return ""
    shown = max(0, md.count("\n") - 1)               # lines = N data + header + sep
    noun = (noun or "").strip()
    lead = f"showing {shown:,} of {t:,}" if t > shown else f"{t:,}"
    head = f"{lead} {noun}".rstrip()
    return f"{head}:\n\n{md}"


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
