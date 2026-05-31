"""UniProt CC (`comments`) narrative helpers.

biobtree's 2026-05-31 refresh exposed the CC block on uniprot entries
(BIOBTREE_ISSUES.md #9 RESOLVED). The text contains UniProt evidence
codes (`(PubMed:NNNN, PubMed:NNNN, ...)`, `{ECO:0000269|PubMed:NNNN}`,
`(By similarity)`) that are noise for page rendering. This module:
  - fetches the comments + isoforms block,
  - strips evidence codes for human-facing display,
  - extracts the lead sentence for declarative-lead + JSON-LD description use.

Shared between gene + disease pipelines.
"""
import re
from atlas.biobtree import entry

# Evidence-code patterns to strip from CC text:
#   (PubMed:11111, PubMed:22222, ...)
#   (PubMed:11111)
#   {ECO:0000269|PubMed:NNNN}
#   (By similarity)
_PUBMED_GROUP = re.compile(r"\s*\(\s*(?:PubMed:\d+\s*,?\s*)+\)")
_ECO          = re.compile(r"\s*\{ECO:[^}]+\}")
_BY_SIMILARITY = re.compile(r"\s*\((?:By similarity|Probable)\)")

# Note that some CCs use bracketed reference codes like '[MIM:151623]'.
# These ARE useful (link to OMIM) so we keep them.

# Subcellular_location often has '. Note=...' suffixes — those are useful
# context; keep. Tissue_specificity may end in periods + literature trailers
# (PubMed:NNNN). We strip evidence trailers only.

def strip_evidence_codes(text: str) -> str:
    """Remove UniProt evidence-code clutter for human display."""
    if not text:
        return ""
    t = _PUBMED_GROUP.sub("", text)
    t = _ECO.sub("", t)
    t = _BY_SIMILARITY.sub("", t)
    # Collapse double spaces, fix space-before-punct introduced by stripping.
    t = re.sub(r"\s+", " ", t)
    t = re.sub(r"\s+([.,;:])", r"\1", t)
    return t.strip()


# Split CC text into sentences for first-sentence extraction. UniProt's
# function blocks are sentence-rich; we want the LEAD only for the
# declarative-lead + JSON-LD description.
_SENT_RX = re.compile(r"(?<=[.!?])\s+(?=[A-Z(])")

def first_sentence(text: str, max_chars: int = 280) -> str:
    """Return the first sentence of `text` (post-evidence-strip), truncated
    if longer than max_chars. Falls back to the whole text if no sentence
    boundary is found within the cap."""
    s = strip_evidence_codes(text)
    if not s:
        return ""
    parts = _SENT_RX.split(s, maxsplit=1)
    first = parts[0].strip()
    if len(first) > max_chars:
        return first[: max_chars - 1].rstrip() + "…"
    return first


def fetch_cc(uniprot_acc: str) -> dict:
    """Returns the raw UniProt comments + isoforms dict from biobtree's entry.

    Shape:
      {
        "comments": {"function": "...", "subunit": "...", ...},
        "isoforms": [{"id": "P04637-1", "names": [...], "is_canonical": True}, ...],
        "alternative_names": [...],
        "name": "Cellular tumor antigen p53",  # primary name
      }
    Returns an empty dict when the accession is unreviewed (biobtree returns
    empty Attributes) or on fetch error."""
    try:
        en = entry(uniprot_acc, "uniprot")
    except Exception:
        return {}
    attrs = ((en.get("Attributes") or {}).get("Uniprot") or {})
    return {
        "comments": attrs.get("comments") or {},
        "isoforms": attrs.get("isoforms") or [],
        "alternative_names": attrs.get("alternative_names") or [],
        "name": (attrs.get("names") or [None])[0],
    }
