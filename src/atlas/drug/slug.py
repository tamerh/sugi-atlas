"""Drug-name → URL/file slug. Deterministic, lowercase, ASCII-only.

Used for:
- body_gate snapshot filename (snapshots/drug/<slug>.json)
- dist output directory (sugi-atlas-dist/atlas/drug/<slug>/)
- Hugo URL (sugi.bio/atlas/drug/<slug>/)

Drugs slug from their canonical ChEMBL name (e.g. "IMATINIB" -> "imatinib").
Mirrors disease/slug.py; kept separate per the parallel-not-generalized
entity philosophy.
"""
import re
import unicodedata


def slugify(name: str) -> str:
    """Reproducible slug. Stable across runs — only renames if `name` changes."""
    s = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode("ascii")
    s = s.lower()
    s = re.sub(r"['’]", "", s)
    s = re.sub(r"[^a-z0-9]+", "-", s)
    s = re.sub(r"-+", "-", s).strip("-")
    return s
