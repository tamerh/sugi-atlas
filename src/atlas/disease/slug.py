"""Disease-name → URL/file slug. Deterministic, lowercase, ASCII-only.

Used for:
- body_gate snapshot filename (snapshots/disease/<slug>.json)
- dist output directory (sugi-atlas-dist/atlas/disease/<slug>/)
- Hugo URL (sugi.bio/atlas/disease/<slug>/)

Examples:
  "age-related macular degeneration" -> "age-related-macular-degeneration"
  "Alzheimer's disease"              -> "alzheimers-disease"
  "Type 2 Diabetes Mellitus"         -> "type-2-diabetes-mellitus"
"""
import re
import unicodedata


def slugify(name: str) -> str:
    """Reproducible slug. Stable across runs — only renames if `name` changes."""
    # Strip accents/diacritics
    s = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode("ascii")
    s = s.lower()
    s = re.sub(r"['’]", "", s)        # apostrophes drop (alzheimer's -> alzheimers)
    s = re.sub(r"[^a-z0-9]+", "-", s)      # everything non-alnum becomes a separator
    s = re.sub(r"-+", "-", s).strip("-")   # collapse + trim
    return s
