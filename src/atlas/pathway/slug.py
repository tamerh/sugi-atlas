"""Reactome pathway → URL/file slug. Deterministic, lowercase, ASCII-only.

Pathways slug from their Reactome name (e.g. "Apoptosis" -> "apoptosis"). The
stable-rename property comes from the name; the Reactome id (R-HSA-…) is the
durable identity and is carried in the page body for disambiguation. Mirrors
drug/slug.py per the parallel-not-generalized entity philosophy.
"""
import re
import unicodedata


def slugify(name: str) -> str:
    s = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode("ascii")
    s = s.lower()
    s = re.sub(r"['’]", "", s)
    s = re.sub(r"[^a-z0-9]+", "-", s)
    s = re.sub(r"-+", "-", s).strip("-")
    return s
