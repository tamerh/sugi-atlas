"""Sugi Variant cross-links — Atlas → the sibling product at sugi.bio/variant.

Sugi Variant is a deterministic per-variant reference over the pathogenic-gated
slice of ClinVar (Pathogenic / Likely pathogenic / Conflicting — ~464k variants,
6,537 genes; incl. ncRNA and mitochondrial). Each covered gene has a per-gene
variant index at /variant/gene/{SYMBOL}; a gene outside that gate has no page.

Coverage is NOT checked here (no cross-service manifest): the caller gates on the
gene's OWN Pathogenic/Likely-pathogenic ClinVar count — which is exactly Sugi
Variant's corpus gate — so a gene that clears the gate always resolves. The
reverse direction (Variant → Atlas gene/disease) already exists in
sugivariant/links.py, so this closes the loop.
"""
from urllib.parse import quote

_BASE = "https://sugi.bio/variant"


def gene_variants_url(symbol):
    """Sugi Variant per-gene variant index URL for a gene symbol, or None for an
    empty symbol. Gene symbols are canonical upper-case (TP53); Sugi Variant keys
    its index on the same symbol."""
    s = (symbol or "").strip()
    return f"{_BASE}/gene/{quote(s)}" if s else None
