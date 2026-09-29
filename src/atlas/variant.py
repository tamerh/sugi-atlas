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

# The variant-slug contract is SHARED with Sugi Variant via the `sugislug` package
# (github.com/tamerh/sugi-slug, editable-installed like sugibiobtree) — so a slug
# Atlas builds is byte-identical to Sugi Variant's own and a deep link always
# resolves. Sugi Variant owns the contract + versions it; Atlas consumes.
from sugislug import variant_slugs, parse_hgvs, name_gene

_BASE = "https://sugi.bio/variant"

# ClinVar germline classifications Sugi Variant builds a per-variant page for (its
# pathogenic-gated corpus). A variant outside this set has NO page — gate on it
# before emitting a per-variant deep link, or the link 404s.
_VARIANT_CORPUS_CLASSES = {
    "pathogenic", "likely pathogenic", "pathogenic/likely pathogenic",
    "conflicting classifications of pathogenicity",
}


def in_variant_corpus(classification):
    """True if a ClinVar classification is in Sugi Variant's pathogenic-gated corpus."""
    return (classification or "").strip().lower() in _VARIANT_CORPUS_CLASSES


def gene_variants_url(symbol):
    """Sugi Variant per-gene variant index URL for a gene symbol, or None for an
    empty symbol. Gene symbols are canonical upper-case (TP53); Sugi Variant keys
    its index on the same symbol."""
    s = (symbol or "").strip()
    return f"{_BASE}/gene/{quote(s)}" if s else None


def variant_page_url(clinvar_name, gene=None):
    """Sugi Variant per-variant page URL, derived from a ClinVar `name`
    (e.g. 'NM_002485.5(NBN):c.1889C>A (p.Ser630Ter)') via the shared sugislug
    contract → byte-identical to Sugi Variant's own slug, so the link resolves.
    Uses the transcript gene from the name (name_gene) for correctness, falling
    back to `gene`. None when the name yields no slug. The CALLER gates on
    classification (in_variant_corpus) — Sugi Variant only builds pathogenic-tier
    pages, so link only those."""
    name = clinvar_name or ""
    g = name_gene(name) or gene
    if not g:
        return None
    c_form, p_form = parse_hgvs(name)
    slug, _ = variant_slugs(g, c_form, p_form)
    return f"{_BASE}/{slug}" if slug else None
