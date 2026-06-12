"""§8 — interactions: PPIs (STRING/IntAct/BioGRID via interaction records,
score-bearing), SIGNOR signaling, ESM2/Diamond similarity.

Hot-path note: STRING/IntAct return rows already sorted by score desc
(verified live). BioGRID has no score column (no ranking claim made on
that source). Render only shows top 30 per source — so paginating
thousands of rows just to slice top 30 is pure waste.

We cap each PPI fetch to 1 page (100 rows) and source the *true*
per-source totals from the canonical UniProt entry's xref counts
(string_interaction|N, intact|N, biogrid_interaction|N). One entry call
per gene instead of up to 60 pages of pagination per source per gene.
For TP53 (14,764 STRING rows): ~150 → 3 biobtree calls for the PPI block.
No content change — page renders the same top-30 + true total count.
"""
from atlas.biobtree import map_all, entry, xref_counts
from atlas.gene.sections.base import Section

CHAINS = (
    ">>uniprot>>string_interaction",
    ">>uniprot>>intact",
    ">>uniprot>>biogrid_interaction",
    ">>uniprot>>signor",
    ">>uniprot>>esm2_similarity>>uniprot",
    ">>uniprot>>diamond_similarity>>uniprot",
)
DATASETS = ("uniprot", "string_interaction", "intact", "biogrid_interaction",
            "signor", "esm2_similarity", "diamond_similarity", "corum")

# Page cap for PPI sources. STRING/IntAct are score-desc sorted by biobtree,
# so first 100 dominates the top-30 we render with comfortable margin.
# BioGRID is not score-sorted but the rendered view makes no top-N claim.
_PPI_CAP = 1


def _f(x):
    try: return float(x)
    except (TypeError, ValueError): return 0.0


_SYM_CACHE = {}


def _uniprot_symbol(uni):
    """Resolve a UniProt accession → HGNC gene symbol (readable STRING partner),
    cached process-wide. None when unresolvable (keep the accession then)."""
    if uni in _SYM_CACHE:
        return _SYM_CACHE[uni]
    sym = None
    try:
        h = map_all(uni, ">>uniprot>>hgnc")
        if h:
            he = entry(h[0]["id"], "hgnc")
            syms = ((he.get("Attributes") or {}).get("Hgnc") or {}).get("symbols") or []
            sym = syms[0] if syms else None
    except Exception:
        pass
    _SYM_CACHE[uni] = sym
    return sym


def collect(a):
    bundle = {"section": "08_interactions", "symbol": a.symbol}
    uni = a.canonical_uniprot

    # True PPI totals from the uniprot entry's xref counts (one extra
    # entry call per gene; bounded). Avoids paginating thousands of rows
    # just to compute *_count for the render header.
    string_n = intact_n = biogrid_n = 0
    if uni:
        try:
            xc = xref_counts(entry(uni, "uniprot"))
            string_n  = xc.get("string_interaction", 0)
            intact_n  = xc.get("intact", 0)
            biogrid_n = xc.get("biogrid_interaction", 0)
        except Exception:
            pass

    # PPIs via the interaction RECORDS (NOT >>...>>uniprot, which collapses to
    # bare partner ids) — these carry the per-edge scores/evidence.
    #
    # The map projection's `uniprot_b` is a FIXED side — for ~half the rows it's
    # the query protein itself, not the partner (biobtree #34), which both looks
    # like a self-loop and loses the real partner (uniprot_a). So read each top
    # record's entry() (it carries uniprot_a + uniprot_b + score) and take the
    # non-query side as the partner. entry() is local-cache fast (~0.3ms each).
    st = map_all(uni, ">>uniprot>>string_interaction", cap=_PPI_CAP) if uni else []
    st.sort(key=lambda t: _f(t.get("score")), reverse=True)
    string = []
    for t in st[:30]:
        attrs = (((entry(t["id"], "string_interaction") or {}).get("Attributes")
                  or {}).get("StringInteraction") or {})
        ua, ub = attrs.get("uniprot_a"), attrs.get("uniprot_b")
        partner = ua if ub == uni else ub
        if not partner or partner == uni:    # skip blanks + self-interactions
            continue
        psym = _uniprot_symbol(partner)
        # Also skip same-GENE self-edges: a secondary/isoform accession that maps
        # back to this gene's own symbol (MYO18A, MIEF1) — a self-loop by gene
        # identity even though the accession differs.
        if psym and a.symbol and psym.upper() == a.symbol.upper():
            continue
        string.append({"partner": partner, "partner_symbol": psym,
                       "score": attrs.get("score") or t.get("score")})
    bundle["string"] = string
    bundle["string_count"] = string_n or len(st)

    ia = map_all(uni, ">>uniprot>>intact", cap=_PPI_CAP) if uni else []
    ia.sort(key=lambda t: _f(t.get("confidence_score")), reverse=True)
    bundle["intact"] = [{"a": t.get("protein_a_gene"), "b": t.get("protein_b_gene"),
                         "type": t.get("interaction_type"),
                         "score": t.get("confidence_score")} for t in ia[:30]]
    bundle["intact_count"] = intact_n or len(ia)
    # Full distinct IntAct (physical) partner symbols — the gene set for
    # interaction-partner ORA on the page (render-time, atlas.ora). Self-edges
    # excluded; symbols are looked up against the precomputed membership table.
    _self = (a.symbol or "").upper()
    bundle["interaction_partners"] = sorted(
        {g for t in ia for g in (t.get("protein_a_gene"), t.get("protein_b_gene"))
         if g and g.upper() != _self})

    bg = map_all(uni, ">>uniprot>>biogrid_interaction", cap=_PPI_CAP) if uni else []
    bundle["biogrid"] = [{"partner": t.get("interactor_b_symbol"),
                          "method": t.get("experimental_system")} for t in bg[:30]]
    bundle["biogrid_count"] = biogrid_n or len(bg)

    # ESM2 / Diamond similarity partners — render shows top 20 each.
    # First page (100 rows) more than covers it; deeper pages are noise.
    def partners(ds):
        return [t["id"] for t in (map_all(uni, f">>uniprot>>{ds}>>uniprot", cap=1)
                                  if uni else [])]
    bundle["esm2_similar"] = partners("esm2_similarity")
    bundle["diamond_similar"] = partners("diamond_similarity")

    # SIGNOR signaling — render shows top 30. First page covers it for any
    # cohort gene (TP53 paginates ~hundreds; bounded already cap=1).
    sig = map_all(uni, ">>uniprot>>signor", cap=1) if uni else []
    bundle["signor"] = [{"a": t.get("entity_a"), "b": t.get("entity_b"),
                         "effect": t.get("effect"), "mechanism": t.get("mechanism")} for t in sig]
    bundle["signor_count"] = len(sig)

    # CORUM — named, experimentally-characterized protein complexes this protein
    # is a subunit of. A curated complement to the score-ranked STRING/IntAct
    # firehose: a discrete, biologically-defined complex membership. Schema:
    # id|name|organism|subunit_count|subunit_genes|has_drug_targets. Human only;
    # largest complexes first.
    cor = [t for t in (map_all(uni, ">>uniprot>>corum") if uni else [])
           if (t.get("organism") or "Human").lower() == "human"]
    cor.sort(key=lambda t: int(t.get("subunit_count") or 0), reverse=True)
    bundle["corum"] = [{"id": t.get("id"), "name": t.get("name"),
                        "subunit_count": t.get("subunit_count"),
                        "subunits": t.get("subunit_genes"),
                        "has_drug_targets": t.get("has_drug_targets") == "true"}
                       for t in cor]
    bundle["corum_count"] = len(cor)
    return bundle

SECTION = Section(
    id="8", name="interactions",
    description="Protein-protein interactions (STRING/IntAct/BioGRID with scores), SIGNOR signaling, ESM2/Diamond structural similarity",
    needs=("canonical_uniprot",),
    produces=("string", "intact", "biogrid", "signor", "esm2_similar", "diamond_similar",
              "interaction_partners", "corum", "corum_count"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
