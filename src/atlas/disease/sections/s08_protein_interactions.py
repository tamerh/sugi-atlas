"""§8 — protein interactions among cohort genes (STRING / IntAct / BioGRID).
REUSE wrapper over gene §8 + inter-cohort edge graph.

Key value-add: surface UNDRUGGED cohort genes that interact with DRUGGED ones
(indirect druggability path — §17 feeds off this).

NOTE: undrugged_to_drugged_bridges stays empty here — it requires the §10
(drug_targets) bundle, which isn't a dependency of §8. The render-only §17
joins §8's `per_gene_interactions` with §10's drugged-target set at compose
time and fills the bridge list."""
from collections import Counter

from atlas.section import Section
from atlas.disease.cohort import fan
from atlas.gene.sections import s08_interactions

CHAINS   = (">>uniprot>>string_interaction", ">>uniprot>>intact",
            ">>uniprot>>biogrid_interaction")
DATASETS = ("uniprot", "string_interaction", "intact", "biogrid_interaction")


def _f(x):
    try: return float(x)
    except (TypeError, ValueError): return 0.0


def _canonical_edge(a, b):
    """Undirected edge in alphabetical order."""
    return (a, b) if a <= b else (b, a)


def collect(a):
    g8 = fan(s08_interactions.SECTION.collect_fn, a.cohort)

    cohort_symbols = {ga.symbol for ga in a.cohort}
    # Map UniProt id -> cohort symbol (for STRING partners, which are uniprot ids).
    uni_to_sym = {ga.canonical_uniprot: ga.symbol for ga in a.cohort
                  if ga.canonical_uniprot}

    per_gene_interactions = []
    # edge -> set of source dataset names
    edge_sources: dict = {}

    for b in g8:
        sym = b.get("symbol")
        if not sym or b.get("_error"):
            per_gene_interactions.append({
                "symbol": sym, "hgnc_id": b.get("hgnc_id"),
                "interactor_count": 0, "top_interactors": [],
            })
            continue

        string_list = b.get("string", []) or []
        intact_list = b.get("intact", []) or []
        biogrid_list = b.get("biogrid", []) or []

        # Collect partner symbols + per-source counts. STRING partners are
        # uniprot ids — resolve to cohort symbols when possible; otherwise
        # keep the uniprot id so the count isn't zero for non-cohort hits.
        all_partners = []  # ordered, with potential dupes; gene §8 already sorts string by score
        for t in string_list:
            partner_uni = t.get("partner")
            partner_sym = uni_to_sym.get(partner_uni, partner_uni)
            if partner_sym:
                all_partners.append(partner_sym)
                if partner_sym in cohort_symbols and partner_sym != sym:
                    e = _canonical_edge(sym, partner_sym)
                    edge_sources.setdefault(e, set()).add("string_interaction")

        for t in intact_list:
            pa, pb = t.get("a"), t.get("b")
            # IntAct rows give gene symbols on both ends; pick the partner.
            partner = pb if pa == sym else pa
            if partner:
                all_partners.append(partner)
                if partner in cohort_symbols and partner != sym:
                    e = _canonical_edge(sym, partner)
                    edge_sources.setdefault(e, set()).add("intact")

        for t in biogrid_list:
            partner = t.get("partner")
            if partner:
                all_partners.append(partner)
                if partner in cohort_symbols and partner != sym:
                    e = _canonical_edge(sym, partner)
                    edge_sources.setdefault(e, set()).add("biogrid_interaction")

        # Dedup partner list preserving order (STRING-by-score wins ties).
        seen = set()
        unique_partners = []
        for p in all_partners:
            if p not in seen and p != sym:
                seen.add(p)
                unique_partners.append(p)

        interactor_count = (b.get("string_count", 0) + b.get("intact_count", 0)
                            + b.get("biogrid_count", 0))
        per_gene_interactions.append({
            "symbol": sym,
            "hgnc_id": b.get("hgnc_id"),
            "interactor_count": interactor_count,
            "top_interactors": unique_partners[:5],
        })

    cohort_edges = [
        {"a": e[0], "b": e[1], "sources": sorted(srcs)}
        for e, srcs in sorted(edge_sources.items())
    ]

    hub_genes = sorted(per_gene_interactions,
                       key=lambda g: g["interactor_count"], reverse=True)[:10]
    hub_genes = [{"symbol": g["symbol"],
                  "interactor_count": g["interactor_count"]} for g in hub_genes]

    return {
        "section": "08_protein_interactions",
        "mondo_id": a.mondo_id,
        "per_gene_interactions": per_gene_interactions,
        "cohort_edges": cohort_edges,
        "cohort_edge_count": len(cohort_edges),
        "hub_genes": hub_genes,
        # Placeholder — §17 render-only join with §10 drug_targets fills this.
        "undrugged_to_drugged_bridges": [],
    }

SECTION = Section(
    id="8", name="protein_interactions",
    description=("Cohort-internal interaction graph (STRING / IntAct / BioGRID): "
                 "hub genes, undrugged-↔-drugged bridges for indirect "
                 "druggability."),
    needs=("cohort",),
    produces=("per_gene_interactions", "cohort_edges", "hub_genes",
              "undrugged_to_drugged_bridges"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
