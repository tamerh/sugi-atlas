"""§14 — pathway enrichment: which Reactome pathways are most over-represented
across the disease's associated genes.

This is an AGGREGATE section, so it runs over the WIDE enrichment cohort
(~250 evidence-ranked genes) rather than the small display cohort — pathway
enrichment is far sharper over 250 genes than 75. It fans ONE cheap chain
(`>>hgnc>>ensembl>>reactome`) per gene, not the full gene §7 plan, so wide
stays cheap. (The direct `>>hgnc>>reactome` edge is empty; the gene route is
via Ensembl.)"""
from collections import Counter, defaultdict
from atlas.section import Section
from atlas.disease.cohort import enrichment_fan
from atlas.ora import enrich, reactome_background

CHAINS   = (">>hgnc>>ensembl>>reactome",)
DATASETS = ("hgnc", "ensembl", "reactome")

TOP_N = 30


def collect(a):
    rows_per_gene = enrichment_fan(a.enrichment_cohort, ">>hgnc>>ensembl>>reactome")

    pathway_gene_counts: Counter = Counter()
    pathway_name: dict = {}
    pathway_to_genes: dict = defaultdict(list)
    genes_with_pathways = 0

    for _hgnc, sym, reactome in rows_per_gene:
        if reactome:
            genes_with_pathways += 1
        seen = set()
        for p in reactome:
            pid = p.get("id")
            if not pid or pid in seen:           # one count per (gene, pathway)
                continue
            seen.add(pid)
            pathway_gene_counts[pid] += 1
            if p.get("name") and pid not in pathway_name:
                pathway_name[pid] = p.get("name")
            if sym and sym not in pathway_to_genes[pid]:
                pathway_to_genes[pid].append(sym)

    # Over-representation analysis vs the genome-wide Reactome background: rank by
    # statistical enrichment (BH-FDR, then fold), NOT raw count — so size-biased
    # umbrella pathways (Signal Transduction, Disease) stop floating to the top and
    # the disease-specific pathways surface. Counts + gene lists are KEPT (the
    # consumer, human or agent, cites them as ground-truth); enrichment is added.
    universe_n, sizes = reactome_background()
    items = enrich([{"id": pid, "name": pathway_name.get(pid), "k": cnt,
                     "K": sizes.get(pid, 0), "gene_symbols": pathway_to_genes[pid]}
                    for pid, cnt in pathway_gene_counts.items()],
                   cohort_n=genes_with_pathways, universe_n=universe_n)

    def _rank(it):                       # tested (lowest FDR) first; else by count
        tested = it["fdr"] is not None
        return (0 if tested else 1, it["fdr"] if tested else 0.0,
                -(it["fold"] or 0.0), -it["k"], it["id"])
    items.sort(key=_rank)
    top_pathways = [{"id": it["id"], "name": it["name"], "gene_count": it["k"],
                     "background_size": it["K"], "fold": it["fold"], "fdr": it["fdr"],
                     "gene_symbols": it["gene_symbols"]}
                    for it in items[:TOP_N]]

    return {
        "section": "14_pathways",
        "mondo_id": a.mondo_id,
        "top_pathways": top_pathways,
        "pathway_gene_counts": pathway_gene_counts,
        "pathway_count": len(pathway_gene_counts),
        # enrichment-cohort provenance — the breadth the stats were computed over
        "enrichment_size": len(a.enrichment_cohort),
        "genes_with_pathways": genes_with_pathways,
    }


SECTION = Section(
    id="14", name="pathways",
    description=("Cohort-level Reactome pathway enrichment over the wide "
                 "enrichment cohort: top over-represented pathways by the "
                 "number of associated genes touching each."),
    needs=("enrichment_cohort",),
    produces=("top_pathways", "pathway_gene_counts", "pathway_count",
              "enrichment_size", "genes_with_pathways"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
