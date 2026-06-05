"""§14 — pathway / GO enrichment: which Reactome pathways and GO biological
processes are most OVER-REPRESENTED across the disease's associated genes.

Aggregate section over the WIDE enrichment cohort (~250 evidence-ranked genes).
Ranking is real over-representation analysis (ORA): a hypergeometric test of the
cohort's overlap against a genome-wide background (atlas.ora + the precomputed
data/background/{reactome,go}.json), Benjamini-Hochberg FDR, sorted by enrichment
— so size-biased umbrella categories no longer float to the top by raw count.
Raw counts + gene members are kept (consumer ground-truth); enrichment is added."""
from collections import Counter, defaultdict
from atlas.section import Section
from atlas.disease.cohort import enrichment_fan
from atlas.ora import enrich, background

CHAINS   = (">>hgnc>>ensembl>>reactome", ">>hgnc>>uniprot>>go")
DATASETS = ("hgnc", "ensembl", "reactome", "uniprot", "go")

TOP_N = 30


def _rank(it):                           # tested (lowest FDR) first; else by count
    tested = it["fdr"] is not None
    return (0 if tested else 1, it["fdr"] if tested else 0.0,
            -(it["fold"] or 0.0), -it["k"], it["id"])


def _tally(rows_per_gene, predicate=None):
    """Count a fanned chain into (counts, names, to_genes, genes_with_any) — one
    count per (gene, category)."""
    counts: Counter = Counter()
    names: dict = {}
    to_genes: dict = defaultdict(list)
    with_any = 0
    for _hgnc, sym, rows in rows_per_gene:
        rows = [r for r in rows if predicate is None or predicate(r)]
        if rows:
            with_any += 1
        seen = set()
        for r in rows:
            cid = r.get("id")
            if not cid or cid in seen:
                continue
            seen.add(cid)
            counts[cid] += 1
            if r.get("name") and cid not in names:
                names[cid] = r.get("name")
            if sym and sym not in to_genes[cid]:
                to_genes[cid].append(sym)
    return counts, names, to_genes, with_any


def _ora_table(counts, names, to_genes, cohort_n, bg_name):
    universe_n, sizes = background(bg_name)
    items = enrich([{"id": cid, "name": names.get(cid), "k": cnt,
                     "K": sizes.get(cid, 0), "gene_symbols": to_genes[cid]}
                    for cid, cnt in counts.items()],
                   cohort_n=cohort_n, universe_n=universe_n)
    items.sort(key=_rank)
    return [{"id": it["id"], "name": it["name"], "gene_count": it["k"],
             "fold": it["fold"], "fdr": it["fdr"], "gene_symbols": it["gene_symbols"]}
            for it in items[:TOP_N]]


def collect(a):
    rc, rn, rg, genes_with_pathways = _tally(
        enrichment_fan(a.enrichment_cohort, ">>hgnc>>ensembl>>reactome"))
    top_pathways = _ora_table(rc, rn, rg, genes_with_pathways, "reactome")

    gc, gn, gg, genes_with_go = _tally(
        enrichment_fan(a.enrichment_cohort, ">>hgnc>>uniprot>>go"),
        predicate=lambda r: r.get("type") == "biological_process")
    top_go = _ora_table(gc, gn, gg, genes_with_go, "go")

    return {
        "section": "14_pathways",
        "mondo_id": a.mondo_id,
        "top_pathways": top_pathways,
        "pathway_gene_counts": rc,
        "pathway_count": len(rc),
        "top_go": top_go,
        "go_count": len(gc),
        "genes_with_go": genes_with_go,
        # enrichment-cohort provenance — the breadth the stats were computed over
        "enrichment_size": len(a.enrichment_cohort),
        "genes_with_pathways": genes_with_pathways,
    }


SECTION = Section(
    id="14", name="pathways",
    description=("Cohort-level Reactome pathway + GO biological-process "
                 "over-representation analysis (hypergeometric vs a genome-wide "
                 "background, BH-FDR) over the wide enrichment cohort."),
    needs=("enrichment_cohort",),
    produces=("top_pathways", "pathway_gene_counts", "pathway_count",
              "top_go", "go_count", "genes_with_go",
              "enrichment_size", "genes_with_pathways"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
