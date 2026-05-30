"""§14 — pathway analysis: per-cohort-gene Reactome pathways + cohort-level
aggregation (which pathways are most over-represented?). REUSE wrapper over
gene §7."""
from collections import Counter, defaultdict
from atlas.section import Section
from atlas.disease.cohort import fan
from atlas.gene.sections import s07_pathways

CHAINS   = (">>hgnc>>reactome", ">>hgnc>>ensembl>>reactome")
DATASETS = ("hgnc", "ensembl", "reactome")

TOP_N = 30


def collect(a):
    bundles = fan(s07_pathways.SECTION.collect_fn, a.cohort)

    per_gene_pathways = []
    pathway_gene_counts: Counter = Counter()
    pathway_name: dict = {}
    pathway_to_genes: dict = defaultdict(list)

    for b in bundles:
        sym = b.get("symbol")
        hgnc_id = b.get("hgnc_id")
        reactome = b.get("reactome") or []
        ids, names = [], []
        for p in reactome:
            pid = p.get("id")
            if not pid:
                continue
            pname = p.get("name")
            ids.append(pid)
            names.append(pname)
            if pname and pid not in pathway_name:
                pathway_name[pid] = pname
            pathway_gene_counts[pid] += 1
            if sym and sym not in pathway_to_genes[pid]:
                pathway_to_genes[pid].append(sym)
        per_gene_pathways.append({
            "symbol": sym,
            "hgnc_id": hgnc_id,
            "reactome_ids": ids,
            "reactome_names": names,
        })

    top_pathways = []
    for pid, count in pathway_gene_counts.most_common(TOP_N):
        top_pathways.append({
            "id": pid,
            "name": pathway_name.get(pid),
            "gene_count": count,
            "gene_symbols": list(pathway_to_genes[pid]),
        })

    return {
        "section": "14_pathways",
        "mondo_id": a.mondo_id,
        "per_gene_pathways": per_gene_pathways,
        "pathway_gene_counts": pathway_gene_counts,
        "top_pathways": top_pathways,
        "pathway_count": len(pathway_gene_counts),
    }


SECTION = Section(
    id="14", name="pathways",
    description=("Cohort-level Reactome pathway aggregation: top enriched "
                 "pathways, per-pathway gene counts, druggability per pathway."),
    needs=("cohort",),
    produces=("per_gene_pathways", "top_pathways", "pathway_gene_counts",
              "pathway_count"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
