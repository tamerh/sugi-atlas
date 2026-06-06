"""§8 — target pathways. Reuses the gene §7 pathways collector over the drug's
target genes (cohort fan) and aggregates: Reactome pathways shared across the
target set (ranked by how many targets touch each) + dominant GO
biological-process terms. Shows the biology the drug's targets sit in —
e.g. Imatinib's ABL1/DDR1/DDR2 → kinase-signaling pathways."""
from collections import Counter
from atlas.drug import cohort
from atlas.gene.sections import s07_pathways
from atlas.section import Section


def collect(a):
    gas = cohort.target_gene_anchors(a)
    bundles = cohort.over_targets(s07_pathways.collect, gas)

    # Reactome: pathway id -> {name, set(gene symbols touching it)}
    pw = {}
    for b in bundles:
        for p in (b.get("reactome") or []):
            e = pw.setdefault(p["id"], {"name": p.get("name"), "genes": set()})
            if p.get("name"):
                e["name"] = p["name"]
            e["genes"].add(b.get("symbol"))

    # GO biological_process: term -> count of targets annotated with it
    go_bp = Counter()
    go_names = {}
    for b in bundles:
        for t in ((b.get("go") or {}).get("biological_process") or []):
            go_bp[t["id"]] += 1
            go_names[t["id"]] = t.get("name")

    top_pathways = sorted(pw.items(), key=lambda kv: (-len(kv[1]["genes"]), kv[0]))
    return {
        "section": "08_target_pathways",
        "target_genes": [ga.symbol for ga in gas],
        "pathway_count": len(pw),
        "top_pathways": [{"id": pid, "name": e["name"],
                          "gene_count": len(e["genes"]),
                          "genes": sorted(g for g in e["genes"] if g)}
                         for pid, e in top_pathways[:50]],
        "top_go_bp": [{"id": gid, "name": go_names.get(gid), "target_count": n}
                      for gid, n in go_bp.most_common(15)],
    }


SECTION = Section(
    id="8", name="target_pathways",
    description=("Reactome + GO biological-process aggregation over the drug's "
                 "target genes (reuses gene §7 pathways collector via fan-out)"),
    needs=("targets",),
    produces=("target_genes", "pathway_count", "top_pathways", "top_go_bp"),
    datasets=("reactome", "go", "uniprot", "ensembl", "hgnc"),
    chains=(">>uniprot>>reactome", ">>uniprot>>go"),
    collect_fn=collect,
)
