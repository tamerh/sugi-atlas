"""§7 — pathways: Reactome (union over all reviewed UniProt + ensembl route),
MSigDB gene sets, GO terms (BP/MF/CC unioned UniProt-GOA + Ensembl annotation),
plus top-level GO + Reactome parent rollups for hierarchical navigation."""
from collections import Counter
from atlas.biobtree import map_all, xref_counts
from atlas.gene.sections.base import Section

CHAINS = (
    ">>uniprot>>reactome", ">>ensembl>>reactome",
    ">>hgnc>>msigdb",
    ">>uniprot>>go", ">>ensembl>>go",
    ">>go>>goparent", ">>reactome>>reactomeparent",
)
DATASETS = ("reactome", "msigdb", "go", "uniprot", "ensembl", "hgnc")

# Top-N terms to expand into parent rollups. Bounded so we don't pay a
# parent call per GO term on TP53-style high-GO-count genes — the top-N
# already captures the dominant categories.
_PARENT_ROLLUP_TOP_N = 20

def collect(a):
    bundle = {"section": "07_pathways", "symbol": a.symbol}
    xc = xref_counts(a.hgnc_entry)

    # Reactome: union all reviewed uniprots + ensembl gene-level route.
    # Dual-product genes (CDKN2A p16+p14ARF) and pathways the canonical uniprot
    # alone misses both land here.
    rx = {}
    for u in a.reviewed_uniprots:
        for t in map_all(u, ">>uniprot>>reactome"):
            rx[t["id"]] = t.get("name") or rx.get(t["id"])
    for t in (map_all(a.ensembl_id, ">>ensembl>>reactome") if a.ensembl_id else []):
        rx[t["id"]] = t.get("name") or rx.get(t["id"])
    bundle["reactome"] = [{"id": k, "name": v} for k, v in rx.items()]
    bundle["reactome_count"] = len(bundle["reactome"])

    msig = map_all(a.hgnc_id, ">>hgnc>>msigdb")
    bundle["msigdb"] = [{"id": t["id"], "name": t.get("standard_name"),
                         "collection": t.get("collection")} for t in msig]
    bundle["msigdb_total"] = xc.get("msigdb", len(msig))

    # GO: UniProt-GOA + Ensembl annotation diverge ~20%; per-product GO differs
    # (CDKN2A p14ARF mitophagy terms absent from p16/ensembl).
    go_map = {}
    for u in a.reviewed_uniprots:
        for t in map_all(u, ">>uniprot>>go"):
            go_map[t["id"]] = t
    for t in (map_all(a.ensembl_id, ">>ensembl>>go") if a.ensembl_id else []):
        go_map[t["id"]] = t
    grouped = {"biological_process": [], "molecular_function": [], "cellular_component": []}
    for t in go_map.values():
        grouped.setdefault(t.get("type"), []).append({"id": t["id"], "name": t.get("name")})
    bundle["go"] = grouped
    bundle["go_counts"] = {k: len(v) for k, v in grouped.items()}

    # GO parent rollup — for each of the top-N GO terms in each category,
    # fetch its `goparent` to map fine-grained terms to L1/L2 categories.
    # Counts give the page a "what high-level process dominates" view.
    go_parent_counts = Counter()
    go_parent_names = {}
    for ns, terms in grouped.items():
        for t in terms[:_PARENT_ROLLUP_TOP_N]:
            for p in map_all(t["id"], ">>go>>goparent"):
                pid = p.get("id")
                if not pid:
                    continue
                go_parent_counts[pid] += 1
                go_parent_names[pid] = p.get("name")
    bundle["go_parent_rollup"] = [
        {"id": pid, "name": go_parent_names.get(pid), "term_count": n}
        for pid, n in go_parent_counts.most_common()
    ]

    # Reactome parent rollup — same idea, but Reactome's hierarchy is
    # tighter; one parent per pathway is the norm.
    reactome_parent_counts = Counter()
    reactome_parent_names = {}
    for rid in list(rx.keys())[:_PARENT_ROLLUP_TOP_N]:
        for p in map_all(rid, ">>reactome>>reactomeparent"):
            pid = p.get("id")
            if not pid:
                continue
            reactome_parent_counts[pid] += 1
            reactome_parent_names[pid] = p.get("name")
    bundle["reactome_parent_rollup"] = [
        {"id": pid, "name": reactome_parent_names.get(pid), "pathway_count": n}
        for pid, n in reactome_parent_counts.most_common()
    ]

    return bundle

SECTION = Section(
    id="7", name="pathways",
    description="Reactome pathways (union over all reviewed uniprots + ensembl), MSigDB gene sets, GO terms (BP/MF/CC)",
    needs=("hgnc_id", "hgnc_entry", "ensembl_id", "reviewed_uniprots"),
    produces=("reactome", "msigdb", "go", "go_counts",
              "go_parent_rollup", "reactome_parent_rollup"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
