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

# GO ECO evidence (now in the lite >>uniprot>>go projection). The experimental
# subtree = directly-assayed / mutant-phenotype / interaction evidence — the
# high-confidence annotations, vs phylogenetic/computational/electronic (IEA).
_ECO_EXPERIMENTAL = {
    "ECO:0000314", "ECO:0000315", "ECO:0000316", "ECO:0000353", "ECO:0000270",
    "ECO:0000269", "ECO:0006056", "ECO:0007005", "ECO:0007007", "ECO:0007003",
}
_ECO_LABEL = {
    "ECO:0000314": "direct assay", "ECO:0000315": "mutant phenotype",
    "ECO:0000316": "genetic interaction", "ECO:0000353": "physical interaction",
    "ECO:0000270": "expression pattern", "ECO:0000269": "experiment",
    "ECO:0000304": "traceable author", "ECO:0000303": "author statement",
    "ECO:0000250": "sequence similarity", "ECO:0000266": "sequence orthology",
    "ECO:0000318": "phylogenetic", "ECO:0000247": "sequence alignment",
    "ECO:0000501": "electronic (IEA)", "ECO:0007669": "electronic (IEA)",
    "ECO:0000305": "curator inference", "ECO:0000307": "no biological data",
}
_NS_ABBR = {"biological_process": "BP", "molecular_function": "MF",
            "cellular_component": "CC"}


def _eco_tier(eco):
    """Coarse evidence tier for the per-gene GO summary."""
    if eco in _ECO_EXPERIMENTAL:
        return "experimental"
    if eco in ("ECO:0000501", "ECO:0007669"):
        return "electronic"
    if eco in ("ECO:0000318", "ECO:0000250", "ECO:0000266", "ECO:0000247"):
        return "computational"
    if eco in ("ECO:0000303", "ECO:0000304", "ECO:0000305"):
        return "author/curator"
    return "other"

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
    # (CDKN2A p14ARF mitophagy terms absent from p16/ensembl). UniProt rows carry
    # the ECO evidence code (ensembl rows don't), so UniProt wins on shared terms
    # — ensembl only ADDS terms UniProt lacks — to preserve evidence.
    go_map = {}
    for u in a.reviewed_uniprots:
        for t in map_all(u, ">>uniprot>>go"):
            go_map[t["id"]] = t
    for t in (map_all(a.ensembl_id, ">>ensembl>>go") if a.ensembl_id else []):
        go_map.setdefault(t["id"], t)
    grouped = {"biological_process": [], "molecular_function": [], "cellular_component": []}
    experimental = []
    ev_counts = Counter()
    for t in go_map.values():
        eco = t.get("evidence")
        ev_counts[_eco_tier(eco)] += 1
        row = {"id": t["id"], "name": t.get("name"), "evidence": eco}
        grouped.setdefault(t.get("type"), []).append(row)
        if eco in _ECO_EXPERIMENTAL:
            experimental.append({"id": t["id"], "name": t.get("name"),
                                 "namespace": _NS_ABBR.get(t.get("type"), t.get("type")),
                                 "evidence_label": _ECO_LABEL.get(eco, eco)})
    experimental.sort(key=lambda r: (r["namespace"], r["name"] or r["id"]))
    bundle["go"] = grouped
    bundle["go_counts"] = {k: len(v) for k, v in grouped.items()}
    bundle["go_experimental"] = experimental
    bundle["go_experimental_count"] = len(experimental)
    bundle["go_total"] = len(go_map)

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
              "go_experimental", "go_experimental_count", "go_total",
              "go_parent_rollup", "reactome_parent_rollup"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
