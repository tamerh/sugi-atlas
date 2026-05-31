#!/usr/bin/env python3
"""Map biobtree's queryable datasets against what Atlas collectors use.

Source of truth for the queryable graph: biobtree's `/api/help?topic=edges`
— the authoritative graph of which dataset connects to which. (NOT
`/api/meta`: that lists 500+ derived xref-identifier types, including many
that aren't real first-class datasets.)

Source of truth for Atlas usage: each entity's section REGISTRY. Every
Section dataclass declares its `datasets` and `chains` metadata, so we union
those (no source-text scraping — used to break whenever a section file
moved).

  python3 -m atlas.bench.dataset_coverage              # gene + disease
  python3 -m atlas.bench.dataset_coverage --entity gene
  python3 -m atlas.bench.dataset_coverage --entity disease
  python3 -m atlas.bench.dataset_coverage --json       # machine-readable
"""
import argparse, json, re, urllib.request

API = "http://127.0.0.1:8000/api"


def edge_graph():
    """Parse biobtree's edges help into {src_dataset: [target_datasets, ...]}.
    Returns also the flat set of all node names."""
    d = json.load(urllib.request.urlopen(f"{API}/help?topic=edges", timeout=20))
    g = {}
    for line in d.get("edges", "").splitlines():
        m = re.match(r"^([a-z0-9_]+):\s*(.+)$", line.strip())
        if m:
            # targets may carry inline "# comment" annotations — strip them.
            targets = [t.split("#")[0].strip() for t in m.group(2).split(",")]
            g[m.group(1)] = [t for t in targets if re.fullmatch(r"[a-z0-9_]+", t)]
    return g


def datasets_used_by_registry(registry):
    """Union of every section's declared `datasets` + chain-extracted nodes.

    chains are strings like '>>uniprot>>interpro' — we extract the lowercase
    identifier tokens. Also catches any filter-suffix removal (e.g.
    `chembl_molecule[highestDevelopmentPhase>=1]` → chembl_molecule)."""
    out = set()
    for sec in registry.values():
        out.update(sec.datasets or ())
        for chain in (sec.chains or ()):
            for token in re.findall(r"([a-z0-9_]+)", chain or ""):
                out.add(token)
    return out


def _coverage_for(entity_name, registry, edge_nodes):
    used = datasets_used_by_registry(registry) & edge_nodes
    uncov = sorted(edge_nodes - used)
    return {
        "entity": entity_name,
        "edge_nodes_total": len(edge_nodes),
        "used_count": len(used),
        "used": sorted(used),
        "uncovered_count": len(uncov),
        "uncovered": uncov,
    }


def gather():
    g = edge_graph()
    edge_nodes = set(g) | {t for ts in g.values() for t in ts}
    from atlas.gene.sections import REGISTRY as GENE_REG
    from atlas.disease.sections import REGISTRY as DISEASE_REG
    return {
        "gene":    _coverage_for("gene",    GENE_REG,    edge_nodes),
        "disease": _coverage_for("disease", DISEASE_REG, edge_nodes),
        "edge_nodes_total": len(edge_nodes),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--entity", choices=["gene", "disease", "both"], default="both")
    ap.add_argument("--json", action="store_true", help="emit machine-readable JSON")
    args = ap.parse_args()
    data = gather()
    if args.json:
        print(json.dumps(data, indent=2))
        return

    entities = ("gene", "disease") if args.entity == "both" else (args.entity,)
    print(f"queryable datasets (edge nodes): {data['edge_nodes_total']}\n")
    for ent in entities:
        c = data[ent]
        print(f"=== {ent.upper()} entity ===")
        print(f"used: {c['used_count']} / {c['edge_nodes_total']}\n")
        print(f"USED ({c['used_count']}):\n  " + ", ".join(c['used']) + "\n")
        print(f"UNCOVERED ({c['uncovered_count']}):\n  " + ", ".join(c['uncovered']) + "\n")
    if args.entity == "both":
        # Diff: which uncovered-by-gene are now covered-by-disease (or vice versa)?
        gene_used = set(data["gene"]["used"])
        dis_used = set(data["disease"]["used"])
        union_uncov = sorted((set(data["gene"]["uncovered"]) &
                              set(data["disease"]["uncovered"])))
        print(f"NEVER USED BY EITHER ({len(union_uncov)}):\n  " + ", ".join(union_uncov))


if __name__ == "__main__":
    main()
