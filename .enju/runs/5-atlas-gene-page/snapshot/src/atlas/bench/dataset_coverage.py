#!/usr/bin/env python3
"""Map biobtree's queryable datasets against what the gene collector uses.

Source of truth = the MCP map EDGES (`/api/help?topic=edges`) — the authoritative
graph of what-connects-to-what. (Do NOT use /api/meta: it lists 500+ derived
xref-identifier types, e.g. a "cosmic" id that is NOT integrated somatic data.)

Re-run to see coverage and to spot NEWLY ADDED datasets: incoming somatic
(civic / intogen / + one more, in progress) will appear here as new edge nodes
once integrated.

  python3 biobtree/collector/dataset_coverage.py
"""
import os, re, json, urllib.request

API = "http://127.0.0.1:8000/api"
HERE = os.path.dirname(os.path.abspath(__file__))

def edge_graph():
    d = json.load(urllib.request.urlopen(f"{API}/help?topic=edges", timeout=20))
    g = {}
    for line in d.get("edges", "").splitlines():
        m = re.match(r"^([a-z0-9_]+):\s*(.+)$", line.strip())
        if m:
            # targets may carry inline "# comment" annotations — strip them
            targets = [t.split("#")[0].strip() for t in m.group(2).split(",")]
            g[m.group(1)] = [t for t in targets if re.fullmatch(r"[a-z0-9_]+", t)]
    return g

def main():
    g = edge_graph()
    nodes = set(g) | {t for ts in g.values() for t in ts}
    src = open(os.path.join(HERE, "collect.py")).read()
    used = ({d for d in re.findall(r">>([a-z0-9_]+)", src)} |
            {d for d in re.findall(r'"([a-z0-9_]+)"', src)}) & nodes
    uncov = sorted(nodes - used)
    print(f"queryable datasets (edge nodes): {len(nodes)}  |  used by collector: {len(used)}\n")
    print(f"USED ({len(used)}):\n  " + ", ".join(sorted(used)) + "\n")
    print(f"UNCOVERED ({len(uncov)}):\n  " + ", ".join(uncov))

if __name__ == "__main__":
    main()
