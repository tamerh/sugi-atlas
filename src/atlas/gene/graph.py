#!/usr/bin/env python3
"""Generate Mermaid diagrams of Atlas's gene pipeline from section metadata.

Three views — pick with --view:

  overview    Anchors → 12 sections → render → summary → publish (the outer DAG).
  datasets    Biobtree dataset graph traversed by Atlas: nodes are datasets,
              edges are the multi-hops declared by sections.
  section     Multi-hop walks for one section's chains (--section <id>).

Output is Mermaid (markdown-fenced or raw with --raw). Paste into GitHub,
mermaid.live, the preprint figure, anywhere mermaid renders.

  python -m atlas.gene.graph --view overview
  python -m atlas.gene.graph --view datasets > docs/atlas_dataset_graph.md
  python -m atlas.gene.graph --view section --section 12
"""
import argparse, re
from atlas.gene.sections import REGISTRY

def _hops(chain):
    """Parse '>>hgnc>>clinvar>>mondo>>gwas_study' -> ['hgnc','clinvar','mondo','gwas_study'].
    Filters like clinvar[...] are stripped; tokens are the dataset names only."""
    parts = [p for p in chain.split(">>") if p]
    return [re.sub(r"\[.*", "", p) for p in parts]

def _overview():
    L = ["flowchart TD",
         "  symbol([symbol]) --> resolve[resolve anchors]",
         "  resolve --> A[(hgnc · ensembl · uniprots · canonical_tx)]"]
    for sid in REGISTRY:
        s = REGISTRY[sid]
        L.append(f"  A --> S{sid}[§{sid} {s.name}]")
        L.append(f"  S{sid} --> R{sid}[render §{sid}]")
        L.append(f"  R{sid} --> body[(body.md)]")
    L += [
        "  body --> summary[LLM summary (Qwen3)]",
        "  body --> gate[body_gate]",
        "  summary --> page[(page.md)]",
        "  gate --> page",
        "  page --> publish[publish to dist]",
    ]
    return "\n".join(L)

def _datasets():
    """Dataset-level multi-hop graph: union of every chain's hops across sections.
    Edge labels include the section ids that traverse that hop."""
    edges = {}  # (src, dst) -> set(section_ids)
    for sid, s in REGISTRY.items():
        for ch in s.chains:
            hops = _hops(ch)
            for i in range(len(hops) - 1):
                edges.setdefault((hops[i], hops[i + 1]), set()).add(sid)
    nodes = sorted({n for e in edges for n in e})
    L = ["flowchart LR"]
    for n in nodes:
        L.append(f"  {n}({n})")
    for (src, dst), secs in sorted(edges.items()):
        label = "§" + ",".join(sorted(secs, key=int))
        L.append(f"  {src} -- {label} --> {dst}")
    return "\n".join(L)

def _section(sid):
    s = REGISTRY[sid]
    L = [f"flowchart LR", f"  %% §{sid} {s.name} — {s.description}"]
    seen_edges = set()
    for ch in s.chains:
        hops = _hops(ch)
        for i in range(len(hops) - 1):
            e = (hops[i], hops[i + 1])
            if e in seen_edges: continue
            seen_edges.add(e)
            L.append(f"  {e[0]} --> {e[1]}")
    return "\n".join(L)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--view", choices=("overview", "datasets", "section"), default="overview")
    ap.add_argument("--section", help="section id (1-12), required for --view section")
    ap.add_argument("--raw", action="store_true", help="emit raw mermaid (no ```mermaid fence)")
    args = ap.parse_args()
    if args.view == "overview":
        body = _overview()
    elif args.view == "datasets":
        body = _datasets()
    elif args.view == "section":
        if not args.section:
            ap.error("--view section requires --section <id>")
        body = _section(args.section)
    if args.raw:
        print(body)
    else:
        print("```mermaid"); print(body); print("```")

if __name__ == "__main__":
    main()
