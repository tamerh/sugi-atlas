#!/usr/bin/env python3
"""Derive the canonical per-section query plan from committed *_calls.json.

Across all genes that have per-section data, group calls by SECTION and reduce
each to a template (tool + chain/source, inputs dropped). Report per section:
the chains/sources used, and coverage (how many genes used each) -> the stable
plan vs per-gene variation.
"""
import json, glob, re, collections, os

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BASE = os.path.join(REPO, "biobtree/gene/sections")

def tmpl(call):
    tool = call.get("tool","?").replace("biobtree_","")
    a = call.get("args",{})
    if tool == "map":
        chain = a.get("chain","")
        return f"map {chain}"
    if tool == "entry":
        return f"entry {a.get('dataset','?')}"
    if tool == "search":
        return f"search {a.get('dataset','(default)')}"
    return tool

# section files look like  <N>_<name>_calls.json
sec_calls = collections.defaultdict(lambda: collections.Counter())  # section -> tmpl -> gene-coverage
sec_total = collections.defaultdict(lambda: collections.Counter())  # section -> tmpl -> total count
sec_genes = collections.defaultdict(set)

genes = [d for d in os.listdir(BASE) if os.path.isdir(os.path.join(BASE,d))]
for g in genes:
    for f in glob.glob(os.path.join(BASE,g,"*_calls.json")):
        base = os.path.basename(f)
        if base == "_calls.json": continue
        m = re.match(r"(\d+)_(.+)_calls\.json", base)
        if not m: continue
        section = f"{int(m.group(1)):02d}_{m.group(2)}"
        try: calls = json.load(open(f))
        except Exception: continue
        sec_genes[section].add(g)
        seen=set()
        for c in calls:
            t=tmpl(c); sec_total[section][t]+=1
            if t not in seen: sec_calls[section][t]+=1; seen.add(t)

print(f"{len(genes)} genes with per-section data\n")
for section in sorted(sec_calls):
    n=len(sec_genes[section])
    print(f"{'='*64}\n{section}  ({n} genes)\n{'='*64}")
    for t,cov in sorted(sec_calls[section].items(), key=lambda kv:(-kv[1],-sec_total[section][kv[0]])):
        core = "CORE" if cov==n else ("    " if cov>=0.6*n else "  ~ ")
        print(f"  [{core}] {cov:>2}/{n}  x{sec_total[section][t]:<4} {t}")
    print()
