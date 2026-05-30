#!/usr/bin/env python3
"""V2: collapse calls to their TEMPLATE VOCABULARY to size the deterministic core.

We don't care how many times the model drilled into individual IDs (that's
mechanical fan-out code can do in a loop). We care about the distinct *chains*
and *entry/search sources* that produce the data. So: ignore inputs entirely,
key only on (call-type, chain-or-source). Report how concentrated call volume
is in the top templates -> tells us how much is templatable vs long-tail.
"""
import os, re, glob, collections, statistics
import yaml

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = os.path.join(REPO, "hugo/content/biobtree")

def load_fm(path):
    t = open(path).read()
    if not t.startswith("---"): return {}
    e = t.find("\n---", 3)
    if e < 0: return {}
    try: return yaml.safe_load(t[3:e]) or {}
    except Exception: return {}

def template(desc):
    """Reduce a desc to its call-type + chain/source, dropping all inputs."""
    desc = desc.strip()
    m = re.match(r"^(\w+)\((.*)\)$", desc, re.S)
    if not m: return ("?", desc[:40])
    typ, args = m.group(1), m.group(2)
    if typ == "map":
        # chain is the part starting at the first '>>'
        i = args.find(">>")
        chain = args[i:] if i >= 0 else args
        # normalize numeric thresholds inside filters to keep them (they're fixed)
        return ("map", chain)
    # entry/search: signature is the source (last bare arg) if present
    parts = [p.strip() for p in args.split(",")]
    src = parts[-1] if len(parts) > 1 and re.match(r"^[a-z_]+$", parts[-1]) else "(default)"
    return (typ, src)

def analyze(cat):
    files = sorted(glob.glob(f"{ROOT}/{cat}/*.md"))
    if not files: return
    n = 0
    cov = collections.Counter(); tot = collections.Counter()
    per_page = []
    for f in files:
        fm = load_fm(f); calls = fm.get("api_calls") or []
        if not isinstance(calls, list): continue
        n += 1; per_page.append(len(calls))
        seen = set()
        for c in calls:
            d = (c or {}).get("desc","") if isinstance(c,dict) else ""
            if not d: continue
            t = template(d); tot[t]+=1
            if t not in seen: cov[t]+=1; seen.add(t)
    total_calls = sum(tot.values())
    print(f"\n{'='*70}\n{cat}: {n} pages, {total_calls} total calls, "
          f"{len(tot)} distinct TEMPLATES (vs raw signatures)\n{'='*70}")
    print(f"calls/page min={min(per_page)} median={int(statistics.median(per_page))} max={max(per_page)}")
    # concentration: cumulative % of call volume by top templates
    ranked = sorted(tot.items(), key=lambda kv:-kv[1])
    cum=0
    print(f"\n{'cov':>4} {'tot':>5} {'cum%':>5}  template")
    print("-"*70)
    for t,c in ranked:
        cum += c
        core = "*" if cov[t]==n else (" " if cov[t]>=0.5*n else ".")
        print(f"{cov[t]:>4} {c:>5} {100*cum/total_calls:>4.0f}%  {core} {t[0]:5} {t[1][:50]}")
    print(f"\nlegend: * = on every page (deterministic core), "
          f"blank = >=50% pages, . = sparse/long-tail")

for cat in ("gene","disease"):
    analyze(cat)
