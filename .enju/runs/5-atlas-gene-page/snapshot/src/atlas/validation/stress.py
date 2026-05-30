#!/usr/bin/env python3
"""Stress the collector across a structurally diverse gene set NOT in the 25
reference set. We have no historical pages for these — so we can't run
fact-coverage. Instead surface ANOMALIES we can act on:

  - resolution / collector EXCEPTIONS
  - sections returning empty (no data at all)
  - structural surprises in counts (e.g. a "protein-coding" gene with 0 proteins
    OR an ncRNA with proteins, etc.)

Usage: python3 stress.py
"""
import json, sys, os, time, traceback
from atlas.gene import collect as C

# Diverse stress set covering categories under-represented in the 25:
#   lncRNA, mtDNA-encoded, HLA polymorphic, Y-linked, mega-protein, pharmacogenes
GENES = [
    ("MALAT1",   "lncRNA"),
    ("XIST",     "lncRNA / female-specific"),
    ("TTN",      "largest protein"),
    ("HLA-A",    "HLA polymorphic"),
    ("MT-ND1",   "mitochondrial-encoded"),
    ("CYP2C19",  "pharmacogene"),
    ("VKORC1",   "pharmacogene / tiny"),
    ("KRAS",     "small oncogene"),
    ("NF1",      "large tumor suppressor"),
    ("SRY",      "Y-linked"),
]

SECTIONS = [
    ("1",  C.collect_gene_ids),
    ("2",  C.collect_transcripts),
    ("3",  C.collect_protein_ids),
    ("4",  C.collect_structure),
    ("5",  C.collect_orthologs),
    ("6",  C.collect_variants),
    ("7",  C.collect_pathways),
    ("8",  C.collect_interactions),
    ("9",  C.collect_tf_regulation),
    ("10", C.collect_drugs),
    ("11", C.collect_expression),
    ("12", C.collect_diseases),
]

def size_of(v):
    if isinstance(v, list): return len(v)
    if isinstance(v, dict): return len(v)
    if v is None or v == "": return 0
    return 1

def section_signal(b):
    """Compact key=size summary, skipping meta + empty/zero keys."""
    skip = {"section", "symbol"}
    items = []
    for k, v in b.items():
        if k in skip: continue
        s = size_of(v)
        if s == 0: continue
        items.append(f"{k}={s}")
    return " ".join(items) if items else "(all-empty)"

def main():
    out = []
    for sym, note in GENES:
        print(f"\n=== {sym} ({note}) ===")
        per = {}
        for sec, fn in SECTIONS:
            t0 = time.time()
            try:
                b = fn(sym)
                sig = section_signal(b)
                dt = time.time() - t0
                per[sec] = (sig, dt, None)
                print(f"  §{sec:>2} ({dt:4.1f}s) {sig[:140]}")
            except SystemExit as e:
                per[sec] = (None, time.time()-t0, f"SystemExit: {e}")
                print(f"  §{sec:>2} SystemExit: {e}")
            except Exception as e:
                per[sec] = (None, time.time()-t0, f"{type(e).__name__}: {e}")
                print(f"  §{sec:>2} ERROR {type(e).__name__}: {e}")
        out.append((sym, note, per))
    print("\n" + "=" * 80 + "\nANOMALIES (errors or sections returning ALL-EMPTY)\n")
    for sym, note, per in out:
        flags = []
        for sec, fn in SECTIONS:
            sig, dt, err = per[sec]
            if err: flags.append(f"§{sec}={err[:40]}")
            elif sig == "(all-empty)": flags.append(f"§{sec}_EMPTY")
        print(f"  {sym:<8} ({note:<26}) {' | '.join(flags) or '(clean)'}")

if __name__ == "__main__":
    main()
