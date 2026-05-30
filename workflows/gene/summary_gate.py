#!/usr/bin/env python3
"""Hardened two-pass LLM judge over the executive summary."""
import json, os
from atlas.validation import summary_gate

ctx = json.load(open(os.path.join(os.environ["ENJU_RUN_DIR"], "context.json")))
symbol = ctx["iteration"]["symbol"]
skip = str(ctx["params"].get("skip_summary", "false")).lower() == "true"

if skip:
    print(f"summary_gate {symbol}: SKIPPED")
    open(f"build/{symbol}/judge.json", "w").write('{"verdict":"skipped"}')
    raise SystemExit(0)

body = open(f"build/{symbol}/body.md").read()
summary = open(f"build/{symbol}/summary.md").read()
r = summary_gate.check_summary(body, summary)
with open(f"build/{symbol}/judge.json", "w") as f:
    json.dump(r, f, indent=2)
print(f"summary_gate {symbol}: {r['verdict']} both={len(r.get('both', []))} single={len(r.get('single', []))}")
