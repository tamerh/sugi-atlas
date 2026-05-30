#!/usr/bin/env python3
"""Body-gate regression check. Exits non-zero on regression so Enju marks the task failed_retryable."""
import json, os, sys
from atlas.validation import body_gate

ctx = json.load(open(os.path.join(os.environ["ENJU_RUN_DIR"], "context.json")))
symbol = ctx["iteration"]["symbol"]

bundle = json.load(open(f"build/{symbol}/bundle.json"))
r = body_gate.check(symbol, bundle)
with open(f"build/{symbol}/body_gate.json", "w") as f:
    json.dump(r, f, indent=2)
print(f"body_gate {symbol}: {r['verdict']} ({r['summary']})")
sys.exit(0 if r["verdict"] in {"clean", "drift", "first_run"} else 2)
