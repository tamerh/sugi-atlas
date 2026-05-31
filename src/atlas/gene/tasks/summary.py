#!/usr/bin/env python3
"""LLM executive summary via OpenRouter. Skips on params.skip_summary=true."""
import json, os
from atlas.bench import summary as B

ctx = json.load(open(os.path.join(os.environ["ENJU_RUN_DIR"], "context.json")))
symbol = ctx["iteration"]["symbol"]
skip = str(ctx["params"].get("skip_summary", "false")).lower() == "true"
model = ctx["params"].get("summary_model") or "qwen/qwen3-235b-a22b-2507|Together"

if skip:
    print(f"summary {symbol}: SKIPPED")
    open(f"build/{symbol}/summary.md", "w").write("")
    raise SystemExit(0)

body = open(f"build/{symbol}/body.md").read()
key = B.api_key()
prompt = B.instruction(symbol, "gene") + "\n\nBODY:\n" + body
d, dt = B.call(model, prompt, key, max_tokens=600)
if "choices" not in d:
    raise SystemExit("API: " + str(d)[:200])
txt = str(d["choices"][0]["message"].get("content") or "").strip()
open(f"build/{symbol}/summary.md", "w").write(txt + "\n")
print(f"summary {symbol}: {len(txt)}c in {dt:.1f}s ({model})")
