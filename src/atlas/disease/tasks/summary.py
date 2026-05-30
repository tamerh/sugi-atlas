#!/usr/bin/env python3
"""LLM executive summary via OpenRouter. Skips on params.skip_summary=true.
Reuses atlas.bench.summary like the gene side; the prompt label is the
canonical disease name."""
import json, os
from atlas.bench import summary as B
from atlas.disease.slug import slugify

ctx = json.load(open(os.path.join(os.environ["ENJU_RUN_DIR"], "context.json")))
disease = ctx["iteration"]["disease"]
name = disease["name"] if isinstance(disease, dict) else str(disease)
slug = (disease.get("slug") if isinstance(disease, dict) else None) or slugify(name)
skip = str(ctx["params"].get("skip_summary", "false")).lower() == "true"
model = ctx["params"].get("summary_model") or "qwen/qwen3-235b-a22b-2507|Together"

if skip:
    print(f"summary {slug}: SKIPPED")
    open(f"build/{slug}/summary.md", "w").write("")
    raise SystemExit(0)

# Prefer the resolved canonical name from collect_render for the prompt label.
am_path = f"build/{slug}/anchors_meta.json"
if os.path.exists(am_path):
    am = json.load(open(am_path))
    label = am.get("canonical_name") or name
else:
    label = name

body = open(f"build/{slug}/body.md").read()
key = B.api_key()
prompt = B.INSTR.format(g=label) + "\n\nBODY:\n" + body
d, dt = B.call(model, prompt, key, max_tokens=600)
if "choices" not in d:
    raise SystemExit("API: " + str(d)[:200])
txt = str(d["choices"][0]["message"].get("content") or "").strip()
open(f"build/{slug}/summary.md", "w").write(txt + "\n")
print(f"summary {slug}: {len(txt)}c in {dt:.1f}s ({model})")
