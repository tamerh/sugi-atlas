#!/usr/bin/env python3
"""Benchmark summary models via OpenRouter — same deterministic body, swap models.

Builds the §1-6 rendered body for each gene (deterministic, no model), then asks
each model for a grounded executive summary. Reports latency, token usage, and a
ROUGH grounding check (IDs/large numbers in the summary not present in the body).

  export OPENROUTER_API_KEY=...        # or put the key in ~/.secrets/openrouter
  python3 biobtree/collector/bench_summary.py TP53 BRCA1 AR
  OR_MODELS="anthropic/claude-haiku-4.5,deepseek/deepseek-chat" python3 ... TP53

Summaries are saved to /tmp/sumbench/<gene>__<model>.md for quality review.
"""
import os, sys, re, json, time, urllib.request, urllib.error
from atlas.gene import collect as C
from atlas.gene import render as R

ENDPOINT = "https://openrouter.ai/api/v1/chat/completions"
MODELS = (os.environ.get("OR_MODELS") or
          "anthropic/claude-haiku-4.5,"
          "deepseek/deepseek-chat,"
          "google/gemini-2.0-flash-001,"
          "openai/gpt-4o-mini").split(",")
INSTR = ("Write a concise 3-6 sentence executive summary for the gene {g} using STRICTLY "
         "AND ONLY the facts in the body below. Rules: (1) every statement must be directly "
         "supported by the body; (2) do NOT add any biological function, role, mechanism, "
         "pathway, or significance from your own knowledge -- even if true -- unless it "
         "appears in the body; (3) do NOT infer or editorialize; (4) prefer restating "
         "concrete facts (IDs, counts, names, classifications) over interpretation. Output "
         "only the paragraph, no preamble.")

# atoms we can verify against the body: exact IDs + 3+ digit numbers (counts)
ATOM = re.compile(r"HGNC:\d+|ENS[A-Z]*\d+|NM_\d+|NP_\d+|NR_\d+|IPR\d+|PF\d{5}|CCDS\d+"
                  r"|[OPQ][0-9][A-Z0-9]{3}[0-9]|\b\d{3,}\b")

def api_key():
    k = os.environ.get("OPENROUTER_API_KEY")
    if k:
        return k.strip()
    for p in ("~/.secrets/openrouter", "~/.openrouter_key"):
        p = os.path.expanduser(p)
        if os.path.exists(p):
            return open(p).read().strip()
    sys.exit("Set OPENROUTER_API_KEY env var or put the key in ~/.secrets/openrouter")

def body_for(gene):
    return "\n\n".join(R.RENDER[s](C.SECTIONS[s](gene))
                       for s in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"])

def call(model, prompt, key, max_tokens=400, response_format=None):
    # model may be "slug" or "slug|Provider" to force a guardrail-compliant route
    provider = None
    if "|" in model:
        model, provider = model.split("|", 1)
    payload = {"model": model, "temperature": 0.0, "max_tokens": max_tokens,
               "messages": [{"role": "user", "content": prompt}]}
    if provider:
        payload["provider"] = {"only": [provider], "allow_fallbacks": False}
    if response_format is not None:
        payload["response_format"] = response_format
    data = json.dumps(payload).encode()
    req = urllib.request.Request(ENDPOINT, data=data, method="POST",
                                 headers={"Authorization": f"Bearer {key}",
                                          "Content-Type": "application/json",
                                          "X-Title": "biobtree-summary-bench"})
    t0 = time.time()
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read().decode()), time.time() - t0

def ungrounded(summary, body):
    return sorted({a for a in ATOM.findall(summary) if a not in body})

def main():
    genes = sys.argv[1:] or ["TP53"]
    key = api_key()
    os.makedirs("/tmp/sumbench", exist_ok=True)
    print(f"{'gene':<8}{'model':<36}{'p_tok':>7}{'c_tok':>7}{'lat_s':>7}{'ungr':>6}  ungrounded-atoms")
    print("-" * 100)
    for g in genes:
        body = body_for(g)
        prompt = INSTR.format(g=g) + "\n\nBODY:\n" + body
        for m in MODELS:
            m = m.strip()
            try:
                d, dt = call(m, prompt, key)
                if "choices" not in d:
                    print(f"{g:<8}{m:<36}  API: {json.dumps(d)[:70]}"); continue
                txt = str(d["choices"][0]["message"].get("content") or "").strip()
                u = d.get("usage", {})
                ug = ungrounded(txt, body)
                open(f"/tmp/sumbench/{g}__{m.replace('/', '_')}.md", "w").write(txt)
                print(f"{g:<8}{m:<36}{u.get('prompt_tokens', 0):>7}{u.get('completion_tokens', 0):>7}"
                      f"{dt:>7.1f}{len(ug):>6}  {', '.join(ug[:6])}")
            except urllib.error.HTTPError as e:
                print(f"{g:<8}{m:<36}  HTTP {e.code}: {e.read().decode('utf-8','replace')[:60]}")
            except Exception as e:
                print(f"{g:<8}{m:<36}  ERR {str(e)[:60]}")
    print("\nSummaries saved to /tmp/sumbench/*.md (read them to judge prose quality).")

if __name__ == "__main__":
    main()
