#!/usr/bin/env python3
"""LLM-judge faithfulness check over bench summaries in /tmp/sumbench.

For each <gene>__<model>.md, regenerate the deterministic body and ask a judge
model which statements in the SUMMARY are NOT supported by the BODY. The judge
scores against the BODY only — a statement true in general but absent from the
BODY counts as UNSUPPORTED. This catches the prose drift the atom-check misses
(e.g. "tumor suppressor / regulates the cell cycle" added from training).

Four hardenings vs. the original judge:
  1. **Atom-list hint** — every id/accession/symbol/large-number is pre-extracted
     from the body and given to the judge as "GROUNDED-ATOMS". Kills the big
     false-positive class where the judge flagged body-derived ids/counts.
  2. **JSON mode** — OpenRouter `response_format=json_object` so the model returns
     pure JSON. Regex fallback retained for providers that ignore it.
  3. **Neutral judge default** — DeepSeek-v3.1 via DeepInfra (different family
     from the summary models so Haiku doesn't judge Haiku).
  4. **Two-pass adversarial verify** — set OR_JUDGE2 to a second model and only
     claims flagged by BOTH judges are reported as high-confidence unsupported;
     claims flagged by one are surfaced as "single-vote" (for inspection).

  export OPENROUTER_API_KEY=...                # or ~/.secrets/openrouter
  OR_JUDGE=deepseek/deepseek-chat|DeepInfra    # primary judge (default)
  OR_JUDGE2=google/gemini-2.0-flash-001        # optional second judge
  python3 judge.py [gene ...]                  # default: all summaries in /tmp/sumbench
  python3 judge.py --dry TP53                  # offline structural check (no API)
"""
import os, sys, re, json, glob
from atlas.bench import summary as B

JUDGE  = os.environ.get("OR_JUDGE",  "claude:claude-sonnet-4-6")
JUDGE2 = os.environ.get("OR_JUDGE2", "")  # empty = single-judge mode

# Two backends:
#   "claude:<model>"   → invokes the local `claude -p` CLI (uses the operator's
#                        Anthropic subscription; strong reasoning, well-suited to
#                        strict-grounding. Empirically: Sonnet returns {"unsupported": []}
#                        on the same prompt where DeepSeek-v3.1 flags 11/11 false positives).
#   "<vendor/slug>[|Provider]" → OpenRouter via B.call (multi-vendor, cost-tunable for scale).

# Atoms grounded in the body — high recall (better to over-include than to let
# a body-derived id get flagged as a hallucination).
ATOM_RE = re.compile(r"""
    HGNC:\d+ | ENS[A-Z]{0,4}\d+ | NM_\d+ | NR_\d+ | NP_\d+ | CCDS\d+ |
    [OPQ][0-9][A-Z0-9]{3}[0-9] | P\d{5} |
    AF-[A-Z0-9]+-F\d+ |
    IPR\d+ | PF\d{5} |
    R-HSA-\d+ | GO:\d+ | M\d{4,} |
    CHEMBL\d+ | NCT\d+ | MA\d+\.\d+ |
    MIM:\d+ | MONDO:\d+ | HP:\d+ | Orphanet:\d+ | GCST\d+ |
    rs\d+ | SIGNOR-\d+ | EBI-\d+ |
    \d{2,}\.\d+ | \b\d{3,}\b
""", re.X)
SYMBOL_RE = re.compile(r"\b[A-Z][A-Z0-9]{1,7}\b")

def body_atoms(body):
    return set(ATOM_RE.findall(body)) | set(SYMBOL_RE.findall(body))

PROMPT = """You are a strict faithfulness checker. You are given a BODY of factual data, a SUMMARY written from it, and a GROUNDED-ATOMS list pre-extracted from the body (every id, accession, number, gene symbol, and dataset term that appears in the body).

List every statement in the SUMMARY that is NOT supported by the BODY. Any fact, number, name, or claim (function, mechanism, role, biological role) that does not appear in or directly follow from the BODY counts as UNSUPPORTED.

DO NOT flag a claim simply because exact words differ. If its identifiers, numbers, gene symbols, or dataset terms are in GROUNDED-ATOMS and the claim restates body data, accept it. Only flag when the SUMMARY adds factual content the BODY does not contain. A statement true in general but absent from the BODY is UNSUPPORTED.

Return STRICT JSON only, no prose: {"unsupported": ["<short verbatim offending phrase>", ...]}. Empty list if fully grounded.

GROUNDED-ATOMS (sample, body is authoritative):
%s

BODY:
%s

SUMMARY:
%s"""

def normalize(s):
    return re.sub(r"\W+", " ", (s or "").lower()).strip()[:120]

def _parse_unsupported(txt):
    """Extract {"unsupported": [...]}. Try direct JSON first, then regex-extract
    a {...} blob, then return (None, 'parse-fail')."""
    try:
        return json.loads(txt).get("unsupported", []), None
    except Exception:
        pass
    m = re.search(r"\{.*\}", txt, re.S)
    if not m:
        return None, "parse-fail"
    try:
        return json.loads(m.group(0)).get("unsupported", []), None
    except Exception:
        return None, "parse-fail"

def judge_once(body, atoms, summary, key, model):
    import subprocess, time
    sample = ", ".join(sorted(atoms)[:300])
    prompt = PROMPT % (sample, body, summary)
    t0 = time.time()

    if model.startswith("claude:"):
        # Local `claude -p` CLI — operator's Anthropic subscription, strong judge.
        # `key` is unused in this branch (Claude auth is via the CLI's own config).
        clm = model.split(":", 1)[1]
        try:
            res = subprocess.run(
                ["claude", "-p", "--model", clm],
                input=prompt, capture_output=True, text=True, timeout=180,
            )
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            return None, time.time() - t0, f"claude-cli: {e}"
        if res.returncode != 0:
            return None, time.time() - t0, f"claude-cli rc={res.returncode}: {res.stderr[:200]}"
        uns, err = _parse_unsupported(res.stdout)
        return uns, time.time() - t0, err

    # OpenRouter path (kept for cost-tunable scale + multi-judge diversity).
    d, dt = B.call(model, prompt, key, max_tokens=800,
                   response_format={"type": "json_object"})
    if "choices" not in d:
        return None, dt, f"API: {json.dumps(d)[:80]}"
    txt = str(d["choices"][0]["message"].get("content") or "")
    uns, err = _parse_unsupported(txt)
    return uns, dt, err

def check_summary(body, summary, key=None):
    """In-process API used by pipeline.py. Runs the configured judges (one or
    two via OR_JUDGE / OR_JUDGE2), returns {verdict, both, single, raw}.

    verdict:
      clean       — both judges (or the only judge) returned []
      flagged     — at least one high-confidence (both-vote) unsupported claim
      single_vote — only single-vote claims, no high-confidence ones
      error       — judge parse/API failure
    """
    if key is None:
        key = B.api_key()
    atoms = body_atoms(body)
    judges = [JUDGE] + ([JUDGE2] if JUDGE2 and JUDGE2 != JUDGE else [])
    flagsets, err_seen = [], None
    for j in judges:
        uns, dt, err = judge_once(body, atoms, summary, key, j)
        if uns is None:
            err_seen = err; break
        flagsets.append({normalize(u): u for u in uns})
    if err_seen:
        return {"verdict": "error", "error": err_seen, "judges": judges}
    if len(flagsets) == 1:
        both, single = list(flagsets[0].values()), []
    else:
        ks = [set(fs.keys()) for fs in flagsets]
        common = ks[0] & ks[1]
        both = [flagsets[0][k] for k in common]
        single = sorted({v for fs in flagsets for k, v in fs.items() if k not in common})
    verdict = "flagged" if both else ("single_vote" if single else "clean")
    return {"verdict": verdict, "both": both, "single": single, "judges": judges}

def main():
    args = [a for a in sys.argv[1:] if a != "--dry"]
    dry = "--dry" in sys.argv
    only = set(args)
    files = sorted(glob.glob("/tmp/sumbench/*.md"))
    if not files:
        print("no summaries in /tmp/sumbench — run bench_summary.py first"); return
    judges = [JUDGE] + ([JUDGE2] if JUDGE2 and JUDGE2 != JUDGE else [])
    print(f"primary judge: {JUDGE}")
    if len(judges) > 1:
        print(f"second judge:  {JUDGE2}  (two-pass intersect mode)")
    print()

    if dry:
        for g in sorted({os.path.basename(f).split('__')[0] for f in files}):
            if only and g not in only: continue
            body = B.body_for(g); atoms = body_atoms(body)
            print(f"  {g}: body={len(body)}c atoms={len(atoms)} sample={sorted(atoms)[:5]}")
        return

    key = B.api_key()
    header = f"{'gene':<7}{'model':<42}{'jdg':>4}{'both':>5}{'single':>7}  high-confidence unsupported"
    print(header); print("-" * (len(header) + 30))
    for f in files:
        base = os.path.basename(f)[:-3]
        gene, _, model = base.partition("__")
        if only and gene not in only: continue
        summary = open(f).read().strip()
        if not summary:
            print(f"{gene:<7}{model:<42}  EMPTY"); continue
        body = B.body_for(gene); atoms = body_atoms(body)
        flagsets, err_seen = [], None
        for j in judges:
            uns, dt, err = judge_once(body, atoms, summary, key, j)
            if uns is None:
                err_seen = err; break
            flagsets.append({normalize(u): u for u in uns})
        if err_seen:
            print(f"{gene:<7}{model:<42}  ?  ({err_seen})"); continue
        if len(flagsets) == 1:
            both = list(flagsets[0].values()); single = []
        else:
            keys = [set(fs.keys()) for fs in flagsets]
            common = keys[0] & keys[1]
            both = [flagsets[0][k] for k in common]
            single = sorted({v for fs in flagsets for k, v in fs.items() if k not in common})
        msg = " | ".join(both[:4])
        if single: msg += f"  [single-vote: {' | '.join(single[:2])}]"
        print(f"{gene:<7}{model:<42}{len(judges):>4}{len(both):>5}{len(single):>7}  {msg}")

if __name__ == "__main__":
    main()
