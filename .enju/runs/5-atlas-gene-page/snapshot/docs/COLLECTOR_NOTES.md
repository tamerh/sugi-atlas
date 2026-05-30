# Deterministic biobtree collector + summary layer

Generate grounded gene reference pages by **separating mechanical work from
judgment**: code gathers the data and renders the tables (deterministic, no
model); an LLM is used *only* for the executive-summary prose.

```
biobtree (local REST API)  ──►  collect.py   ──►  render.py   ──►  page tables
  127.0.0.1:8000/api            (data, no model)   (markdown, no model)
                                      │
                                      └──►  body  ──►  LLM summary  ──►  exec summary
                                                       (the ONLY model step)
```

## Why this exists

The original pipeline let an LLM drive biobtree agentically (decide calls, make
them, write prose). Mining the published pages showed that for these structured
sections the tool-call trajectory is **mostly deterministic** — 98–426 calls/page
collapse to ~3–6 canonical chains per section; the rest is agentic waste (syntax
drift, alt-paths, redundant fan-out). So the data-gathering and table-formatting
need **no model at all** — which is cheaper, reproducible, and *more accurate*.

## Files

| File | Role |
|---|---|
| `collect.py` | 12 section collectors → structured JSON bundle (no model). `python3 collect.py TP53 6` |
| `render.py` | bundle → markdown tables, HTML-unescaped (no model). `python3 render.py TP53 all` |
| `coverage.py` | validates a collector vs the historical per-section output across the 25 committed genes. `python3 coverage.py 2` |
| `bench_summary.py` | compares summary models via OpenRouter on the full body (latency/tokens/atom-check). `OPENROUTER_API_KEY=… python3 bench_summary.py TP53 BRCA1 AR` |
| `judge.py` | LLM faithfulness gate — flags summary claims not supported by the body. `python3 judge.py TP53` |
| `../../BIOBTREE_ISSUES.md` | biobtree API issues found (shared with the local-LLM effort) |

## Sections (gene)

1 gene_ids · 2 transcripts · 3 protein_ids · 4 structure · 5 orthologs ·
6 variants · 7 pathways · 8 interactions · 9 tf_regulation · 10 drugs ·
11 expression · 12 diseases — all built and rendering.

## Key findings

### Collector is *more correct* than the agentic model
Validated §1–6 against all 25 committed genes (FULL fact coverage). Where the
deterministic collector differs from the old pages, it's because it's right:
- picks the true **MANE-Select** transcript (model chose ad-hoc ones)
- stays **human-only** (model pulled rat/mouse orthologs)
- unions domains across **dual-product** genes (CDKN2A p16+p14ARF)
- direction-filtered **CollecTRI** so high-degree TFs (TP53) get downstream targets

The 25-gene coverage check caught a real bug each section (ID-resolution
pagination/ambiguous symbols; canonical-transcript definition; dual-product
domains; AlphaFold id format).

### Reliable anchors / techniques
- `>>hgnc>>uniprot` = the canonical reviewed Swiss-Prot product(s).
- MANE transcript via `>>ensembl>>refseq[is_mane_select==true]` → `>>refseq>>transcript`.
- The **hgnc `entry` xref table** gives exact totals (clinvar/spliceai/msigdb/hpo/
  gwas/collectri…) in one call — used instead of counting (which the 100-cap breaks).
- **Chain FILTERS** beat the broken pagination (see BIOBTREE_ISSUES.md #8).

### Summary model comparison (full 12-section body, TP53/BRCA1/AR)
Cost is negligible at any scale (≈$12 for all 20k human genes even on Haiku), so
the choice is **faithfulness + latency**, not price.

| Model (route) | latency | grounding | note |
|---|---|---|---|
| **Qwen3-235B-A22B** (Together) | 3–7s | excellent | comprehensive + grounded; ~16× cheaper than Haiku — **the pick** |
| Claude Haiku 4.5 (Anthropic) | 3–6s | excellent | quality anchor |
| DeepSeek-v3.1 (DeepInfra) | 1.3–3s | tightest | terse but fastest/cheapest |
| Gemini 2.0 Flash (Google) | 2–9s | good | one intermittent truncation observed |
| GPT-4o-mini (OpenAI) | 3–7s | weakest | invents an occasional count (e.g. AR "1,053 variants") |

DeepSeek V4 (Pro/Flash) are **reasoning models** — overkill/slow for a summary;
V4-Flash returns empty unless given a large reasoning budget.

### The prompt matters more than the model
With a loose prompt ("highlight what's notable"), even strong models injected
**external knowledge** ("tumor suppressor… regulates the cell cycle… DNA damage
response") — true but **not in the body**. A **strict prompt** (forbid any
function/role/mechanism claim not in the body; restate concrete facts only)
**eliminated the drift** while keeping coverage. That strict prompt is the
production prompt (see `INSTR` in `bench_summary.py`).

### Faithfulness checking
- **Atom-check** (IDs/large numbers in summary present in body): reliable for
  IDs, blind to prose drift — scored everything 0 even when models editorialized.
- **LLM-judge** (`judge.py`): catches prose drift the atom-check misses, BUT a
  Haiku judge has a **high false-positive rate** — it flags body-derived
  restatements (aliases, transcript counts, PDB methods) as "unsupported". Use it
  as a **screen** (read the flagged claims), not a verdict. TODO: neutral/stronger
  judge, robust JSON, feed it the body's term list to cut false positives.

## Dataset-wiring backlog (from `dataset_coverage.py`)

Reference = the MCP **map edges** (`/api/help?topic=edges`), the real
what-connects-to-what graph — **117 queryable datasets**, collector uses **41**.
(Do NOT use `/api/meta`'s 500+ list — most are derived xref-id types, e.g. a
"cosmic" id that is NOT integrated somatic data.) Run `dataset_coverage.py` to
see coverage and to spot newly added datasets. Status:
- ✅ **§7** GO + Reactome unioned across **every reviewed uniprot** + the
  gene-level ensembl route. UniProt-GOA vs Ensembl annotation diverge ~20%, and
  dual-product genes carry distinct terms per product (CDKN2A p14ARF mitophagy GO
  absent from p16/ensembl; the canonical `cu[0]` alone missed them). §7 fact-cov
  24/25 FULL; the lone residual (AKT1 R-HSA-8949469) is in neither reactome route
  — a biobtree reactome-edge gap, not collector-fixable.
- ✅ **§11** `bgee_evidence` — per-tissue expression scores (fixed the gene-level gap)
- ✅ **§11** `fantom5_gene` — FANTOM5 CAGE gene-level expression (tpm_avg/max,
  samples_expressed, breadth). The model's other "expression" sources are
  **xref-only, NOT integrated** in biobtree (not edge nodes; `entry` returns
  `Attributes:Empty`): `hpa`, `gtex`, `expressionatlas`, `proteomicsdb`,
  `tabula_sapiens`. Wiring them would just emit bare ids — not done (no data to
  cover). §11 fact-coverage is FULL (0 missing IDs across 25 genes).
  `scxa_expression` is a real edge but the graph is experiment-centric, so the
  route returns an all-genes list, not gene-scoped per-cell expression.
- ✅ **§5** `paralog` (e.g. TP53 → TP63, TP73)
- ✅ **§8** `signor` — directed signaling (effect + mechanism)
- ✅ **§10** `bindingdb` affinities; molecules now filtered to phased DRUGS
  (`highestDevelopmentPhase>=1`) and ranked by phase (was: thousands of
  ID-ordered screening compounds)
- ✅ **§6** `dbsnp` — via `>>hgnc>>entrez>>dbsnp` (sampled rsIDs). Direct
  `>>hgnc>>dbsnp` is empty: dbsnp records xref entrez not hgnc, yet the edge
  guide declares hgnc→dbsnp (reported to biobtree dev).
- ⏳ **§10** `gtopdb` — PENDING a biobtree fix (do NOT skip). Data exists, but the
  `map` lite-projection emits empty fields (`>>uniprot>>gtopdb` → `1797|||`); the
  full format exists at the backend (`mode=full`) but isn't exposed via REST.
  Reported to biobtree dev; fix expected ~2026-05-30. **Revisit when it lands**:
  wire the GtoPdb target classification (type/family) into §10, and use the
  interaction record for affinities.
- ✅ **§10** `clinical_trials` — via the DISEASE route `>>hgnc>>{gencc,clinvar}>>
  mondo>>clinical_trials` (biobtree's intended pattern; clinical_trials xrefs
  MONDO/ChEMBL and there's no gene→drug→trials path). Disease-level, not
  drug-specific, but clean + relevant (BRCA1 → breast-cancer/HBOC trials). The
  `chembl_molecule>>clinical_trials` route is avoided — ChEMBL target→molecule is
  bioactivity-based, so off-target drugs (Levodopa@phase4 vs EGFR) pollute it.
  No clean small-molecule mechanism-of-action edge exists in biobtree; gene→drugs
  uses bioactivity + phase filter (biobtree's own JAK2 test pattern).
- TODO: `ctd_disease_association` (§12), `scxa_expression` (§11),
  `pharmgkb_guideline` (§10 dosing), `cellphonedb`, `fantom5_*`, `ufeature`.

**Somatic (in progress, wait):** CIViC + intOGen + one more are being
implemented in biobtree. They are NOT yet edges — they'll appear as new nodes in
`dataset_coverage.py` once integrated; wire a somatic block then.

## Open refinements
- ~~map pagination~~ RESOLVED — the cursor param is `p` not `page` (FastAPI
  silently dropped the wrong param). Collector now paginates fully; interaction
  counts, CollecTRI breadth, and ClinVar per-class breakdown are now real (e.g.
  TP53 ClinVar: 749 Pathogenic / 954 VUS / …). See BIOBTREE_ISSUES.md #8.
- ~~Interaction per-edge scores~~ RESOLVED — query the interaction *record*
  (`>>uniprot>>string_interaction` etc.), not `>>...>>uniprot` (which collapses
  to bare partner ids). STRING score (×1000), IntAct confidence, BioGRID method.
- **Per-tissue** Bgee expression not collected (the map is gene-level).
- §7–12 not yet run through the 25-gene coverage harness (§1–6 are).
- Harden `judge.py`; consider an automated grounding gate in the pipeline.

## Parallel local-LLM effort
A separate investigation (`/data/localllm`) independently reached the **same
hybrid architecture** (deterministic lookups + local model as formatter-only +
programmatic grounding validation); its leading local model is
Qwen3-30B-A3B-Instruct-2507.
