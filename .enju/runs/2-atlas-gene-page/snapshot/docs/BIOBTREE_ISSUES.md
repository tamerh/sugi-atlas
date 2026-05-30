# biobtree MCP — Issues & Improvement Requests

**Date:** 2026-05-28
**Context:** Found while evaluating local LLMs (Qwen3 family) as agents driving the biobtree MCP server (`http://localhost:8000/mcp`, tools: `biobtree_search`, `biobtree_map`, `biobtree_entry`) to generate grounded gene/disease reference pages. Benchmark gene: **TP53**.

**Headline:** Most of these are *error-/empty-response signaling* problems. They are survivable for a strong model (Claude Haiku copes, at the cost of extra calls), but they actively **break weaker/local models**, which cannot tell "I queried wrong" from "the data does not exist." In one case below, a capable local model (Qwen3-30B-A3B-Instruct-2507) **looped 25 times and failed to produce any answer** purely because of issue #2.

biobtree has all the needed data — every TP53 field (incl. GRCh38 coordinates) is retrievable with the right call. The problems are about *guiding the caller to the right call and signaling failure clearly*.

---

## Issue #1 — `not_found` is overloaded (highest impact)

A query that is **malformed / resolves to an empty path** returns the **same** `not_found` as a genuinely **nonexistent entity**.

**Repro:**
```
map(terms="TP53", chain=">>hgnc>>ensembl>>entrez>>mim>>refseq")
→ { "stats": {"queried":1, "total":0}, "mappings": null, "not_found": ["TP53"] }
```
…even though each shorter hop works fine:
```
map("TP53", ">>hgnc")          → 2 mappings
map("TP53", ">>ensembl")       → 3 mappings
map("TP53", ">>hgnc>>ensembl") → 1 mapping
```

**Impact:** A model reads `not_found:[TP53]` as "TP53 has no data" and either gives up (writes "Not found") or falls back to hallucinated knowledge. This single ambiguity caused most of the field omissions/hallucinations we observed.

**Suggested fix:** Distinguish the two cases. e.g. for a recognized input whose chain yields nothing: `"status": "empty_path", "message": "Input recognized but the requested chain produced no mappings; try a shorter chain."` vs. `"status": "unknown_input"` for a truly unknown term.

---

## Issue #2 — Unknown dataset name returns a bare `API error 400` (caused a real failure)

Using a wrong/aliased dataset name returns an opaque error with no hint.

**Repro:**
```
map(terms="HGNC:11998", chain=">>hgnc>>omim")
→ { "error": "Biobtree API error: 400" }
```
The correct dataset is **`mim`**, not `omim`:
```
map(terms="HGNC:11998", chain=">>hgnc>>mim")
→ mappings: ["191170"]   ✅ (the OMIM gene ID)
```

**Impact (observed):** Qwen3-30B-A3B-Instruct-2507 — otherwise our best-behaved local model, using correct idioms and refusing to hallucinate — got stuck repeatedly calling `>>hgnc>>omim`, receiving `API error 400`, trying side-paths (mondo/gencc/clinvar), and looping. It hit the 25-step cap **without ever producing an answer**. The bare 400 gave it nothing to learn from, and it never discovered `mim`.

**Suggested fix:** Replace bare 400s with actionable text: `"unknown dataset 'omim' — did you mean 'mim'? Valid datasets: ..."`. Optionally accept `omim` as an alias for `mim`.

---

## Issue #3 — Some bad chains return an empty-but-statless blob (no error, no `not_found`)

A third distinct empty-response shape — neither an error nor a `not_found`, and missing the usual `stats`.

**Repro:**
```
map(terms="TP53", chain=">>omim")        → { } with stats=None, mappings=None, no not_found
map(terms="TP53", chain=">>hgnc>>omim")  → { } with stats=None, mappings=None, no not_found
```
(compare: `map("TP53", ">>mim")` → proper `not_found:["TP53"]`)

**Impact:** Inconsistent response schemas make it impossible for a caller (or our parser) to reliably detect "this failed." Three different "nothing came back" shapes exist (`not_found`, `error 400`, statless empty).

**Suggested fix:** One consistent response envelope with an explicit status field for every outcome (ok / empty / unknown_input / unknown_dataset / error).

---

## Issue #4 — Long multi-hop chains fail silently / intersectively

Every individual hop can succeed while the combined chain returns empty (see Issue #1 repro). There is no signal indicating *which* hop broke the chain.

**Suggested fix:** Either return partial results up to the break (with a note of where it stopped), or document prominently that **short single-hop chains are the reliable pattern** and multi-hop chains are best-effort.

---

## Issue #5 — `map` source semantics are non-obvious (doc clarity)

`map` looks the input up in the **first** dataset of the chain. Mapping *from a resolved ID* only works if the chain starts at that ID's own dataset:
```
map("HGNC:11998", ">>ensembl")        → not_found    (chain starts at ensembl; HGNC:11998 isn't an ensembl term — correct, but confusing)
map("HGNC:11998", ">>hgnc>>ensembl")  → works        (chain starts at hgnc ✅)
map("7157",       ">>ensembl")        → not_found
```
The `biobtree_map` tool description says *"HGNC:* or gene symbols → >>hgnc"*, which implies HGNC IDs are valid map inputs but does **not** make clear the chain must *start* at `>>hgnc`. Right after a model gets `HGNC:11998` from `search`, its natural next step `map(HGNC:11998, >>ensembl)` fails.

**Note:** This `not_found` is technically *correct* behavior — it's a documentation/affordance gap, not a bug.

**Suggested fix:** Clarify in the tool description that the chain must begin with the input's own dataset; ideally auto-detect the source dataset from the ID prefix (`HGNC:` → hgnc, `ENSG` → ensembl, numeric → entrez, etc.).

---

## Issue #6 — `entry` xrefs return counts, not values

`entry` shows how many links exist to each dataset, but not the linked IDs.

**Repro:**
```
entry("HGNC:11998", "hgnc")
→ xrefs: ["ensembl|1", "entrez|1", "mim|1", "refseq|93", ...]   (counts only)
```
So you know an OMIM/Ensembl/Entrez link exists, but must issue a separate `map` per dataset to get the actual value. There is no single "give me all the cross-references" call, which multiplies round-trips and the chance a weak model drops a field.

**Suggested fix:** An option to resolve xrefs to their target IDs (e.g. `entry(..., expand_xrefs=true)`), or a convenience that returns all standard gene IDs in one call.

---

## Issue #7 — `map` mixes species with no human filter

For a gene symbol, `map` returns orthologs across species with no way to restrict to the query organism.

**Repro:**
```
map("TP53", ">>ensembl")
→ ENSG00000141510 (homo_sapiens), ENSRNOG00000010756 (rat), ENSDARG00000035559 (zebrafish)
```

**Impact:** The caller must post-filter to `homo_sapiens`; a weak model may report a non-human row.

**Suggested fix:** Optional `species`/`taxon` filter, or default to the species implied by the query.

---

## Priority for us

1. **#2** (bare 400 on bad dataset name) — directly broke a capable model; cheapest, highest-impact fix.
2. **#1 / #3** (ambiguous & inconsistent empty/`not_found` signaling) — root cause of most omissions and hallucinations.
3. **#5** (map source-dataset doc clarity).
4. **#4, #6, #7** — efficiency / ergonomics.

The first three would let a small local model recover from its own mistakes instead of giving up or hallucinating — which is the difference between "usable" and "not" for local-model automation.

---

# Additional issues — deterministic REST collector (Claude, 2026-05-28)

Found while building a deterministic data collector that drives the biobtree
**REST API** (`http://127.0.0.1:8000/api`, `search`/`entry`/`map`) directly from
code (no model) to gather gene-page data. Distinct from the MCP-agent issues
above; no overlap with #1–#7.

## Issue #8 — ~~`map` pagination is broken~~ RESOLVED: wrong param name (`p` not `page`)

**RESOLVED 2026-05-28** — pagination is NOT broken. The map cursor param is **`p`**,
not `page`. Feeding the `next_token` back as `&p=<token>` advances correctly:
BRCA1 → refseq returns **759 unique IDs over 15 pages** and the MANE `NM_007294`
*is* reachable. The original report (and our collector) passed `&page=` — and
**FastAPI silently drops unknown query params**, so the cursor was ignored and
every request re-served page 1 with the same `next_token`, giving the illusion of
a stuck cursor.

```
curl "http://127.0.0.1:8000/api/map?i=ENSG00000012048&m=>>ensembl>>refseq"   # next_token=...
curl "...&p=<next_token>"   # advances ✓   (vs &page=<token> which is ignored)
```

**The real (legitimate) issue is ergonomics / silent failure** (same class as
#1–#3): the response field is `next_token` but the param is `p` — nothing tells
the caller to feed one into the other, and `page` is the natural wrong guess.
FastAPI swallowing the unknown param means a wrong guess fails *silently* instead
of erroring.

**Suggested fix:** name the request param `next_token` (match the response field)
or accept `page` as an alias; reject unknown query params with a 400 so wrong
guesses surface instead of silently re-serving page 1.

**Collector:** fixed to use `p=` — high-fan-out chains (full refseq/transcript/
interaction lists, real ClinVar per-class breakdown) are now retrieved completely.

## Issue #9 — No `reviewed` flag on UniProt (map or entry)

`>>ensembl>>uniprot` returns bare accessions (`schema: id`) mixing reviewed
(Swiss-Prot) and unreviewed (TrEMBL) with no way to tell them apart; a
`[reviewed==true]` filter returns nothing, and the uniprot `entry` exposes only
`names/alternative_names/sequence/id/name`.

**Workaround:** treat `>>hgnc>>uniprot` as the canonical reviewed set (HGNC
xrefs only curated Swiss-Prot product(s); >1 for dual-product genes like CDKN2A
→ P42771 + Q8N726).

**Suggested fix:** expose a `reviewed` boolean on uniprot map rows / entry.

## Issue #10 — `>>uniprot>>alphafold` is empty for very large proteins

For proteins >~2700 aa (ATM Q13315, BRCA2 P51587, DMD P11532) the alphafold map
returns 0 rows, although AlphaFold DB has fragmented models (AF-<acc>-F1..Fn) —
i.e. no pLDDT/id for exactly the proteins AF had to fragment.

**Workaround:** construct `AF-<acc>-F1` for every reviewed protein; attach pLDDT
only when the map provides it (F2..Fn not enumerated).

**Suggested fix:** index fragmented AlphaFold models with per-fragment metrics.

## Issue #11 — `>>uniprot>>ufeature` leaks ortholog features (distinct from #7)

Querying `ufeature` from a **unique, unambiguous UniProt accession** returns
features keyed by *other* accessions (cross-species orthologs), not just the
queried one. This is distinct from #7 (gene-symbol species mixing): a UniProt
accession like `P04637` is uniquely human — there is no ambiguity to resolve.

**Repro (TP53 / P04637):**
```
map(i="P04637", m=">>uniprot>>ufeature")   # paginated to depth, then count by prefix
→ total fetched: 2278
  prefix counts: P04637=1518 (human, correct)
                 P13481=40   (chimp Tp53)
                 P56423=40   (cotton-rat Tp53)
                 P56424=40   (woodchuck Tp53)
                 P61260=40   (cat Tp53)
                 O09185=39   (mouse Tp53)
                 Q00366=39, Q8SPZ3=39, Q9TTA1=39, Q9WUR6=39, ...
```
The result schema is `id|type|description|location_begin|location_end` where
`id` is `{source_uniprot_acc}_F{n}` — i.e. the feature ID itself encodes the
*source* accession, which can be any ortholog, not the one queried.

**Likely root cause:** the `uniprot → ufeature` edge appears to be indexed at
the **gene / ortholog-cluster level**, not per-accession — so any query within
the cluster returns the whole cluster's features.

**Impact:**
1. *Correctness risk if not filtered.* A caller that trusts the result will
   conflate ortholog features with the queried protein (e.g. a "DNA-binding
   region 102–292" from the mouse Tp53 entry rendered as a human TP53 feature).
2. *Indirect coverage loss.* Pagination is finite (we cap at 100 pages × 100
   rows = 10000). For a heavily-orthologized human protein with a very large
   ufeature set, ortholog rows can push some human features past the cap before
   they're returned — and the caller has no signal that this is happening.
3. *Same class as #7's signaling problem* — the response shape gives no hint
   that mixed-species results are being returned.

**Workaround (deterministic collector):** paginate deep, then post-filter
`id.startswith(query_accession + "_")` per reviewed UniProt product. Working,
but unnecessary if the index were per-accession.

**Suggested fix:** index `ufeature` (and any other annotation edge keyed off
UniProt) per-accession; for a unique accession query, return only that
accession's features. If the cluster-level join is intentional, add an
explicit `?species=` / `?accession=<exact>` filter to scope the result, and
document the behavior.
