# How Atlas works

Atlas turns a local [biobtree](https://biobtree.org) instance into deterministic
reference pages for genes, drugs, and diseases. This doc covers the machinery
shared by all three entity types; the per-entity mining is in
[gene.md](gene.md), [drug.md](drug.md), and [disease.md](disease.md), and the
cross-entity graph in [mesh.md](mesh.md).

## The core principle: code gathers, the model only summarizes

Every fact on an Atlas page is gathered and rendered by code — no model decides
what to query or how to format a table. The pipeline is:

```
biobtree ──► collect ──► render ──► body ──► (LLM summary) ──► publish
            (no model)  (no model)          (optional, gated)
```

The deterministic body, the one-line lead sentence, the JSON-LD, and the
cross-entity mesh are all pure functions of the biobtree bundle. The lead
sentence in particular is assembled from collected fields — *not* written by a
model — so it has no hallucination surface. Even the statistics are deterministic
code, not heuristics: disease cohort enrichment (Reactome / GO / protein-family)
and gene interaction-partner enrichment are real over-representation analysis — a
hypergeometric test against a precomputed genome-wide background with
Benjamini-Hochberg FDR (`atlas.ora`, no scipy). The **only** model touch is an
optional executive summary, which is **off by default**, fenced by an
adversarial faithfulness gate when on, and always disclosed as
*"Summary written by {model} from the deterministic data below."*

This makes the whole corpus reproducible: a `generated_at` stamp anchors each
page, and every biobtree query is recorded in a per-call log.

## biobtree transport

All data flows through one client (`biobtree/client.py`) with three primitives —
`search`, `entry`, `bbmap` (the mapping/chain call) — plus pure parsers and a
reproducibility log. There are two interchangeable transports with
byte-identical responses:

| Transport | Path | Use |
|---|---|---|
| **direct** (default) | biobtree Go server `/ws/`, `mode=lite` injected | ~10× faster for bulk page-gen — no uvicorn hop, no double-JSON, no N+1 not-found probes |
| **gate** | FastAPI/MCP `/api/` | fallback for environments that only expose the gate |

Selected with `ATLAS_BIOBTREE_TRANSPORT`; base URL with `ATLAS_BIOBTREE`.

**Chain syntax.** The `bbmap` call walks the biobtree edge graph with
`>>dataset>>dataset…` hops and optional inline filters:
`>>hgnc>>uniprot`, `>>ensembl>>refseq[is_mane_select==true]`,
`>>chembl_target>>chembl_molecule[highestDevelopmentPhase>=1]`. `map_all`
paginates a chain (cursor param `p`) and dedupes, loop-safe with a page cap.

Three reliability details are load-bearing and reused by every section:

- **Universal HTML-unescape.** Every response string is recursively
  `html.unescape`'d exactly once, at the single point all responses flow
  through — so every consumer renders clean text (`3'-end`, not `3&apos;-end`).
  Safe because entities never contain a raw `|`.
- **Column-shift guards.** Responses are pipe-delimited against a `schema`
  header. Any row whose field count ≠ the column count is dropped, because an
  unescaped pipe shifting the columns would silently misalign every later value
  (the "MeSH-row bug class"). The map parser additionally protects escaped `\|`
  before splitting.
- **Retry + terminal errors.** Bounded exponential backoff on 5xx/429/timeout/
  truncated bodies; a 4xx, or biobtree's inline `{"Err": …}` body (the direct Go
  transport returns query errors with HTTP 200), raises a typed `BiobtreeError`
  *terminally* — a bad chain won't fix itself on retry, and the batch driver
  skips+logs that one entity rather than crashing.

## Resolution: seed by ID, not name

Each entity resolves *once* to an immutable anchor record that all its sections
share (eliminating dozens of redundant lookups). The corpus is seeded by **stable
IDs, not names**, across all three types — HGNC symbol, MONDO id, ChEMBL id —
because biobtree's name search ranks by xref count, so a name resolves to the
highest-xref *subtype* rather than the umbrella term (`cardiomyopathy` →
*dilated cardiomyopathy 1G*) or a nameless screening entry outranks the real drug
(`Vemurafenib` → an empty-named ChEMBL row). Seeding by id pins the resolution
decision once; slugs still derive from the canonical name.

## The frozen page contract

Every page of a given type emits the **same canonical sequence of `## H2`
zones**, in the same order, each with a stable lower-kebab `{#anchor}` id — even
when a zone is empty (it gets an *informative placeholder*). The contract is
specified in [PAGE_CONTRACT.md](PAGE_CONTRACT.md) and enforced over the built
corpus by the integration tests. Three helpers implement it:

- `emit_canonical` lays out the frozen H2 sequence, body or placeholder.
- `with_heading_id` stamps an explicit `{#anchor}` on each sub-section's
  heading — a **backend-owned** id, never Hugo's prose-derived auto-id (which
  would break a deep link the moment a heading count changed).
- `demote` nests each renderer's section as `### H3` under its canonical H2; the
  data tables inside a section are `#### H4` with their own stable anchors (the
  `H4_IDS` contract), so the TOC nests H2 → H3 → H4 and any table is deep-linkable.

The anchor ids are treated as a **stable public API**: external links and the
JSON-LD `@id`s depend on them, so they don't drift.

## Frontmatter and key facts

A single builder produces the Hugo frontmatter for both the per-page and batch
paths, so the shape can't diverge. Beyond title/slug it carries machine-facing
fields derived deterministically from the bundle: a typed `identifier` (HGNC
symbol / MONDO id / ChEMBL id — templates key on this, not the display symbol),
`alt_names` for search (deliberately *not* Hugo's reserved `aliases:`, which
would emit 301s), a `tldr` digest, and `section_defaults` (open/collapsed hints
keyed by anchor id). The published `datasets:` list is the *true* per-page source
set, read from the call log — including cohort fan-out. The raw bundle and the
API-call trail are intentionally not published.

## schema.org JSON-LD

Each page emits a typed schema.org node (`Gene` / `MedicalCondition` / `Drug`) in
two forms: an inline `<script type="application/ld+json">` in the body, and a
machine-fetchable `entity.jsonld` sidecar. The highest-leverage field is
**`sameAs`** — the federated-identity cross-refs (HGNC/NCBI/UniProt/Ensembl/OMIM;
Mondo/EFO/MeSH/Orphanet; ChEMBL/PubChem/ChEBI/ATC) that let an agent decide "is
this the same entity?" — which none of the surveyed competitors emit.

Gene nodes go further: one fully-addressable `Protein` node per reviewed UniProt
product (`@id …/gene/SYM/#protein-<acc>`, matching the page anchor), with the
reciprocal `encodesBioChemEntity`/`isEncodedByBioChemEntity` edges and Bioschemas
seeds (residue annotations, AlphaFold/PDB 3D-model representations). Because
schema.org has no forward Gene→Disease/Drug predicate, those cross-entity edges
are emitted under `@reverse` with real predicates, reusing the same curated mesh
edges as the human-facing block (see [mesh.md](mesh.md)). The inline copy caps
every long list (to 15) and points at the sidecar for the full graph, so the
readable lead is never buried.

## Building the corpus

The batch driver builds the whole corpus in three phases, designed so the
cross-entity mesh resolves completely in one pass with no concurrency-unsafe
shared writes:

1. **Collect (parallel).** Each worker resolves + collects one entity, caches
   the bundle under `<dist>/cache/`, and returns its set of resolvable keys.
   Workers touch only their own files — no races.
2. **Merge (single writer).** All key-sets combine into the complete
   `manifest.json`; every resolved forward edge is inverted into
   `reverse_edges.json`; and `evidence.json` is written — per-entity evidence
   scores plus the corpus-wide signal distributions and per-slug component counts
   that power the render-time corpus-relative framing ("top N% corpus-wide") and
   the cross-entity reads (a thin page quoting a related page's numbers).
3. **Render (parallel).** Each worker reloads its cached bundle plus the now
   complete, read-only manifest and writes `page.md` + `entity.jsonld`. No
   re-collect; manifest is read-only — no races.

Parallelism uses a process pool (not threads), because the call log and the
manifest are module globals safe only for process-level parallelism. Failures
are skip-and-log, never fatal. The seed lists themselves are regenerated from
biobtree's own ingest sources (`build_corpus.py`): HGNC for genes, the
signal-ranked Mondo corpus for diseases, and gated ChEMBL molecules for drugs
(dropping non-therapeutic reagents and salt-form children whose parent is also
present).

The everyday loop is wrapped by `atlas.sh` — `test all` builds the dense
reference set and runs the gates; `prod` builds the full corpus. See the
[README](../README.md).

Each build stamps its provenance into every page's frontmatter — `atlas_version`
+ `atlas_commit` (from `git describe`, captured at build start) and the
`biobtree_version` it read from `/ws/meta`. `atlas.sh release` tags and publishes
the **pipeline** only — the corpus is never attached to a release.

## Validation

Three layers guard correctness:

- **Per-page gate (`body_gate`).** A fresh bundle is regression-tested against a
  saved per-entity snapshot: `clean` (identical), `drift` (minor numeric
  movement, or a key declared *shrinkable*), `regression` (a key emptied or a
  count dropped sharply), or `first_run`. A regression refuses to publish; in
  batch it is log-only.
- **Summary gate.** When the optional LLM summary is on, a two-pass adversarial
  judge flags any statement not grounded in the deterministic body; only claims
  both passes flag count, and `--strict-summary` refuses to publish them.
- **Test gates.** Hermetic **unit** tests (transport stubbed, no network) guard
  the logic; **integration** tests run over the *built corpus* — the frozen
  contract (H2 order + the frozen H3/H4 anchor sets), mesh integrity, JSON-LD
  validity, frontmatter schema, and data-quality (no `nan`, no float32 noise, no
  HTML entities, no duplicate rows, no raw ontology ids as labels, no unbalanced
  bold markers). The release rule is: dense build → integration green → only then
  the full corpus.
