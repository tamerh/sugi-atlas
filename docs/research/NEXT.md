# Atlas roadmap — what's left

Synthesized from [`01_landscape_and_ai.md`](01_landscape_and_ai.md),
[`02_biobtree_mining.md`](02_biobtree_mining.md),
[`03_page_audit.md`](03_page_audit.md).

This doc lists only **open items**. Shipped items are recorded in the
`data`-branch git log — search e.g. `git log --grep="Path A"` or
`git log --grep="Path C"`.

Conventions:
- `[ ]` — TODO.
- `⏸️` — paused on a design decision (see bottom).
- `🕓` — waiting on an upstream biobtree fix (see [`BIOBTREE_ISSUES.md`](../BIOBTREE_ISSUES.md)).

---

## Path A — AI-friendliness envelope

Remaining (all cross-repo, deferred to launch):
- [ ] `/llms.txt` → drafted at `docs/site-drafts/llms.txt`, copy to
      `biobtree-content/static/` at launch.
- [ ] `robots.txt` allowlist → drafted at `docs/site-drafts/robots.txt`, copy at launch.
- [ ] Per-gene `<lastmod>` sitemap → Hugo `config.toml` tweak in `biobtree-content`.
- [ ] `Last-Modified` HTTP header from `generated_at` → Hugo theme + nginx pass-through.

## Path B — Provenance moat

- [ ] **Per-fact HTML anchors** on the page (deep-linkable single facts).
      Defer until we see how AI clients use page-level + section-level
      anchors first.
- [ ] **`<link rel="alternate" type="application/ld+json">`** discovery
      hints in `<head>` for `entity.jsonld` / `provenance.json` / `bundle.json`
      → Hugo theme work (pre-launch, cross-repo).

## Path C — Content depth

### From the page audit (biggest user-visible gaps)

**SHIPPED 2026-05-31** — top-of-audit content gaps:
- UniProt CC narratives + named isoforms (`git log --grep="BIOBTREE #9"`) —
  function / subunit / subcellular_location / tissue_specificity / ptm /
  disease / cofactor / domain narratives + named isoforms table; declarative
  lead + JSON-LD use FUNCTION CC; disease §5 cohort_function_summary.
- Drug × variant × indication precision-medicine triple via CIViC
  (`git log --grep="audit gap #6"`) — gene §10 + disease §13 carry the
  clinical-evidence table at CIViC Level A→E, deduped by association,
  with per-row PubMed/CIViC links. Cancer subset only; PGx slice
  (CYP × dosing — not in CIViC's scope) follows later via PharmGKB
  once biobtree #13 lands.
- [ ] **ClinGen dosage sensitivity verdicts** (audit Top-10 #4) —
      haploinsufficiency / triplosensitivity score per gene. Public REST
      API at `https://search.clinicalgenome.org/`; one call per gene.
      Network-in-collect concern same as Open Targets; lean: file as a
      biobtree ingest request.
- [ ] **NCBI RefSeq summary paragraph** (audit Top-10 #5) — the
      `Entrezgene_summary` field on NCBI Gene; an independent curated
      narrative often complementary to UniProt FUNCTION. E-utilities efetch
      pull, one call per gene. Same network concern + same biobtree-ingest
      lean as ClinGen above.
- [ ] **Ensembl Compara gene tree** (audit Top-10 #10) — 200+ species
      orthologs + paralog tree + duplication/speciation node labels.
      Atlas's §5 today surfaces ~4 ortholog rows from biobtree; the full
      gene tree is a ~30 GB Compara source. Cleanest path: biobtree ingest
      → `>>ensembl>>gene_tree` edge → existing §5 collector renders the
      summary. User pause direction: hold for now.
- [ ] **Regulatory build features** (cross-gene weakness #5 in audit) —
      Ensembl Regulatory Build for promoters / enhancers / CTCF / TFBS
      *around* the gene region. Atlas has JASPAR motifs (curated
      transcription-factor binding) but not the per-gene regulatory
      annotations Ensembl ships. Would need a new biobtree dataset.
- [ ] **GeneRIFs** (cross-gene weakness #9 in audit) — NCBI's per-gene
      paper-anchored one-line claims. Thousands per popular gene; high
      value for RAG / citation-aware LLMs. Big literature corpus; needs
      a new biobtree dataset (NCBI provides bulk GeneRIF dumps).

### From the biobtree mining (still open)

- [ ] **chembl_molecule drilldowns** per phased drug (per-drug
      `chembl_activity` / `clinical_trials` / `pubchem` / `chebi`). Turns §10's
      "Molecules" list from "names + phase" into structured records.
      Per-drug fan-out cost grows with phased-drug count — defer until
      asked for.
- [ ] **PMID per ChEMBL activity row.** The chembl_document enrichment we
      shipped attaches title+journal to the **assay sample** rows (3 per page).
      Adding it to every `chembl_activities` row would need the same
      `chembl_assay>>chembl_document` fetch per activity — bounded by the
      `chembl_activity_potent_count` cap (≤30 rendered). Acceptable cost,
      but body-table noise: the same paper recurs across many activities for
      a target. Defer until we have a use case.
- [ ] **`literature_mappings` → PMID resolution.** chembl_document entry
      carries `literature_mappings|1` xref. Route to actual PMID still TBD;
      probe `>>chembl_document>>literature_mappings` next time it's needed.
- [ ] **`brenda_kinetics` / `brenda_inhibitor`.** Empty for all reference
      genes via every probed route. Either genuinely sparse in biobtree or
      route TBD — revisit on next biobtree refresh.
- [ ] **`civic_variant` / `civic_evidence` / `civic_assertion`.** Return
      n=0 via hgnc/uniprot — route TBD. CIViC's gene-level narrative is
      already wired; variant-level would be the next depth tier.
- [ ] 🕓 **Drug patent landscape** (BIOBTREE #25 / #26 / #27). Today the drug
      §1 shows the SureChEMBL mention total + the dominant matched-structure
      share (e.g. Imatinib 109,520; one structure = 98%). A real landscape is
      blocked on three biobtree gaps: assignee + CPC/IPC technology class
      (#25 — attributes documented but not populated), distinct patent
      **families** (#26 — `patent_compound` has no `patent_family` rollup; the
      honest dedup of mention inflation), and jurisdiction / filing-timeline /
      recent-patents (#27 — `>>patent_compound>>patent` is id-ordered with no
      date sort or facets, so any sample is biased). When those land: restore
      a dedicated **`## Patent landscape`** section — distinct families, top
      assignees, CPC technology areas, recent filings, jurisdiction mix. The
      collector already captures `patent_compound_breakdown` (s11); the
      `chembl_molecule>>patent_compound>>patent` hop is verified working, so
      this is a render + bounded-fetch task once the data/facets exist.

### 2026-05-30 biobtree refresh — all wirables SHIPPED

The 117 → 139 dataset-node refresh contributed seven new sections of §-level
content; all shipped end-to-end (collect → render → provenance → body_gate →
snapshot → dist) and visible on the six reference-gene pages. Per `git log`:

- intOGen + CIViC cancer-significance overview (cancer-narrative paragraph
  between page lead and Summary; intOGen driver role + cancer types)
- RNAcentral (§1 — closes the lncRNA/miRNA gap MALAT1/XIST/HOTAIR/MIR21/NEAT1/H19)
- BRENDA EC annotation (§3 — BRCA1=2.3.2.27 E3, KRAS=3.6.5.2 GTPase, TTN=2.7.11.1 kinase)
- ChEMBL assay depth + type breakdown + 3 samples (§10)
- Cellosaurus cell-line resources with category breakdown (§10)
- patent_compound per-molecule count (§10 molecules table)
- chembl_document source paper per assay sample (§10 assay-sample table)

## Hygiene fixes

All hygiene items from this section shipped 2026-05-31. Per `git log`:
- `dataset_coverage.py` rewrite (walks the REGISTRY for both entity types,
  no source-text scraping) — commit `155dbf3`
- `body_gate.verdict()` per-key shrink allow-list via the new
  `Section.shrinkable` field; real-world hints declared on gene §2/§3 —
  commit `155dbf3`

## Open biobtree dependencies

| Issue | Status | Atlas impact |
|---|---|---|
| #4 silent multi-hop failure | open | Largely subsumed by #1/#3 RESOLVED |
| #6 entry xrefs counts-not-values | open | Not blocking — Atlas uses `map_all` |
| ~~#9 UniProt CC + reviewed flag + isoforms~~ | **✅ RESOLVED 2026-05-31** | `comments` block + named `isoforms` + `is_canonical` flag now in uniprot entry. Atlas-side wiring pending — see new Path C item below. |
| ~~#10 AlphaFold empty for >2700 aa~~ | **✅ RESOLVED 2026-05-31** (biobtree extended coverage; remaining empties are AlphaFold-DB upstream gaps for very-large proteins >~3000 aa — Atlas adds graceful footnote) |
| ~~#12 pubchem_activity gap~~ | **✅ RESOLVED 2026-06-01** — KRAS now returns 4927 rows; Atlas's stale workaround comment removed. |
| #13 pharmgkb_guideline / _clinical / _variant empty | 🟡 **partial 2026-06-01** — `pharmgkb_clinical` ✓ and `pharmgkb_variant` ✓ now return data (wired into §10); `pharmgkb_guideline` still empty. |
| ~~#14 reactome pathway entries with empty `name`~~ | **✅ RESOLVED 2026-06-01** — all probed pathways carry names; Atlas's "Unnamed pathway" fallback removed. |
| ~~#15 chembl_molecule parent/child salt-form linkage~~ | **✅ RESOLVED 2026-06-01** — both directions work: `>>chembl_molecule>>chembl_moleculeparent` (child→parent) and `>>chembl_molecule>>chembl_moleculechild` (parent→children). Unblocks drug entity at scale. |
| ~~Mondo OBO cross-ontology xrefs + UBERON anatomy~~ | **✅ RESOLVED 2026-06-01** — `>>mondo>>{doid,sctid,umls,ncit,medgen,icd10cm,icd11,gard,meddra,nord,uberon}` all work. §1 federated identifier table extended; JSON-LD `sameAs` + `code` + `associatedAnatomy` populated. |
| ~~#16 list-ids endpoint~~ | **retracted 2026-05-31** — corpus enumeration belongs upstream (HGNC TSV, Mondo OBO, ChEMBL releases); biobtree shouldn't duplicate. |
| ~~#17 bulk xref-count check~~ | **retracted 2026-05-31** — per-entry xref counts work fine, only paid once per release with local cache. No real bottleneck. |
| #25 patent attrs (assignee/CPC/IPC) not populated | open 🕓 | Blocks drug patent **assignee** breakdown + **technology-class** (CPC/IPC) landscape — fields documented in `patent.md` but empty across CN/EP/US/WO. |
| #26 no `patent_family` rollup on `patent_compound` | open 🕓 | Blocks **distinct-family** count (honest dedup of SureChEMBL mention inflation); only reachable by `entry()`-ing 100k+ patents otherwise. |
| #27 `>>patent_compound>>patent` id-ordered, no date sort / facets | open 🕓 | Blocks **jurisdiction / timeline / recent-patents** — any bounded sample is a sampling artifact. |
| #28 `collectOntologyIDs` over-matches common terms → contaminated `mondo→{clinical_trials,intogen,civic}` | open 🕓, dev fixing (exact-match, needs re-index) | **⚠ HAS AN ACTIVE ATLAS MITIGATION TO REVERT.** Now: s13 title-validates trials (commit 43c5736) — too-strict proxy (cardiomyopathy 317→92). **On resolve:** (1) swap s13 title-match → exact condition-match (self-deactivating; needs `conditions` added to the trial map compact_fields — dev offered) or remove it; (2) re-check disease §13 civic / cohort CIViC route for residual contamination. **Blast-radius check done (2026-06-02): intogen blast radius = 0 (Atlas reaches it gene-first via `hgnc→intogen`, never the contaminated `mondo→intogen`); `mondo→civic_evidence` ~0 contamination on the head (294/303 rows match canonical; the 9 others are real glioma evidence) → no interim civic guard, wait for upstream.** Re-check if the corpus expands into the rare-disease tail before #28 lands. |

## Page output contract — web team requests (2026-06-02)

The frontend does passage/chunk indexing (Google, Perplexity, Pagefind), so the
*output contract* (H2 set, anchor IDs, frontmatter schema) is an SEO/AI lever, not
just cosmetics. **P0 is a freeze-before-launch decision** — changing the H2 set or
anchor IDs later 404s every deep link in the wild. Sequenced as the web team
proposed (P0 contract → P1 lede/links → P2 schema → P3 polish).

- [x] **P0 · Canonical H2 taxonomy ✅, same set + order on every entity, emit even if
      empty** (placeholder = informative, e.g. "non-coding RNA — no protein
      product", not bare "no data"). *Design + freeze jointly before launch.*
      **Reconcile with layer-B:** the web team's flat 6–7 H2 set (Summary /
      Identifiers / Gene structure / Protein / Function / Disease & clinical /
      Drugs) supersedes the current 3-zone nesting (Gene-locus / Protein-product /
      Clinical) — and it also fixes the "fewer sections but each too dense"
      feedback. The gene↔protein split survives as sibling H2s (Gene structure /
      Protein / Function); layer-A typed-Protein JSON-LD (`@id`/encodes) is
      independent of body nesting and stays. Need parallel frozen sets for
      disease + drug. **Mostly our work (render order + headers + placeholders).**
- [x] **P0 · Stable explicit anchor IDs ✅** via `## Heading {#kebab-id}` (goldmark
      `parser.attribute.title=true` — web team to enable). Backend-owned IDs
      decoupled from prose, locked across regens (`#summary #identifiers
      #transcripts #expression #variants #protein #function #drugs #trials
      #generif #pathways`). Replaces the current `<a id>` HTML anchors. Freeze
      with the taxonomy.
- [x] **P1 · `## Summary` H2 as the lede ✅** — wrap the existing intro region
      (lead + RefSeq + At-a-glance) in one `## Summary {#summary}` so it's a
      passage-indexable chunk. ⚠ Reconcile with the earlier "intro is *not* a
      section" preference — one Summary section wrapping the whole intro (not
      peer headers) satisfies both; confirm before doing.
- [x] **P2 · `tldr:` frontmatter ✅** — 3–5 structured key-fact bullets from the
      bundle (same source as At-a-glance) → frontend "Key facts" box + clean
      snippets.
- [x] **P2 · Decouple symbol/title ✅ (added typed `identifier`; `symbol`=slug kept) + canonical-case** — `symbol` is the slug for
      genes but a kebab slug for diseases (`ataxia-telangiectasia`), so templates
      can't trust `.Params.symbol`. Make a typed identifier field (HGNC / MONDO /
      ChEMBL id) + `title` = display. Drug titles already de-SHOUTed (#12).
- [x] **P2 · Search-alias frontmatter ✅ (`alt_names`)** — prev/alias symbols (genes), xrefs
      (diseases), synonyms/brands (drugs) for Pagefind ("Her2" → ERBB2).
      ⚠ **`aliases:` is Hugo-reserved (generates URL redirects)** — use a
      different field name (`alt_names:` / `synonyms:`), NOT `aliases:`.
- [x] **P3 · `section_defaults:` ✅ open/collapsed hints** — skippable for v1; the
      frontend can apply its own rules once the H2 set + IDs are frozen.
- [x] **P3 · Emit `index.md` not `page.md` ✅ ready via `ATLAS_PAGE_FILENAME` (default page.md; flip on web-team go)** — one-liner in `batch.render_one`;
      lets `biobtree-content` mount the dist as a Hugo module (no sync wrapper).
      Coordinate timing (breaks their current sync script).
- [ ] 🕓 **DEFERRED — P1 · Per-table source link ("showing 40 of 312 · view all
      →")**. Filed for separate structural design (per 2026-06-02 decision). Why
      it's not a quick win: needs a uniform per-table `(total, canonical
      source-URL)` contract, and **neither is uniform** — totals aren't reliably
      exposed for every table (ties into biobtree #6 entry-counts) and several
      sources lack a clean per-entity "view all" deep link (patents #25–27,
      BindingDB, BRENDA). Also must align with our existing "section-level source
      link only, no per-row links" convention. Design the per-dataset
      (total + deep-link) map once, then apply uniformly — don't bolt on ad-hoc.

## Pre-launch / cross-repo work (defer to launch day)

These live in **`biobtree-content`**, not this repo. Sending AI bots toward
the site too early caches an incomplete corpus.

- [ ] Copy `docs/site-drafts/llms.txt` → `biobtree-content/static/llms.txt`
- [ ] Copy `docs/site-drafts/robots.txt` → `biobtree-content/static/robots.txt`
- [ ] `curl -I https://sugi.bio/{llms,robots}.txt` to verify
- [ ] Hugo `config.toml`: `<lastmod>` per-page sitemap
- [ ] Hugo theme: `<link rel="alternate" type="application/ld+json">` hints
- [ ] Hugo theme: `Last-Modified` HTTP header from `generated_at`
- [ ] After deploy: resubmit sitemap to Google Search Console; watch
      `chatgpt_report.py` for new AI-fetch patterns

## Disease entity — status

**SHIPPED 2026-05-30:**
- 14-section deterministic collector (§1–§14) + 3 render-only derived
  views (§15 drug_repurposing, §16 druggability_pyramid, §17 undrugged_target_profiles)
- `DiseaseAnchors` + `resolve(name|mondo_id)` with 4-route gene-cohort
  union (GWAS, GenCC, ClinVar, CIViC-evidence) capped at top-50 by
  evidence-route count
- `cohort.fan()` helper that reuses 9 of the gene §-collectors over the
  disease cohort — massive code reuse, no duplication
- Render parity with gene side; shared `atlas.render_common.table()`
- Pipeline integration (`run_disease()` in atlas.pipeline)
- Workflow + 4 Enju task scripts (collect_render / body_gate / summary / publish)
- schema.org/MedicalCondition JSON-LD sidecar + provenance.json sidecar
  (per-section dataset + chain + upstream URL trail)
- 18 disease pages built deterministically (dev backlog; full 61-disease
  backlog can resume any time via `python -m atlas.disease.corpus run`)
- 9-item polish pass: §1 empty-row omission + monarchinitiative Mondo URL,
  §4 CIViC link + dual-evidence enrichment, §13 true-set phase distribution
  + parent/child salt-form drug dedupe, §14 unnamed-pathway fallback,
  §16 "+N more" overflow indicator

**ALSO SHIPPED 2026-05-31** (per `git log`):
- Disease declarative-lead sentence (`atlas.page.disease_declarative` —
  Mondo-anchored, with cohort + GWAS + somatic + ClinVar + trials clauses,
  dominant-pathway clause, top-3 drug names)
- schema.org/MedicalCondition gained the `drug` property (top-5
  disease-scoped trial drugs as shallow Drug nodes with ChEMBL URLs)
- Slug-override pattern documented in `atlas/disease/__init__.py`
- Disease body_gate now uses the same `shrinkable` mechanism via the
  shared `Section` dataclass (gene §2/§3 declare hints; disease sections
  can declare theirs as fluctuations appear)

**OPEN — deferred to future iteration:**
- [ ] **LLM executive summaries** for the 18 disease pages (Task #42 — paused per dev-phase direction)
- [ ] **Full 61-disease backlog** (43 more diseases to run; use
      `python -m atlas.disease.corpus` or the Enju workflow at
      `src/atlas/disease/enju.yaml`)
- [ ] **All-diseases scale-out** beyond the curated 61. Approach:
      Atlas-side discovery script parses Mondo's OBO (from obofoundry.org)
      to enumerate all MONDO IDs + canonical names; optionally filters by
      branch (e.g. only `human disease` descendants) to drop irrelevant
      ontology nodes. Same pattern for genes via HGNC's
      `hgnc_complete_set.txt`. Feeds the resulting list into the existing
      `python -m atlas.disease.corpus` driver. No biobtree changes needed —
      this is normal downstream-consumer pattern. (Earlier framing as
      "blocked on a biobtree list-ids endpoint" was retracted on
      reflection; see BIOBTREE_ISSUES.md Retracted section.)
- [ ] **Schema.org MedicalCondition additional fields** —
      `epidemiology` (needs Orphanet prevalence numbers — biobtree exposes
      Orphanet xref but not prevalence body); `signOrSymptom` (HPO
      phenotype-disease links); `associatedAnatomy` (UBERON anatomy
      mapping); `riskFactor`; `cause` (etiology). None blocked on biobtree
      core — would need new ingest sources.

## Out of scope for now

- Open Targets-style genetic associations (L2G scores) — would require
  their GraphQL API. Revisit after Path A+B+C land.
- Multi-language pages — defer until corpus is committed.
- Drug entity — same Section/anchors pattern; tackle once disease+gene
  pages have run at scale.

## V1 release-ready page checklist

In-repo items all satisfied for the 6 reference genes
(TP53/BRCA1/CDKN2A/KRAS/TTN/EGFR). Remaining are cross-repo Hugo work:

- [ ] discoverable in sitemap with `<lastmod>`
- [ ] `<link rel="alternate" ...>` discovery hints in page head
- [ ] `Last-Modified` HTTP header from `generated_at`
