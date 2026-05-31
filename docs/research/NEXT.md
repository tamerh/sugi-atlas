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

UniProt CC narratives + named isoforms — **SHIPPED 2026-05-31** (search
`git log --grep="BIOBTREE #9"`). Closes the audit's #1 gap: gene §3 now
carries `function/subunit/subcellular_location/tissue_specificity/ptm/
disease/cofactor/domain/...` paragraphs + isoforms table; declarative
lead + JSON-LD use FUNCTION CC; disease §5 has cohort_function_summary
("CFHR1: complement regulation; CFB: catalytic C3/C5 convertase; ...").

- **Drug → indication → biomarker triples** in §10 (sotorasib + KRAS-G12C +
      NSCLC; osimertinib + EGFR-T790M; olaparib + BRCAness) — **decided: CIViC.**
      Probe (2026-05-31) found the triple is *already in biobtree*, fully
      normalized: `civic_evidence` (dataset 754) carries
      `molecular_profile | disease | therapies | evidence_type | evidence_level
      | significance` with xrefs to hgnc/mondo/chembl_molecule/pubmed. Reachable
      today via `>>hgnc>>civic_evidence` (gene §10) and `>>mondo>>civic_evidence`
      (disease §13). The earlier Open Targets / FDA-label / biobtree-feature-request
      analysis was a false alarm caused by trying the wrong route
      (`hgnc>>civic_variant`, n=0) — all rejected (wrong layer; data was in
      `civic_evidence` all along). Cancer subset ships in the §10 collector
      (this commit). PGx slice (CYP × dosing — not in CIViC) follows later via
      PharmGKB once biobtree #13 lands. The deleted `SPEC_drug_indication_biomarker.md`
      held the superseded 4-option analysis; full probe trail is in git history.
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
| #12 pubchem_activity gap | **🟡 mostly resolved 2026-05-31** — EGFR/AKT1/NR3C1/AR/ESR1/TP53/CAMK2A all populated. KRAS specifically still empty; upstream root-cause found, fix ETA tomorrow or late afternoon. Atlas's chembl_activity workaround stays in place. |
| #13 pharmgkb_guideline / _clinical / _variant empty | 🕓 fix in tomorrow's release | §10 PharmGKB block can only state existence, not contents; deeper PGx narrative blocked |
| #14 reactome pathway entries with empty `name` | 🕓 fix in tomorrow's release | Disease §14 renders "Unnamed pathway (R-HSA-N…)" for 1-2 pathways per cohort; graceful fallback in place |
| #15 chembl_molecule parent/child salt-form linkage | 🟡 partial 2026-05-31 — child entries now expose `parent` field; remaining ask is a forward map edge. Atlas's per-entry workaround acceptable at disease scale, won't at drug-page scale |
| ~~#16 list-ids endpoint~~ | **retracted 2026-05-31** — corpus enumeration belongs upstream (HGNC TSV, Mondo OBO, ChEMBL releases); biobtree shouldn't duplicate. |
| ~~#17 bulk xref-count check~~ | **retracted 2026-05-31** — per-entry xref counts work fine, only paid once per release with local cache. No real bottleneck. |

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
  backlog can resume any time via `/tmp/run_disease_backlog.py`)
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
      `bin/run_disease_backlog.py` or the Enju workflow at
      `src/atlas/disease/enju.yaml`)
- [ ] **All-diseases scale-out** beyond the curated 61. Approach:
      Atlas-side discovery script parses Mondo's OBO (from obofoundry.org)
      to enumerate all MONDO IDs + canonical names; optionally filters by
      branch (e.g. only `human disease` descendants) to drop irrelevant
      ontology nodes. Same pattern for genes via HGNC's
      `hgnc_complete_set.txt`. Feeds the resulting list into the existing
      `bin/run_disease_backlog.py` driver. No biobtree changes needed —
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
