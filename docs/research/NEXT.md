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

- ⏸️ **UniProt CC narratives** (FUNCTION / SUBUNIT / SUBCELLULAR LOCATION /
      TISSUE SPECIFICITY / DISEASE / PTM) — audit's #1 gap. Paused on the
      design decision below. 🕓 BIOBTREE_ISSUES #9 (v2 expansion forwarded).
- ⏸️ **Named isoforms** (p53α/β/γ, K-Ras4A/4B, p16γ, titin N2A/N2B/N2BA/novex).
      Same UniProt dependency. 🕓 BIOBTREE_ISSUES #9.
- [ ] **Drug → indication → biomarker triples** in §10 (sotorasib + KRAS-G12C +
      NSCLC; osimertinib + EGFR-T790M; olaparib + BRCAness). Complex — needs
      ChEMBL indication + FDA biomarker list bridging. Defer until V1 corpus.

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

- [ ] **Fix `src/atlas/bench/dataset_coverage.py`** — still references the
      pre-refactor `collect.py` shape (broke during the package
      reorganization).
- [ ] **Tune `body_gate.verdict()` threshold** — the "count dropped >50% →
      regression" rule misfired when biobtree's RefSeq REVIEWED-only refresh
      legitimately dropped TP53 mRNA 46→25. Consider a `schema_changed`
      signal that demotes count drops to `drift`, or a per-key
      expected-direction hint in section metadata.

## Open biobtree dependencies

| Issue | Status | Atlas impact |
|---|---|---|
| #4 silent multi-hop failure | open | Largely subsumed by #1/#3 RESOLVED |
| #6 entry xrefs counts-not-values | open | Not blocking — Atlas uses `map_all` |
| **#9 UniProt CC + reviewed flag + isoforms** | **open — biggest** | Blocks Path C UniProt CC + named isoforms |
| #10 AlphaFold empty for >2700 aa | 🕓 fix tomorrow | §4 currently constructs `AF-<acc>-F1` heuristically; pLDDT missing for ATM/BRCA2/DMD |
| #12 pubchem_activity KRAS gap | 🕓 fix tomorrow | Workaround in place via chembl_activity; no functional gap |
| #13 pharmgkb_guideline / _clinical / _variant empty | open | §10 PharmGKB block can only state existence, not contents; deeper PGx narrative blocked |
| #14 reactome pathway entries with empty `name` | open (just filed) | Disease §14 renders "Unnamed pathway (R-HSA-N…)" for 1-2 pathways per cohort; graceful fallback in place |
| #15 chembl_molecule parent/child salt-form linkage exposed only via `childs` on parent | open (just filed) | Disease §13 pays ~30 extra entry calls per disease to dedupe salt-form drugs (TAMOXIFEN + TAMOXIFEN CITRATE, DOCETAXEL + DOCETAXEL ANHYDROUS); workaround acceptable now, won't scale to drug pages |
| #16 No `list-ids` endpoint for a dataset (corpus enumeration) | open (just filed) | **Blocks all-diseases / all-genes scale-out.** Today we depend on external Mondo .obo / HGNC TSV dumps for full enumeration; brittle and defeats the "biobtree is single source" model. |
| #17 No bulk xref-count check (one /entry per id to filter by signal) | open (just filed) | Even with #16 resolved, filtering ~25k Mondo nodes by Atlas-relevant signal (gwas/civic/clinvar/gencc xrefs) needs 25k entry calls. Suggest adding xref counts to search/list response schemas. |

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

## Important — needs design decision

### UniProt CC narratives: how do we source them?

Audit's #1 gap. biobtree's uniprot `entry` exposes only
`names / alternative_names / sequence / id / name` — no CC narrative. Same
gap is filed as BIOBTREE_ISSUES #9 (v2 forwarded to dev).

Three implementation paths:
1. **Forward to biobtree** (already filed). Most disciplined; slowest.
2. **Direct UniProt REST per accession** — *rejected*: per-gene network
   load + reliability liability in the hot path at full-corpus scale.
3. **Bulk UniProt flat-file parse locally** — download Swiss-Prot
   (~80MB compressed), parse CC blocks once, store per-accession locally.
   No network per gene, rebuilds on UniProt release.

**Lean:** (3) short-term, with (1) as the long-term path once biobtree
ships the fields. Parser pattern lives in `/data/biobtree` (Go) or use
Biopython's `Bio.SwissProt`.

**Update 2026-05-30 (biobtree refresh):** Decision became *less acute*.
Wiring **CIViC's gene-level narrative paragraph** (now live in biobtree —
see Path C "new datasets" above) gives Atlas a real first-paragraph
narrative for cancer-relevant genes today, without #9 being fixed. UniProt
CC remains needed for non-cancer genes (titin, housekeeping, etc.), but
the path-of-least-resistance v1 release can ship with CIViC paragraphs as
the headline narrative for the cancer-genome subset and defer the broader
UniProt CC question.

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

**OPEN — deferred to future iteration:**
- [ ] **LLM executive summaries** for the 18 disease pages (Task #42)
- [ ] **Full 61-disease backlog** (43 more diseases to run; use
      `bin/run_disease_backlog.py` or the Enju workflow at
      `src/atlas/disease/enju.yaml`)
- [ ] **All-diseases scale-out** beyond the curated 61: 🕓 blocked on
      BIOBTREE_ISSUES #16 (no `list-ids` endpoint) + #17 (no bulk
      xref-count check). Today's external-dump workaround is brittle;
      proper fix is upstream so biobtree remains the single source.
- [ ] **Disease body_gate threshold tuning** — first_run vs drift mechanics
      are gene-tuned; disease may need its own thresholds
- [ ] **Disease declarative-lead sentence** (gene side has one; disease
      currently skips it — `assemble_page(bundle=None)` path)
- [ ] **schema.org/MedicalCondition coverage audit** — check whether more
      schema.org fields apply (`epidemiology`, `riskFactor`, `signOrSymptom`)
- [ ] **Slug stability override** — currently slug derives from
      Mondo's canonical_name (e.g. endometrial cancer → "endometrial-carcinoma").
      Caller can pass slug explicitly through the workflow record; document.

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
