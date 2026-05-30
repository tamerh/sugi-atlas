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

Decide before resuming. Becomes the next Path C item by impact.

## Out of scope for now

- Open Targets-style genetic associations (L2G scores) — would require
  their GraphQL API. Revisit after Path A+B+C land.
- Multi-language pages — defer until corpus is committed.
- Drug / disease entities — same shape as gene; do after gene pipeline is
  battle-tested at 100+ genes.

## V1 release-ready page checklist

In-repo items all satisfied for the 5 reference genes
(TP53/BRCA1/CDKN2A/KRAS/TTN). Remaining are cross-repo Hugo work:

- [ ] discoverable in sitemap with `<lastmod>`
- [ ] `<link rel="alternate" ...>` discovery hints in page head
- [ ] `Last-Modified` HTTP header from `generated_at`
