# Sugi Atlas — backend / pipeline issues

Living tracker for the **sugi-atlas generator** (the corpus pipeline) — collectors,
renderers, the mesh, the page contract. Web-side issues live in the web repo's
`issues.md`; upstream data gaps are biobtree's. Cross-reference when relevant.

**Priority:** P0 freeze-before-launch · P1 launch-blocker · P2 important · P3 nice-to-have.
Close = delete the entry; keep the file to what's still open.

---

## Content & rendering

### P1 · Thin disease pages — enrich first, then decide on the remainder
**Status:** partly addressed · **Priority:** P2

~4,601 disease pages (25%) looked empty — all six H2 zones as "No X" placeholders.
But "empty" was partly **us not surfacing data we had**: clinical features were
pulled from Orphanet only, missing OMIM→HPO (and Mondo→HPO) annotations. Fixed in
`74ad9e8` (v1.1.6) — `MONDO:0009056` went 0→2 HPO features; corpus-wide lift TBD.

**Before any condensing, re-check the same way for other under-surfaced data:**
- **OMIM→gene** — does the cohort miss Mendelian genes that hang off the OMIM entry
  rather than GenCC/ClinVar? (cutis's OMIM had no gene, but others may.)
- Re-measure after the v1.1.6 regen: how many pages are *still* genuinely empty
  (no HPO, no gene, no anything)?

Only for that genuinely-data-less remainder, decide presentation: keep the honest
all-zones layout, or **condense** into a cross-reference stub. Not a bug — the
lesson (from the v1.1.5/6 session) is to surface what we have *before* hiding zones.

---

### P4 · Non-oncology approved indications logged at phase 3
**Status:** open · **Priority:** P3

The approved-vs-in-trials tiering recognises anticancer (ATC L01/L02) approvals that ChEMBL
logs at phase 3 (imatinib→CML). **Beyond oncology** the data can't resolve "approved for
*this* disease" — a drug approved for e.g. epilepsy but logged phase 3 reads honestly as
"in trials" (conservative, not false). PubChem confirms molecule approval, not per-disease.
Future: cross-check another source (FDA labels / DrugCentral / openFDA) to widen the
per-disease approval signal beyond oncology. Honest as-is — an enhancement, not a fix.

---

### P5 · OMIM→gene as a cohort source (+ a broader "what are we under-surfacing?" audit)
**Status:** open · **Priority:** P3 · **feasibility checked**

Came out of the v1.1.6 session: the HPO win (clinical features were hiding in OMIM, not
just Orphanet) suggests we may under-surface **genes** too. The disease gene cohort is
built from GenCC + ClinVar + GWAS + CIViC — **not** OMIM's phenotype→gene link.

**Feasibility (measured):** `>>mim>>uniprot>>hgnc` resolves (CF 219700 → CFTR). But the
coverage gain is small: of 40 sampled thin (gene_count 0) diseases with an OMIM id, only
**3 (~7.5%)** gain a gene → ~**96 of 1,284** corpus-wide. Most thin diseases are
*genuinely geneless* old OMIM clinical syndromes (cutis-style), so OMIM correctly adds
nothing there. Real but long-tail.

**Why it's parked, not dropped:** even a small addition is worth it (those ~96 pages go
0→1 gene → full molecular sections unlock), but the cohort is the page's **core** — a new
source must integrate with the evidence ranking (OMIM gene = curated, rank near GenCC),
apply to non-thin diseases too (does it add genes they miss?), and be well-tested. Deserves
a comprehensive design pass, not a rushed add.

**Broader framing:** do a deliberate audit of *all* under-surfaced data — the pattern
(HPO via OMIM, genes via OMIM) suggests other sources may be sitting unused behind a
secondary xref. Worth one systematic pass rather than ad-hoc discovery.

---

### P3 · Score-ranked "firehose" tables need a threshold, not a row cap
**Status:** open · **Priority:** P3

The curated tables now cap at `ROW_CAP = 100` (gene/disease/drug render). But a
family of tables is score/significance-ranked with totals in the thousands, where
*any* row cap is arbitrary (row 400 by raw rank is meaningless): **STRING** (EGFR:
11,600), **IntAct**, gene/disease **GWAS associations**, disease **tiered variants**,
**CollecTRI** targets (master regulators), **ESM2/Diamond** homologs, **SpliceAI**,
**AlphaMissense**, cohort/expression **tissues**. These keep their literal caps for
now (40/30) with honest `*+N more*` disclosure.

The real fix is a **meaningful threshold**, then a generous bounded cap:
- STRING → combined score ≥ 700 (high confidence; ≥900 highest)
- IntAct → MI confidence threshold
- SpliceAI → Δscore ≥ 0.5 · AlphaMissense → ≥ 0.564 (likely-pathogenic)
- GWAS → p ≤ 5e-8 (genome-wide significance)

Mostly a render-side filter (score is in the bundle) + a band-count line
("412 interactions: N≥900, M≥700…"). STRING is the flagship to design first.
Not urgent — flagged during the v1.1.8 consistency pass.

---

## To consolidate (currently in the web tracker)

The web repo's `issues.md` still has a "sugi-atlas (generator / pipeline)" section —
those are backend items that belong here. When convenient, move the still-open ones
in and trim the web tracker to a pointer:
- **Per-table source-link contract** (P1, blocked on biobtree per-dataset totals)
- **Per-section default-open hints** (P3, optional)
- *(already done — drop from the web tracker:* `description` frontmatter ✓, `index.md` output ✓*)*

---

## How to use this file

- One issue per entry; keep the title scannable.
- Priority: P0 freeze-before-launch · P1 launch-blocker · P2 important · P3 nice-to-have.
- When closed, delete the entry.
- Cross-reference web / biobtree issues when relevant.
