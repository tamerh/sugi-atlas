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
