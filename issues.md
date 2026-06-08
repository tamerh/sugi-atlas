# Sugi Atlas — backend / pipeline issues

Living tracker for the **sugi-atlas generator** (the corpus pipeline) — collectors,
renderers, the mesh, the page contract. Web-side issues live in the web repo's
`issues.md`; upstream data gaps are biobtree's. Cross-reference when relevant.

**Priority:** P0 freeze-before-launch · P1 launch-blocker · P2 important · P3 nice-to-have.
Close = delete the entry; keep the file to what's still open.

---

## Content & rendering

### P1 · Thin disease pages — condense the data-less ~25%
**Status:** open · **Priority:** P2

~4,601 disease pages (25%) have identifiers + ontology placement but **no molecular
evidence** — no associated gene, variant, trial, or drug exists for them in the
integrated sources (verified: `MONDO:0009056` "cutis verticis gyrata and intellectual
disability" has 5 xrefs, none to a gene). They render all six H2 zones, each with a
"No X" placeholder, so they read as broken/empty.

Decide presentation: keep the honest all-zones layout, or **condense** the empty
molecular zones into one cross-reference-stub line (identifiers + synonyms +
OMIM/Orphanet links + "no curated molecular evidence in integrated sources"). Leaning
condense — make it read as a deliberate stub, not a broken page. The page is genuinely
useful as a cross-reference hub; only the empty-zone clutter is the problem. Not a bug —
purely presentation. *(parked from the v1.1.4 session)*

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
