# Atlas roadmap — derived from research reports

Synthesized from [`01_landscape_and_ai.md`](01_landscape_and_ai.md),
[`02_biobtree_mining.md`](02_biobtree_mining.md),
[`03_page_audit.md`](03_page_audit.md).

Conventions:
- `[x]` items have shipped; commit SHA appended where known.
- `[ ]` items are still TODO.
- `⏸️` items are paused on a decision (see "Important — needs decision" below).
- `🕓` items wait on an upstream biobtree fix (see `docs/BIOBTREE_ISSUES.md`).

---

## Path A — AI-friendliness envelope (do first)

Cheap. Decisive lead vs everyone. Without these, content depth goes
undiscovered by AI crawlers.

- [x] **schema.org `Gene`/`BioChemEntity` JSON-LD** in every page `<head>`,
      with `sameAs` to NCBI Gene / UniProt / Ensembl / HGNC / OMIM. Inline
      `<script>` block at top of body + `entity.jsonld` sidecar. Shipped
      `e61bf9e`.
- [x] **Declarative first sentence** on every page — deterministic, sourced
      from anchors (symbol + name + location + class + canonical UniProt).
      Shipped `14b84df`.
- [x] **Visible "Updated: YYYY-MM-DD"** at the top of every page. Shipped
      with the declarative-lead commit. *(Server-side `Last-Modified` HTTP
      header still pending — Hugo theme work, see Pre-launch section.)*
- [ ] **`/llms.txt`** — drafted at `docs/site-drafts/llms.txt`; copy into
      `biobtree-content/static/` at launch (see Pre-launch section).
- [ ] **`robots.txt` allowlist** — drafted at `docs/site-drafts/robots.txt`;
      copy at launch.
- [ ] **Per-gene `<lastmod>` sitemap** — Hugo config tweak in
      `biobtree-content/config.toml`; pre-launch.

## Path B — Provenance moat (do second)

Establishes Atlas's unique position: the only gene resource where an agent can
cite a single fact with its upstream source.

- [x] **`provenance.json` sidecar** in `dist/atlas/gene/<SYM>/`. schema.org
      `Dataset` mapping each section to its biobtree chains + upstream
      sources (NCBI/UniProt/EBI/etc.) + anchors record. Shipped `6a592a6`.
- [ ] **Per-fact HTML anchors** on the page. Each numeric claim gets an `id`
      so AI agents can deep-link a single fact. Defer until we see how AI
      clients actually use page-level + section-level anchors first.
- [ ] **`<link rel="alternate" type="application/ld+json" href="entity.jsonld">`**
      and similar for `provenance.json` / `bundle.json` so machine clients
      discover sidecars without scraping. Hugo theme work — pre-launch.

## Path C — Content depth

From the biobtree mining + page audit. Each is independently shippable.

### From the page audit (biggest user-visible gaps)

- ⏸️ **UniProt CC narratives** (FUNCTION / SUBUNIT / SUBCELLULAR LOCATION /
      TISSUE SPECIFICITY / DISEASE / PTM) — *audit's #1 gap*. Paused on
      design decision. The expanded ask was forwarded to biobtree dev as the
      v2 of BIOBTREE_ISSUES #9 (now scoped to: CC blocks + `reviewed` flag
      + named isoforms on the uniprot entry payload). 🕓 BIOBTREE_ISSUES #9.
- ⏸️ **Named isoforms** (p53α/β/γ, K-Ras4A/4B, p16γ, titin N2A/N2B/N2BA/novex).
      Source from UniProt's `ALTERNATIVE PRODUCTS` section. Same dependency
      as UniProt CC. 🕓 BIOBTREE_ISSUES #9.
- [ ] **Drug → indication → biomarker triples** in §10 (sotorasib + KRAS-G12C
      + NSCLC; osimertinib + EGFR-T790M; olaparib + BRCAness). Sources to
      bridge: ChEMBL indication + FDA biomarker list. Complex; deferred.

### From the biobtree mining (curl-verified, ready to wire)

- [x] **mirdb** via `>>hgnc>>refseq>>mirdb` — top-30 by max_score, §9
      (post-transcriptional regulators alongside CollecTRI TFs). Shipped
      `989d5d0`. Per-gene: TP53=110, BRCA1=121, CDKN2A=51, KRAS=327, TTN=217.
- [x] **pubchem_activity** via `>>uniprot>>pubchem_activity` — Active rows
      sorted by potency, with clickable CID/AID links derived from the
      activity_id (PMID enrichment per row deferred — entry fetches would
      be expensive at scale). Shipped `2f5ce92`.
- [x] **chembl_activity** (added when KRAS PubChem gap was discovered;
      `BIOBTREE_ISSUES #12`). pChembl-ranked, top 30 at pChembl≥5. Closes
      the KRAS gap on the Atlas side (5,239 activities → 4,825 potent).
      Shipped `189f98a`.
- [x] **ctd_gene_interaction** via `>>hgnc>>entrez>>ctd_gene_interaction` —
      Homo sapiens only, top 30 by PubMed-support count, CV verbs surfaced
      (`increases^expression`, etc.), clickable CTD chemical pages. Shipped
      `c9d4cc2`. Per-gene: TP53=1242, BRCA1=189, CDKN2A=193, KRAS=135, TTN=43.
- [x] **MeSH disease descriptors** via `>>hgnc>>clinvar>>mondo>>mesh` (∪ via
      gencc) — Main descriptors first, tree numbers split, clickable MeSH
      UI links. Shipped `3ef39bd`. Per-gene: TP53=35, BRCA1=16, CDKN2A=11,
      KRAS=20, TTN=39.
- [x] **Anatomy/cell-type names + IDs** — `bgee_evidence` row id encodes the
      UBERON/CL id; we now surface that id (with bioregistry.io link) next
      to each tissue name. Shipped `96155b3`. *(The mining report's
      original "replace bare IDs" was already moot — bgee_evidence had
      shipped names from day one; we repurposed the item to add the
      federated-identity anchor.)*
- [ ] **chembl_molecule drilldowns** per phased drug (per-drug
      `chembl_activity` / `clinical_trials` / `pubchem` / `chebi`). Would
      turn §10's "Molecules" list from "names + phase" into structured
      records. Adds per-drug fan-out (cost grows with phased-drug count) —
      defer until we know it's wanted.

## Hygiene fixes

- [ ] **Fix `src/atlas/bench/dataset_coverage.py`** — still references the
      pre-refactor `collect.py` shape (broke during the package
      reorganization). Noted in `02_biobtree_mining.md`.
- [ ] **Tune `body_gate.verdict()` threshold** — discovered when biobtree's
      RefSeq REVIEWED-only refresh tripped a `regression` verdict for what
      was actually a quality improvement (TP53 mRNA 46→25). The current
      "count dropped >50% → regression" rule can't distinguish "fewer
      predicted XM_/XR_ accessions" from "data unexpectedly missing".
      Consider: a `schema_changed` signal (when the response schema field
      set changes) demotes count drops to `drift`; or a per-key "expected
      direction" hint in section metadata.

## Open biobtree dependencies (waiting upstream)

Tracked in `docs/BIOBTREE_ISSUES.md`. Affecting current Atlas behavior:

| Issue | Status | Atlas impact |
|---|---|---|
| #4 silent multi-hop failure | Open | Largely subsumed by #1/#3 RESOLVED |
| #6 entry xrefs counts-not-values | Open | Not blocking — Atlas uses `map_all` to enumerate IDs |
| #9 UniProt CC narrative + reviewed flag + isoforms | Open — **biggest** | Blocks Path C #1 (UniProt CC) + Named isoforms |
| #10 AlphaFold empty for proteins >2700 aa | 🕓 Fix tomorrow (dev) | §4 currently constructs `AF-<acc>-F1` heuristically; pLDDT missing for ATM/BRCA2/DMD until fix |
| #12 pubchem_activity KRAS coverage gap | 🕓 Fix tomorrow (dev) | Workaround in place (chembl_activity covers KRAS); no functional gap on the page |

Resolved upstream 2026-05-30: #1 (not_found signaling), #2 (bare-400 on bad dataset),
#3 (statless empty blob), #5 (map source-dataset clarity), #7 (species filter
`[genome=="homo_sapiens"]`), #8 (pagination — earlier), #11 (ufeature ortholog
leakage). Atlas's `id.startswith(u + "_")` workaround removed (`72eae7e`).

## Pre-launch / cross-repo work (defer to launch day)

These all live in **`biobtree-content`** (the Hugo site) rather than this
repo, so they're deferred until V1 content is done — sending AI bots toward
the site too early caches an incomplete corpus. Drafts kept in
[`docs/site-drafts/`](../site-drafts/) and updated as Atlas evolves; copy
into `biobtree-content/static/` + tweak Hugo config the day of launch.

- [ ] Copy `docs/site-drafts/llms.txt` → `biobtree-content/static/llms.txt`
- [ ] Copy `docs/site-drafts/robots.txt` → `biobtree-content/static/robots.txt`
- [ ] Verify with `curl -I https://sugi.bio/{llms,robots}.txt`
- [ ] Hugo `config.toml`: enable `<lastmod>` per page in `sitemap.xml`
      (drives Perplexity freshness signal — top finding in
      `01_landscape_and_ai.md`).
- [ ] Hugo theme: emit `<link rel="alternate" type="application/ld+json"
      href="entity.jsonld">` and `<link rel="alternate" type="application/ld+json"
      href="provenance.json">` in the page head so machine clients discover
      the sidecars without scraping the body.
- [ ] Hugo serves `Last-Modified` HTTP header from each page's
      `generated_at` frontmatter (Hugo theme config; nginx pass-through).
- [ ] After deploy: resubmit sitemap to Google Search Console; watch
      `chatgpt_report.py` for the new AI-fetch patterns.

## Important — needs design decision (currently paused)

### UniProt CC narratives — how do we source them?

Audit's #1 gap (in `03_page_audit.md`). The curated UniProt narrative
(FUNCTION / SUBUNIT / SUBCELLULAR LOCATION / TISSUE SPECIFICITY / DISEASE /
PTM) is what makes UniProt + NCBI the AI default today. Shipping it as the
top-of-page block would flip that default for Atlas.

**biobtree's uniprot `entry` does NOT expose CC text** — only
`names / alternative_names / sequence / id / name`. Confirmed via curl on
2026-05-30 (after the upstream refresh — CC still absent). Filed as the v2
expansion of BIOBTREE_ISSUES #9.

Three implementation paths considered:

1. **Forward to biobtree** — already filed (BIOBTREE_ISSUES #9 v2). Slowest
   but most disciplined — single source of truth, every downstream user benefits.
2. **Direct UniProt REST per accession** (fetch `rest.uniprot.org/uniprotkb/
   {acc}.txt`, disk-cache). Fast to ship a parser; **rejected** because at
   full-corpus scale the per-gene network load + reliability liability in
   the hot path isn't worth it.
3. **Bulk UniProt parse ourselves** — download Swiss-Prot flat-file dump
   (~80MB compressed), parse CC blocks once, store the per-accession
   parsed result locally. Loads in seconds at pipeline time, no network per
   gene, rebuilds when UniProt releases a new version.

**Lean:** likely (3) for short-term, with (1) as the long-term path once
biobtree ships the CC fields. Existing parser code: biobtree itself parses
UniProt during ingest — the parser pattern is in `/data/biobtree` (Go); a
small Python equivalent would suffice. Could also reuse Biopython's
`Bio.SwissProt` parser.

**Decide before resuming.** Once resolved, this becomes the next Path C
item because of the audit's verdict on impact.

## Out of scope for now

- Open Targets-style genetic associations (L2G scores) — would require their
  GraphQL API. Worth revisiting after Path A+B+C land.
- Multi-language pages — defer until corpus is committed.
- Drug / disease entities — same shape as gene; do after gene pipeline is
  battle-tested at 100+ genes.

## Sample release-ready page checklist (for V1)

Before publishing a gene at sugi.bio/atlas/gene/X/:

- [x] page.md present, has all 12 sections
- [x] frontmatter has title, symbol, generated_at, atlas_version, biobtree_version
- [x] declarative first sentence present
- [x] "Updated: ..." visible
- [x] JSON-LD with `sameAs` links emitted (inline + sidecar)
- [x] body_gate verdict is `clean` or `first_run`
- [x] bundle.json + provenance.json + entity.jsonld sidecars written
- [ ] discoverable in sitemap with `<lastmod>` — Hugo config (pre-launch)
- [ ] `<link rel="alternate" ...>` discovery hints in page head — Hugo theme (pre-launch)
- [ ] `Last-Modified` HTTP header from `generated_at` — Hugo theme (pre-launch)

All 5 reference genes (TP53, BRCA1, CDKN2A, KRAS, TTN) currently satisfy the
in-repo subset of this checklist — see `/data/sugi-atlas-dist/atlas/gene/`.
