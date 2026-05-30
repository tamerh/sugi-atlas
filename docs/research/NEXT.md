# Atlas roadmap — derived from research reports

Synthesized from [`01_landscape_and_ai.md`](01_landscape_and_ai.md),
[`02_biobtree_mining.md`](02_biobtree_mining.md),
[`03_page_audit.md`](03_page_audit.md).

Update items in place as they ship: leave the box checked, append the commit
SHA + date. Add new items only when grounded in a research finding or
field-discovered need.

---

## Path A — AI-friendliness envelope (do first)

Cheap. Decisive lead vs everyone. Without these, content depth goes
undiscovered by AI crawlers.

- [ ] **schema.org `Gene`/`BioChemEntity` JSON-LD** in every page `<head>`,
      with `sameAs` to NCBI Gene / UniProt / Ensembl / HGNC / OMIM. Emit via
      `src/atlas/page/jsonld.py`; lock in with `tests/unit/test_jsonld.py`.
- [ ] **Declarative first sentence** on every page (one well-formed standalone
      sentence stating what the gene is). Generate from anchors (symbol + name
      + location + function class) deterministically; the LLM exec summary
      stays as commentary below.
- [ ] **Visible "Updated: YYYY-MM-DD"** at the top of every page + `Last-Modified`
      HTTP header (Hugo can do this from page frontmatter).
- [ ] **`/llms.txt`** at the Hugo site root listing site purpose, content scope,
      machine-readable endpoints, citation guidance.
- [ ] **`robots.txt` allowlist** for ChatGPT-User, Claude-User, PerplexityBot,
      Google-Extended, Bingbot, OAI-SearchBot.
- [ ] **Per-gene `<lastmod>` sitemap**; Hugo's sitemap supports this — verify
      it's enabled.

## Path B — Provenance moat (do second)

Establishes Atlas's unique position: the only gene resource where an agent can
cite a single fact with its upstream source.

- [ ] **`provenance.json` sidecar** in `dist/atlas/gene/<SYM>/`. For each
      numeric fact in the page, record: bundle key, biobtree call (URL +
      params), upstream source (NCBI / UniProt / etc.) and its URL.
- [ ] **Per-fact HTML anchors** on the page. Each numeric claim gets an `id`
      so AI agents can deep-link a single fact.
- [ ] **`<link rel="alternate" type="application/ld+json" href="entity.jsonld">`**
      and similar for `bundle.json` / `provenance.json` so machine clients
      discover the sidecars without scraping.

## Path C — Content depth (do third, in parallel where possible)

From the biobtree mining + page audit. Each is independently shippable.

### From the page audit (biggest user-visible gaps)
- [ ] **UniProt CC narratives** (FUNCTION / SUBUNIT / SUBCELLULAR LOCATION /
      TISSUE SPECIFICITY / DISEASE / PTM) at top of §3. Verify biobtree exposes
      these — if not, fetch UniProt's `.txt` directly with caching.
- [ ] **Named isoforms** (p53α/β/γ, K-Ras4A/4B, p16γ, titin N2A/N2B/N2BA/novex).
      Source from UniProt `.txt` ALTERNATIVE PRODUCTS section.
- [ ] **Drug → indication → biomarker triples** in §10 (sotorasib + KRAS-G12C +
      NSCLC; osimertinib + EGFR-T790M; olaparib + BRCAness). Sources to bridge:
      ChEMBL indication + FDA biomarker list.

### From the biobtree mining (curl-verified, ready to wire)
- [ ] **mirdb** via `>>hgnc>>refseq>>mirdb` — 77–100 miRNAs/gene, §11.
- [ ] **chembl_molecule drilldowns** per phased drug (chembl_activity,
      clinical_trials, pubchem, chebi) — §10.
- [ ] **pubchem_activity** via `>>uniprot>>pubchem_activity` — Ki/IC50 +
      qualifier + unit + PMID per row, §10.
- [ ] **ctd_gene_interaction** via `>>hgnc>>entrez>>ctd_gene_interaction` —
      literature-mined chem-gene with CV verbs + PubMed counts, §10.
- [ ] **MeSH disease descriptors** via `>>hgnc>>clinvar>>mondo>>mesh` —
      35 rows/gene with MeSH tree numbers + scope_notes, §12.
- [ ] **Anatomy/cell-type names** via `bgee>>uberon` and `bgee>>cl` to replace
      bare UBERON/CL ids in §11 — trivial readability win.

## Hygiene fixes

- [ ] **Fix `src/atlas/bench/dataset_coverage.py`** — still references the old
      `collect.py` shape, broke during the refactor.

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

## Out of scope for now

- Open Targets-style genetic associations (L2G scores) — would require their
  GraphQL API. Worth revisiting after Path A+B+C land.
- Multi-language pages — defer until corpus is committed.
- Drug / disease entities — same shape as gene; do after gene pipeline is
  battle-tested at 100+ genes.

## Sample release-ready page checklist (for V1)

Before publishing a gene at sugi.bio/atlas/gene/X/:

- [ ] page.md present, has all 12 sections
- [ ] frontmatter has title, symbol, generated_at, atlas_version, biobtree_version
- [ ] declarative first sentence present
- [ ] "Updated: ..." visible
- [ ] JSON-LD with `sameAs` links emitted
- [ ] body_gate verdict is `clean` or `first_run`
- [ ] bundle.json + provenance.json + entity.jsonld sidecars written
- [ ] discoverable in sitemap with `<lastmod>`
