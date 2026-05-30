# Sugi Atlas — Landscape, AI Friendliness, and Differentiation

Research conducted 2026-05-30 to inform the design of `sugi.bio/atlas/gene/<SYMBOL>/` pages so they are preferentially fetched and cited by ChatGPT-User, Claude-User, Perplexity-User, Google AI Overview, and Bing Copilot. Grounded against the existing TP53 sample at `/data/sugi-atlas-dist/atlas/gene/TP53/page.md`.

## 1. Competitive landscape — per-gene reference pages

| Source | URL pattern | Unique data | Notable omissions / weak spots | Structured-data exposure | What AI scrapers actually pull |
| --- | --- | --- | --- | --- | --- |
| NCBI Gene | `ncbi.nlm.nih.gov/gene/<entrez>` (e.g. `/gene/7157`) | Canonical Entrez record, RefSeq linkage, GenBank cross-refs, GeneRIFs, HIV interactions table, multi-assembly coordinates | No schema.org / JSON-LD in `<head>`; meta description thin; sections collapse heavily, expression buried | TSV / JSONL / FASTA downloads; E-utilities API; **no** JSON-LD or schema.org | Title + first paragraph summary; rarely cited verbatim because page is dense and JS-heavy |
| UniProt | `uniprot.org/uniprotkb/<acc>/entry` | Curated function paragraphs, 1500+ sequence features per entry, PTMs, isoforms, evidence codes | Anti-bot fallback served to many fetchers (confirmed in WebFetch test); gene-symbol routing requires redirect | **RDF/Turtle export, SPARQL endpoint (232B triples, 2026_01), JSON, XML, GFF**; no schema.org JSON-LD | Mostly cited via JSON/RDF endpoints, not the HTML page; AI agents prefer `rest.uniprot.org` |
| GeneCards | `genecards.org/cgi-bin/carddisp.pl?gene=<SYMBOL>` | Broadest cross-DB aggregation, paralog/disorder panels, "GIFtS" score | 403/paywall to crawlers (confirmed); commercial; section content lifted from upstream sources | None public; commercial bulk licensing only | Often surfaced in result lists but blocked from extraction — AIs cite cached snippets |
| Open Targets | `platform.opentargets.org/target/<ENSG>` | Target-disease association scores, tractability, safety liabilities | Heavy SPA — initial HTML nearly empty for crawlers; data lives behind GraphQL | **GraphQL API, FTP bulk parquet/JSON**; no JSON-LD on the page | AIs almost never cite the HTML; they hit the GraphQL/FTP bundle |
| Open Targets Genetics | `genetics.opentargets.org/gene/<ENSG>` | L2G locus-to-gene scores, GWAS colocalisation | Same SPA opacity as Platform; merging with Platform UI | GraphQL only | Rarely cited directly |
| Ensembl | `ensembl.org/Homo_sapiens/Gene/Summary?g=<ENSG>` | Canonical transcript/exon coordinates, comparative genomics, regulatory build | 302 region-redirects break fetchers (confirmed); page sprawls across tabs | REST API, GFF3, BioMart; no JSON-LD on HTML | AIs fetch the REST API JSON; HTML page is a visual portal |
| MGI | `informatics.jax.org/marker/<MGI:ID>` | Mouse alleles, knockout phenotypes (381 mutations for Trp53), IMSR strains | Mouse-only; sparse human cross-walks; minimal in-page structured data | TSV reports, batch query; no JSON-LD | Cited for mouse phenotype questions when the page is reachable |
| MyGene.info | `mygene.info/v3/gene/<entrez>` | Real-time merged JSON from 30+ sources, lightweight | **API-only** — no HTML reference page exists | Pure JSON REST | High agent uptake as a data feed, but no "page" to cite |
| Human Protein Atlas | `proteinatlas.org/<ENSG>-<SYMBOL>` | IHC images across normal + cancer tissues, subcellular localisation, RNA tissue specificity scores | Image-centric — text content modest; per-tissue narrative is templated | XML/TSV download per protein; no JSON-LD | Cited for tissue/subcellular questions; the textual `<title>` and summary paragraph are the actual citation anchor |
| MalaCards | `malacards.org/card/<disease>` | Disease-centric aggregator (not strictly gene-centric); useful gene→disease panels | Pay-walled for bulk; commercial; many panels lifted upstream | None public | Rarely cited |
| DisGeNET | `disgenet.com/...` (migrated from `.org` in 2024; 301 confirmed) | Gene-disease association scores, variant-disease evidence, score breakdowns | Authentication wall on many endpoints since the commercial split; old URLs broken | REST API with token, RDF dump (Bio2RDF mirror) | Score-tables are cited when reachable; the `.com` migration broke a lot of AI memory |
| AlphaFold DB | `alphafold.ebi.ac.uk/entry/<UniProt>` | Predicted structures + per-residue pLDDT for every UniProt accession | Page is image-first; very little extractable text body | mmCIF / PDB / pLDDT JSON download; structured metadata API | AIs cite the existence of a model and the pLDDT band; rarely the page prose |
| GTEx Portal | `gtexportal.org/home/gene/<SYMBOL>` | Per-tissue eQTLs, sQTLs, isoform-level expression | SPA — initial HTML nearly empty for crawlers; tissue narrative is interactive plots | REST + bulk; no JSON-LD | Bulk files cited; page rarely |
| Pharos | `pharos.nih.gov/targets/<UniProt>` | TDL classification (Tclin / Tchem / Tbio / Tdark), drug-discovery readiness score, IDG illumination | Coverage skewed toward druggable target families | TCRD MySQL dump, REST API; no JSON-LD | TDL label cited; rest of page rarely |

### Cross-cutting observations
- **None of the major per-gene HTML pages publish schema.org `BioChemEntity` / `Gene` / `Protein` JSON-LD** in `<head>`. UniProt's RDF is the closest thing, but it is on a separate endpoint, not in the gene page. This is a wide-open hole.
- **SPA-rendered pages (Open Targets, GTEx, GeneCards, Ensembl tabs) systematically fail crawler extraction.** Static Markdown/HTML wins by default.
- **AI agents already prefer API-shaped citations** (MyGene.info, UniProt REST, Open Targets GraphQL) over HTML when both exist. A static page that *ships its own JSON sidecar* gets both audiences.
- **Anti-bot walls (GeneCards 403, UniProt fallback page, Ensembl region-redirect)** silently exclude these sources from many AI retrieval paths. Plain HTTP 200 + cache-friendly headers is itself a moat.

### Ranked landscape takeaways for Atlas
1. **Ship plain static HTML + Markdown with no JS gate** — already true; preserve it. This alone beats Open Targets / GTEx / GeneCards for crawler reachability.
2. **Publish a discoverable JSON-LD block in `<head>`** typed as `BioChemEntity` (subtype `Gene` where possible) with `associatedDisease`, `taxonomicRange`, `hasBioPolymerSequence`, `isEncodedByBioChemEntity`. None of the 13 competitors do this.
3. **Keep `bundle.json` as a sibling URL** (`/atlas/gene/TP53/bundle.json`) and link it from `<head>` with `<link rel="alternate" type="application/json">`. This mirrors how UniProt and MyGene.info win agent traffic.
4. **Stable, simple URL shape** (`/atlas/gene/<SYMBOL>/`) — outperforms Entrez-number URLs for symbol queries (LLMs almost always know the symbol, not the Entrez ID).
5. **No login, no cookies, no 403, no 302 region redirect** — a literal differentiator vs UniProt/Ensembl/GeneCards today.

## 2. AI-friendliness patterns — what gets cited

Synthesised from 2026 GEO/LLMO research and direct fetch tests above. Core findings:

- **Perplexity, Bing Copilot, Claude-User do real-time HTTP fetch on every query.** ChatGPT-User and Google AI Overview blend retrieval with cached training data — so two paths matter: (a) be in training data, (b) be extractable live.
- **Source-selection signals (live retrieval) cluster on:** structural extractability (clear H2s, definition-first paragraphs, tables), freshness (updated-in-last-N-days flag), authority signals (citations to primary sources, ORCID-style author IDs, methods transparency), and reachability (HTTP 200, sub-second TTFB, no JS gate).
- **Answer-absorption signals (what actually ends up in the synthesised reply):** quotable one-sentence definitions, labeled bullet facts, numeric values with units, and proximity of fact to evidence link.
- **llms.txt** has ~10% adoption, publicly supported by Anthropic and Perplexity, ignored or unused by Google. **Cheap to ship, modest upside, no downside.**
- **Schema.org FAQ/HowTo** has been over-recommended; recent 2026 analyses show *no* correlation between heavy FAQ schema and Google AI Overview citation rate. Don't waste effort on FAQ markup for a reference page — use `BioChemEntity` instead.
- **Healthcare/YMYL** content is held to a higher bar by Google AI Overview; in healthcare, ~24% of AI Overview citations overlap with the organic top-10 (higher coupling than other verticals). Atlas must look like an authoritative reference, not a blog.

### Concrete tactics Atlas should adopt (ranked)
1. **Embed schema.org JSON-LD `Gene` (or `BioChemEntity`) block in every page `<head>`** with `name`, `alternateName` (aliases), `identifier` (multiple, typed as `PropertyValue` with `propertyID: "HGNC"`, etc.), `taxonomicRange: "Homo sapiens"`, `associatedDisease` (one per OMIM/MONDO), `isInvolvedInBiologicalProcess` (GO BP), `bioChemInteraction` (top STRING partners), `sameAs` linking to NCBI, UniProt, Ensembl. **This is the single highest-leverage move.**
2. **Add `<link rel="alternate" type="application/json" href="bundle.json">` and `<link rel="alternate" type="application/ld+json" href="entity.jsonld">`** so agents can pivot from prose to structured form. Generate a clean `entity.jsonld` sidecar alongside `bundle.json`.
3. **Definition-first prose**: open every page with a single declarative sentence ("TP53 (HGNC:11998) encodes the tumour suppressor p53, a 393-aa transcription factor on 17p13.1.") — this is the most quoted line in 2026 GEO studies.
4. **Per-section H2s with self-contained answer tables.** Atlas already does this — keep it. Avoid wrapping facts in long paragraphs.
5. **`Last-Modified`, `ETag`, and visible "Updated: YYYY-MM-DD" in body.** Perplexity cited 30-day-fresh content 82% more often in 2026 testing.
6. **Ship `/llms.txt` at site root** linking to a curated index of gene pages and the data dictionary. Low cost, supported by Anthropic + Perplexity.
7. **Honest, narrow `robots.txt` posture: explicitly *allow* `ChatGPT-User`, `Claude-User`, `PerplexityBot`, `Google-Extended`, `Bingbot`, `Applebot-Extended`.** Many sites accidentally block these via overly broad `Disallow: /`.
8. **Add a `<meta name="citation_*">` (Highwire Press) block** in `<head>` — used by Google Scholar and increasingly by AI retrievers for academic-style attribution. Fields: `citation_title`, `citation_author`, `citation_publication_date`, `citation_doi` (if/when minted).
9. **Author + methods transparency**: a `provenance.md` linked from every page describing "every fact in this page was generated by `<biobtree-version>` from sources `<list>`; AI summary by Qwen3-235B with strict prompt `<hash>`." E-E-A-T loves this; current competitors do not show methods.
10. **Sitemap.xml with `<lastmod>` per gene page** and a separate `news-sitemap` style index for recently-regenerated pages. Helps both Bingbot and the Google freshness pipeline.

## 3. Strategic differentiation — where Atlas wins

Atlas's seven levers, ranked by AI-attraction impact:

| Rank | Lever | Why it wins citations | Counter-source weakness |
| --- | --- | --- | --- |
| 1 | **Deterministic provenance — every fact traceable to biobtree → upstream** | E-E-A-T trust signal; LLMs increasingly weight reproducible sourcing; nothing else in the landscape ships a fact-by-fact provenance manifest | GeneCards/Open Targets are opaque aggregators; AIs are wary of citing them |
| 2 | **`bundle.json` sidecar (per-page machine-readable mirror)** | Agents fetch JSON 5-10× more often than HTML when both exist; no competitor ships a clean per-page JSON | UniProt RDF is whole-database, not per-page; MyGene has no HTML twin |
| 3 | **Unioned GO/Reactome across all reviewed UniProt isoforms** | Dual-product / multi-isoform genes (CDKN2A, INS-IGF2, MLL) are systematically mishandled elsewhere; Atlas gives the union and labels it | Open Targets and GeneCards typically pick the canonical isoform and silently drop the rest |
| 4 | **Multi-hop disease routes (GWAS→MONDO→gene when direct mapping is empty)** | Genes like AR have empty direct GWAS-gene tables in Open Targets; Atlas fills these via mondo bridging | Competitors return "no results" — Atlas returns the indirect route with the chain made explicit |
| 5 | **UniProt sequence-features at scale (1518 features for TP53 with location ranges + mechanism descriptions)** | Highly quotable for "what does the R175H mutation do" style questions; surface these as their own section per residue range | UniProt has them but they are buried in a JS-driven feature viewer that crawlers cannot render |
| 6 | **AlphaFold pLDDT alongside experimental PDB** | "Is the C-terminus disordered" is a common LLM question; Atlas's pLDDT bands answer it inline | AlphaFold DB shows it visually; PDB does not show pLDDT at all |
| 7 | **Per-tissue Bgee scores (numeric, comparable)** | Quotable "TP53 expression is X in liver, Y in testis" answers; GTEx has the data but in a JS plot, HPA has narrative not numbers | GTEx requires API; HPA is qualitative ("low specificity") |

### Three to five new unique-value features Atlas should add (concrete, none of the 13 competitors do them)
1. **`entity.jsonld` sidecar** — schema.org `Gene` JSON-LD per page with `sameAs` linking out to all 12+ external IDs and `citation` linking inward to the section anchors. This makes Atlas the first per-gene resource with a clean public knowledge-graph node.
2. **Per-fact provenance hover-anchors** — every numeric or categorical claim gets an HTML anchor (`<span id="f-clinvar-pathogenic-count">412</span>`) plus an entry in `provenance.json` keying that anchor to `{biobtree_call, upstream_url, fetched_at}`. Lets an agent quote a fact *and* link to its origin in one shot. No competitor offers this.
3. **`/atlas/gene/<SYMBOL>/diff/<prev-version>.md`** — published diff between regenerations of the page. Freshness *and* auditability. Wikipedia-style trust signal nothing in bioinformatics offers.
4. **`/atlas/gene/<SYMBOL>/ask.md`** — pre-rendered, deterministic answers to the top ~30 questions ("What does TP53 do?", "What are common TP53 mutations?", "What drugs target TP53?") with the answer above and the evidence row below. Direct match to the way Perplexity/Google AI assemble responses. Treat it as a static FAQ but typed as `BioChemEntity`/`MedicalCondition`, not `FAQPage` (FAQPage schema is now devalued).
5. **`/atlas/llms.txt` and `/atlas/llms-full.txt`** — root index of every gene page with one-line definitions, machine-curated. Cheap, ~10% of the web has it, and Anthropic + Perplexity actively retrieve it.

### Top 5 actions, ranked
1. Add schema.org `Gene` JSON-LD in `<head>` of every page, with `sameAs` to NCBI/UniProt/Ensembl/HGNC/OMIM/MONDO.
2. Publish `entity.jsonld` and keep `bundle.json` as siblings; link both via `<link rel="alternate">`.
3. Open every page with a one-sentence canonical definition + visible "Updated: <date>" + `Last-Modified` header.
4. Ship `/llms.txt`, a sitemap with `<lastmod>`, and a `robots.txt` that explicitly allows the five major AI user-agents.
5. Add `provenance.json` + per-fact HTML anchors so agents can cite individual numbers, not just the page.

---
*Sources: WebFetch tests against ncbi.nlm.nih.gov, uniprot.org, genecards.org, platform.opentargets.org, proteinatlas.org, informatics.jax.org, pharos.nih.gov, mygene.info, ensembl.org, disgenet.org (May 2026); WebSearch on schema.org BioChemEntity, llms.txt adoption, Perplexity/ChatGPT/AI Overview citation studies (2026).*
