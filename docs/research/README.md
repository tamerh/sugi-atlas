# Atlas research reports

Three parallel research passes (2026-05-30) on what would make Atlas's per-gene
pages attract AI citations and how Atlas compares to the existing gene-info
landscape. Each report stands alone; this index summarizes the top finding(s)
and points at the actionable backlog in [`NEXT.md`](NEXT.md).

> **Start here for the latest, deepest pass:**
> [`05_unique_ai_attraction.md`](05_unique_ai_attraction.md) (2026-06-01) — a
> fresh 4-agent pass focused on what makes Atlas *uniquely* attractive to AI/LLMs.
> It supersedes parts of `01`–`04` (notably it **corrects the "JSON-LD is the #1
> AI-citation lever" and "llms.txt is a priority" claims** with controlled-study
> evidence) and adds the headline strategy 01–04 missed: the **biobtree-MCP ↔
> Atlas grounding flywheel**, plus 12 live-verified computable content moats and a
> shipped-dist audit with 2 launch blockers. Supporting agent drafts live in
> [`_drafts/`](_drafts/).

## [`01_landscape_and_ai.md`](01_landscape_and_ai.md) — Competitive landscape & AI-friendliness
*Surveyed: NCBI Gene, UniProt, GeneCards, Open Targets (+ Open Targets Genetics),
Ensembl, MGI, MyGene.info, Human Protein Atlas, MalaCards, DisGeNET, AlphaFold
DB, GTEx Portal, Pharos.*

- **None of the 13 surveyed competitors ship schema.org `Gene`/`BioChemEntity`
  JSON-LD** — single largest free win for AI citation.
- **Definition-first paragraph + visible "Updated: YYYY-MM-DD"** drives ~82% of
  Perplexity citation lift in 2026 studies. FAQPage schema is now devalued.
- **AI-friendliness essentials** Atlas should adopt: `llms.txt`, allow-listed
  `robots.txt` for ChatGPT-User/Claude-User/PerplexityBot/Google-Extended/Bingbot,
  per-gene `<lastmod>` sitemap, `Last-Modified` headers, JSON-LD with `sameAs`
  to NCBI/UniProt/Ensembl/HGNC/OMIM.
- **The signature moat**: per-fact HTML anchors + `provenance.json` sidecar
  keying each numeric claim to the biobtree call + upstream source. No other
  gene resource exposes per-fact provenance to agents.

## [`02_biobtree_mining.md`](02_biobtree_mining.md) — Underutilized biobtree capabilities
*Empirical curl probes against TP53/BRCA1/EGFR; every recommendation backed by
sample output.*

- **117 edge-graph nodes total. Atlas uses 51, leaves 66 uncovered.** Of those
  66: ~17 are wire-worthy as integrated data, ~14 conditional/hierarchical, the
  rest xref-only or upstream-blocked.
- **Top 5 to wire** (value × effort): **mirdb** (miRNAs targeting the gene,
  scored — Atlas has zero miRNA today); **chembl_molecule drilldowns** (per
  phased drug → activities + trials + pubchem/chebi); **pubchem_activity** (Ki
  /IC50 + qualifier + unit + PMID); **ctd_gene_interaction** (literature-mined
  chem-gene with CV verbs); **MeSH** disease descriptors via `>>hgnc>>clinvar
  >>mondo>>mesh` (tree numbers + scope_notes).
- **Trivial readability wins**: `bgee>>uberon` and `bgee>>cl` give canonical
  anatomy/cell-type names for bare UBERON/CL ids already in §11.
- **Upstream blocked / don't wire**: scxa_expression (projection bug),
  gtopdb_*, pharmgkb_clinical/_guideline, cellphonedb — advertised edges that
  return `not_found` or empty projections.
- **Discovered:** `dataset_coverage.py` is broken since the refactor (still
  references the old `collect.py` shape). Small fix.

## [`03_page_audit.md`](03_page_audit.md) — Atlas vs the world, per-gene
*Side-by-side against NCBI Gene / UniProt / GeneCards / Open Targets / Ensembl
for TP53, BRCA1, EGFR, CDKN2A, KRAS, TTN.*

- **The honest verdict:** an AI agent answering "tell me about TP53" today
  picks UniProt/NCBI first because they lead with declarative narrative.
  Atlas wins decisively on every deep follow-up. Shipping UniProt CC blocks
  would flip the default.
- **Top 3 gaps to close:** UniProt curated free-text CC narratives (FUNCTION,
  SUBUNIT, SUBCELLULAR LOCATION, TISSUE SPECIFICITY, DISEASE, PTM); named
  isoforms (p53α/β/γ, K-Ras4A/4B, p16γ, N2A/N2B/N2BA/novex); drug → indication
  → biomarker triples (sotorasib + KRAS-G12C + NSCLC; osimertinib + EGFR
  T790M; olaparib + BRCAness).
- **Top 3 differentiators** to lean on: paired AlphaMissense + SpliceAI +
  ClinVar at full scale on one page; TF regulatory layer in-line; honest
  reporting (`Antibody resources: 0` when zero, counts-as-floors flags) and
  dual-product correctness on CDKN2A.

## Triangulated takeaway

All three reports converge: Atlas has unmatched **data depth + provenance**;
it loses on **opening narrative + structured-data envelope**. The path is
(1) ship the AI-friendly envelope, (2) add declarative-first paragraphs from
UniProt CC, (3) expose per-fact provenance, (4) wire the new biobtree datasets.
See [`NEXT.md`](NEXT.md) for the prioritized roadmap.
