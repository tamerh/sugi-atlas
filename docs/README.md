# Atlas documentation

How Atlas builds deterministic gene, drug, and disease reference pages from a
local [biobtree](https://biobtree.org) instance.

**Start here:**

- **[how-it-works.md](how-it-works.md)** — the shared architecture: the
  collect→render philosophy, the biobtree transport, the frozen page contract,
  JSON-LD/provenance, corpus building, and the validation gates. Read this first.

**How each page is built** (the mining, section by section):

- **[gene.md](gene.md)** — gene pages: anchors, the eight zones, dual-product
  genes, the two-tier data rule.
- **[drug.md](drug.md)** — drug pages: ChEMBL resolution, the GtoPdb-only curated
  target discipline.
- **[disease.md](disease.md)** — disease pages: the Mondo corpus, the gene-cohort
  fan-out, title-validated trials.

**Cross-cutting:**

- **[mesh.md](mesh.md)** — the cross-entity graph: how a fact on one page becomes
  a navigable link on another.
- **[PAGE_CONTRACT.md](PAGE_CONTRACT.md)** — the frozen H2/anchor reference spec
  (treated as a stable public API; the integration tests encode it).
