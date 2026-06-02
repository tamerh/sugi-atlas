# The cross-entity mesh

The mesh is what turns three independently-built corpora — genes, drugs,
diseases — into one navigable graph. An assertion collected on one page ("TP53's
CIViC evidence names Venetoclax") becomes a link on *both* pages: a forward edge
on TP53, a reverse edge on Venetoclax. It is built entirely from data the
collectors already gathered (no extra network calls), through two build-time JSON
sidecars. For the page machinery the mesh plugs into, see
[how-it-works.md](how-it-works.md).

Two principles drive the design:

- **A page only links to a target that has actually been built.** Until a target
  is in the manifest, the resolver returns nothing and the caller renders plain
  text. The mesh grows monotonically — a second full build completes edges the
  first couldn't yet resolve.
- **Only curated edges are meshed.** Bioactivity assay hits and positional
  ontology overlaps are deliberately excluded. A biomarker is not a target; an
  off-target ChEMBL hit is not a target. This is what makes the labels mean what
  they say.

## The two sidecars

**`manifest.json` — forward resolution.** Maps every *resolvable key* to a slug,
per entity, plus a `canon` slug→display-name map:

| Entity | Keys that resolve to the page |
|---|---|
| gene | symbol, HGNC id |
| disease | MONDO id, EFO id, canonical name + synonyms |
| drug | ChEMBL id, parent ChEMBL, child ChEMBLs, canonical name + alt-names |

IDs are stored verbatim; names are normalised (lowercased, non-alphanumerics
collapsed). Storing a drug's parent *and* child ChEMBLs means a salt-form id and
the parent both resolve to the one canonical slug. The `canon` map exists because
a *synonym* key can resolve to a page whose canonical name differs ("schizoaffective
disorder" → `/schizophrenia/`), so a link can render the name of the page it
actually opens. Lookup tries each key verbatim first, then normalised.

**`reverse_edges.json` — reverse index.** Maps each target URL to the list of
sources that point at it. It is literally the *inversion* of every resolved
forward edge, built in the merge phase after the complete manifest exists (so
every edge resolves). For single-page builds the file is absent and the reverse
mesh degrades silently.

## One source of truth

A single function, `related_targets(entity, bundle)`, computes the forward edges
for a page, returning grouped, deduped Genes / Diseases / Drugs. It is used by
**all three** consumers — the human-facing markdown block, the reverse-edge index,
and the JSON-LD `@reverse` edges — so forward and reverse can never disagree.
Only curated sources feed it:

- **gene** → diseases from GenCC + ClinGen validity; drugs + diseases from CIViC
  predictive evidence. *Not* raw OMIM∩page noise, *not* off-target bioactivity.
- **disease** → cohort genes (ranked by evidence strength, CIViC drivers first);
  drugs from title-validated trial drugs + CIViC therapies. *Not* the bioactivity
  drug cloud.
- **drug** → genes from GtoPdb-curated primary targets only; diseases from
  labelled indications + CIViC.

## The four reverse directions

The key insight: a forward edge from a source of a given type becomes, *on the
target page*, a group labelled by the **predicate viewed from the other end**.
Only directions whose source edge is curated are inverted:

| Forward edge | Reverse label on the target page | Meaning |
|---|---|---|
| gene → drug (CIViC) | **Biomarker genes** | genes whose variants associate this drug — shown on drug pages |
| drug → gene (GtoPdb) | **Targeted by drugs** | drugs that curatedly target this gene — shown on gene pages |
| gene → disease | **Associated genes** | genes asserting association — shown on disease pages |
| drug → disease | **Drugs indicated** | drugs whose indication is this disease — shown on disease pages |

Note the asymmetry: **"Biomarker genes" ≠ "Targeted by drugs"** is the encoded
discipline — being a biomarker for a drug is not the same as being its target.

## The Related block

The `## Related Atlas pages {#related}` section lists forward groups first
(Genes → Diseases → Drugs), then the reverse groups in fixed order. Its rules:

- **Self-links excluded** — a drug's own salt forms or a self-naming CIViC therapy
  never appear; a URL shows under exactly one label.
- **De-SHOUTed drug labels** — CISPLATIN renders as Cisplatin, matching the drug
  page title.
- **Honest disease labelling** — on a disease page the gene group is relabelled
  "Cohort genes", so a polygenic cohort isn't read as causal genes.
- **Capped** at 12 per row with a `(+N more)` suffix.
- **Always emitted** — if nothing is built, a placeholder line appears rather than
  dropping the section, so the table of contents is identical across every page.

## Bidirectionality, worked

TP53's CIViC evidence names Venetoclax → a forward `Drugs` edge on the TP53 page.
The merge phase inverts it, so Venetoclax's page shows TP53 under **Biomarker
genes**. The relationship asserted from one side is navigable from both — and the
JSON-LD carries the same edge for machine traversal. Integrity is enforced over
the built corpus by the integration tests: no dangling links, manifest⇄page
bijection, every reverse URL resolves, destination-canonical labels, and no
self-links.
