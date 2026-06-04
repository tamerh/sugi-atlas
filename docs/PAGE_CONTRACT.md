# Atlas page output contract — H2 taxonomy & anchor IDs

**Status: FROZEN 2026-06-02** (ratified by backend + web team). Treat the anchor
IDs as a **stable API** — same discipline as a public REST endpoint. Rename only
via deprecation + redirect, **never silently** (this is what breaks deep links
in the wild — the GeneRIF "(showing 40)→(showing 50)" class of bug).

## Rules
- Each entity emits its canonical H2 set **in the fixed order below** — never
  alphabetical or data-availability-driven (the frontend TOC follows source
  order).
- Heading IDs are explicit, **lower-kebab**, via goldmark attribute syntax:
  `## Label {#id}`. Requires Hugo `markup.goldmark.parser.attribute.title = true`
  (web team owns this config). IDs are decoupled from prose.
- **Every canonical section is always emitted** — heading + `{#id}` + body, or
  heading + `{#id}` + an *informative* placeholder (cite the expected source,
  e.g. `*No GWAS associations recorded — common-variant studies don't cover this
  disease.*`), never "No data" and never a silently-missing section. This keeps
  the TOC identical across every page of a type.
- `#summary` is always first; `#related` always last.
- Shared IDs mean the same concept across types → cross-entity nav
  ("Same section on EGFR / KRAS →"). Shared: `#summary #identifiers #function
  #disease #drugs #trials #related`. Labels may differ per type; **IDs stay shared**.
- Sub-tables/blocks within a section are H3 with their own stable kebab IDs.

## GENE (8 H2)
| Order | Label | `{#id}` | Absorbs (current renderers) |
|---|---|---|---|
| 1 | Summary | `#summary` | lead sentence + RefSeq + At-a-glance (one wrapper; lead stays the most-indexable line) |
| 2 | Identifiers | `#identifiers` | §1 gene identifiers |
| 3 | Gene structure | `#gene-structure` | §2 transcripts, §11 expression, §9 regulation, functional-genomics, GeneRIFs |
| 4 | Protein | `#protein` | §3 protein ids, §4 structure, residue map (dual-product → H3 `{#protein-<acc>}`, matches JSON-LD `@id`) |
| 5 | Function | `#function` | §7 pathways & GO, §8 interactions |
| 6 | Disease & clinical | `#disease` | **Cancer significance** (intOGen+CIViC, folded here — NOT in Summary), §6 clinical variants, §12 disease associations |
| 7 | Drugs & pharmacology | `#drugs` | §10 drugs (GtoPdb, BindingDB, ChEMBL, PharmGKB) |
| 8 | Related Atlas pages | `#related` | cross-entity mesh (forward + reverse) |

H3 IDs: `#transcripts #expression #regulation #generif #residue-map #pathways #interactions #variants #gtopdb`

## DISEASE (8 H2 — consolidates the current 18)
| Order | Label | `{#id}` | Absorbs |
|---|---|---|---|
| 1 | Summary | `#summary` | lead + At-a-glance |
| 2 | Clinical features | `#clinical` | Epidemiology (prevalence), Signs & symptoms (HPO) — the headline clinical presentation, lifted out of Identifiers |
| 3 | Identifiers | `#identifiers` | disease identifiers, synonyms, data availability |
| 4 | Genetics & variants | `#genetics` | GWAS landscape, Variant details & tiers |
| 5 | Genes & proteins | `#genes` | Mendelian/GenCC overlap & somatic drivers, Cohort genes→proteins, Protein-family classification, Interactions among cohort, Structural data, Expression context |
| 6 | Function | `#function` | Pathway analysis |
| 7 | Therapeutics | `#drugs` | Drug-target analysis, Bioactivity/enzyme, Pharmacogenomics, Chemical tractability, Druggability pyramid, Undrugged target profiles |
| 8 | Clinical trials & evidence | `#trials` | Clinical trials, CIViC |
| 9 | Related Atlas pages | `#related` | mesh |

H3 IDs: `#epidemiology #symptoms #gwas #variant-tiers #gencc #cohort-genes #tractability #civic`

## DRUG (7 H2)
| Order | Label | `{#id}` | Absorbs |
|---|---|---|---|
| 1 | Summary | `#summary` | lead + At-a-glance |
| 2 | Identifiers | `#identifiers` | §1 drug identity & classification |
| 3 | Targets | `#targets` | §2 targets, §3 bioactivity, target pathways |
| 4 | Indications & clinical | `#indications` | §4 indications, §5 clinical trials, CIViC clinical evidence |
| 5 | Pharmacology | `#pharmacology` | pharmacogenomics, pharmacology |
| 6 | Related molecules | `#related-molecules` | §7 related molecules |
| 7 | Related Atlas pages | `#related` | mesh |

H3 IDs: `#bioactivity #indications #trials #civic #pharmacogenomics`

> Drug uses `#indications` (not the shared `#disease`) — reads naturally in a
> shared URL and in Google "Jump to" sitelinks; cross-entity nav uses a
> frontend mapping table, not anchor identity.
