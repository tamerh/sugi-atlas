# Sugi Atlas — corpus statistics & validation (preprint draft)

> **Working/tmp doc** — not committed. Drafted 2026-06-04 for the preprint's
> "corpus statistics + validation" table. Numbers from the completed prod build
> (`853a8aa`, full corpus) and the `c48cb41` archive (section coverage);
> validation counts from the test suite at `7067e31`. All deterministic and
> reproducible from the commit. **Re-measure after the biobtree-release run** —
> coverage (esp. CIViC/trials) will shift up once #32 under-linking is fixed.

---

## Table 1 — Corpus scale

| Entity type | Pages | Share |
|---|---:|---:|
| Genes | 29,338 | 56.4% |
| Diseases | 19,683 | 37.8% |
| Drugs | 2,987 | 5.7% |
| **Total** | **52,008** | — |

**Seeds → pages → skipped: 52,008 → 52,008 → 0.** Full-corpus seeds are
pre-resolved HGNC / Mondo / ChEMBL identifiers, so every seed yields a page.
(Resolution skips occur only in the name-seeded pre-production check — e.g.
Orphanet-only synonyms with no Mondo id — never in the published corpus.)

Build time: full corpus ~43 min (collect+merge+render), 30 workers.

---

## Table 2 — Section-level data coverage

% of pages carrying real data in each evidence-bearing section (vs. an explicit
"no data" placeholder). Computed deterministically from the build (`c48cb41`
archive: 29,338 genes / 19,681 diseases / 2,987 drugs).

### Genes (n = 29,338)
| Section (signal) | Pages w/ data | % |
|---|---:|---:|
| Variants (ClinVar) | 18,863 | 64.3 |
| Protein interactions | 18,526 | 63.1 |
| GWAS associations | 15,175 | 51.7 |
| Clinical trials | 7,204 | 24.6 |
| Drugs | 2,543 | 8.7 |
| CIViC evidence | 333 | 1.1 |

### Diseases (n = 19,681)
| Section (signal) | Pages w/ data | % |
|---|---:|---:|
| Cohort genes | 10,278 | 52.2 |
| Variants | 9,334 | 47.4 |
| Clinical trials | 6,067 | 30.8 |
| Drugs / targets | 3,994 | 20.3 |
| GWAS | 1,165 | 5.9 |
| CIViC evidence | 422 | 2.1 |

### Drugs (n = 2,987)
| Section (signal) | Pages w/ data | % |
|---|---:|---:|
| Indications | 2,536 | 84.9 |
| Clinical trials | 2,268 | 75.9 |
| Patent families | 2,247 | 75.2 |
| Targets | 1,840 | 61.6 |
| Bioactivity | 1,840 | 61.6 |
| CIViC evidence | 177 | 5.9 |

**Framing note:** low absolute rates (CIViC 1–6%, gene–drug 8.7%) are not Atlas
gaps — they reflect how sparse *curated* evidence genuinely is across the whole
genome/disease space. The per-section placeholder makes "no curated data"
explicit rather than hiding it. State this as a design property.

---

## Table 3 — Validation suite (100% pass over all 52,008 pages)

| Suite | Tests | Verifies |
|---|---:|---|
| Unit | 168 | renderers, evidence/percentile logic, declarative leads, JSON-LD builders, helper precision (causal-gene, alias splitting, ChEBI roles) |
| Integration · contract | 6 | frozen H2/H3 anchor set & order, per entity type |
| Integration · frontmatter | 10 | required fields, typed identifiers, `evidence_score`, `section_defaults` |
| Integration · JSON-LD | 5 | schema.org Gene / MedicalCondition / Drug validity |
| Integration · cross-links | 7 | reverse index, no dangling internal links |
| Integration · data quality | 15 | nan/None/boolean leaks, column-shift, duplicate rows, HTML-entity, chemistry fragments, self-edges, float32 artifacts |
| Integration · coverage | 1 | manifest ↔ rendered pages, 1:1 |
| **Integration total** | **44** | run over the **full** corpus, not a sample |

Each integration lint maps to a specific data-quality defect that was found and
fixed (test_quality.py header: "each maps to a fix we shipped").

---

## To do when we return
- Re-measure Tables 1–2 after the biobtree-release run (new IDs in `<HEAD>`).
- Decide: add `atlas.sh stats` to regenerate Tables 1–2 from any corpus so the
  paper's numbers stay reproducible (script exists in spirit — the snippet that
  produced Table 2 streams `evidence_components` from the archive).
- Optional: add the corpus-construction funnel (Mondo universe → admission-gated
  disease seeds) if the reviewer wants the upstream selection numbers too.
