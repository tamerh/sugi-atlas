# Sugi Atlas — corpus statistics & validation (preprint draft)

> **Working/tmp doc** — not committed. Re-measured 2026-06-04 off the comprehensive
> build `f20d676` (admission gate v2 [non-human dropped], phase-3 drugs added,
> audit fixes, index.md layout). Archive: atlas-corpus-20260604-193822-f20d676.
> All deterministic, reproducible from the commit. **Re-measure once more after
> the biobtree-release run** — CIViC/trials coverage will rise when #32
> under-linking is fixed.

---

## Table 1 — Corpus scale

| Entity type | Pages | Share |
|---|---:|---:|
| Genes | 29,338 | 55.7% |
| Diseases | 18,618 | 35.4% |
| Drugs | 4,701 | 8.9% |
| **Total** | **52,657** | — |

**Seeds → pages → skipped: 52,657 → ~52,655 → ~2.** Full-corpus seeds are
pre-resolved HGNC / Mondo / ChEMBL identifiers, so essentially every seed yields
a page. **Admission gates** drop, at seed time: gate-1 disease-characteristic
qualifiers; **gate-2 non-human / veterinary disease terms (~1,065** — e.g.
achondroplasia-in-cattle). **Drug seed = approved (phase 4) ∪ late-stage clinical
(phase 3)**, de-salted + therapeutic-filtered (4,225+1,892 → 4,701).

Build time: full corpus ~44 min (collect+merge+render), 30 workers.

---

## Table 2 — Section-level data coverage

% of pages carrying real data in each evidence-bearing section (vs. an explicit
"no data" placeholder). Computed deterministically from the
build f20d676: 29,338 genes / 18,616 diseases / 4,701 drugs).

### Genes (n = 29,338)
| Section (signal) | % |
|---|---:|
| Variants (ClinVar) | 64 |
| Protein interactions | 63 |
| GWAS associations | 52 |
| Clinical trials | 25 |
| Drugs | 9 |
| CIViC evidence | 1 |

### Diseases (n = 18,616) — note coverage rose vs the pre-gate build (animal empties removed)
| Section (signal) | % |
|---|---:|
| Cohort genes | 55 |
| Variants | 50 |
| Clinical trials | 33 |
| Drugs / targets | 21 |
| GWAS | 6 |
| CIViC evidence | 2 |

### Drugs (n = 4,701) — approved + phase-3
| Section (signal) | % |
|---|---:|
| Indications | 87 |
| Clinical trials | 82 |
| Patent families | 69 |
| Targets | 53 |
| Bioactivity | 53 |
| CIViC evidence | 5 |

(Drug target/bioactivity % dipped vs approved-only — phase-3 candidates more
often lack a fully-resolved target; trials % rose since phase-3 drugs are
trial-heavy and the trials count is now the true total, not the 1,600 cap.)

**Framing note:** low absolute rates (CIViC 1–5%, gene–drug 9%) are not Atlas
gaps — they reflect how sparse *curated* evidence genuinely is across the whole
genome/disease space. The per-section placeholder makes "no curated data"
explicit rather than hiding it. State this as a design property.

---

## Table 3 — Validation suite (100% pass over all 52,657 pages)

| Suite | Tests | Verifies |
|---|---:|---|
| Unit | 170 | renderers, evidence/percentile logic, declarative leads, JSON-LD builders, helper precision (causal-gene, alias splitting, ChEBI roles) |
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
