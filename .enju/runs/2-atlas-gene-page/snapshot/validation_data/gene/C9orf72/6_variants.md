## Clinical variants & AI predictions

### Clinical Variants (ClinVar)

**Summary**
| Metric | Count |
|--------|-------|
| Total variants | 50 |
| Pathogenic | 3 |
| Likely Pathogenic | ~5 |
| Uncertain Significance | ~32 |
| Likely Benign | ~5 |
| Benign | ~5 |

**Top Pathogenic/Likely Pathogenic Variants**

| Variant ID | HGVS Notation | Condition | Classification |
|------------|---------------|-----------|-----------------|
| 31151 | NM_001256054.1:c.-45+163GGGGCC[>24] | Frontotemporal dementia and/or ALS-1 | Pathogenic |
| 183034 | NG_031977.1:g.(5321_5338)ins(60_?) | Frontotemporal dementia and/or ALS-1 | Pathogenic |
| 1343330 | NC_000009.12:g.27573529_27573534GGCCCC[60_?] | Amyotrophic lateral sclerosis | Pathogenic |
| 1192640 | NM_018325.5:c.600+27A>G | FTD/ALS-1 | Uncertain significance |
| 1192639 | NM_018325.5:c.665+115_665+117dup | FTD/ALS-1 | Uncertain significance |
| 366506 | NM_018325.5:c.1426G>C (p.Asp476His) | FTD/ALS-1 | Uncertain significance |
| 4355618 | NM_018325.5:c.505-1G>T | FTD/ALS-1 | Uncertain significance |

*Note: ClinVar primarily contains pathogenic repeat expansions (>20 GGGGCC repeats) and rare variants of uncertain significance. Most variants lack sufficient clinical evidence for definitive pathogenic classification.*

### AI-Based Predictions

#### AlphaMissense Missense Pathogenicity

**Summary**
| Metric | Count |
|--------|-------|
| Total predictions | 600+ |
| Likely Pathogenic | 100+ |
| Ambiguous | 150+ |
| Likely Benign | 350+ |

**Top 30 Likely Pathogenic Missense Variants (AlphaMissense)**

| Genomic Position | Protein Change | AM Pathogenicity | Class |
|------------------|----------------|------------------|-------|
| 9:27548239:A:C | F481L | 0.890 | Likely Pathogenic |
| 9:27548255:T:A | D476V | 0.931 | Likely Pathogenic |
| 9:27548256:C:G | D476H | 0.942 | Likely Pathogenic |
| 9:27548306:A:G | L459P | 0.987 | Likely Pathogenic |
| 9:27548327:G:T | A452D | 0.996 | Likely Pathogenic |
| 9:27548328:C:G | A452P | 0.996 | Likely Pathogenic |
| 9:27548334:C:G | A450P | 0.979 | Likely Pathogenic |
| 9:27548336:A:T | M449K | 0.970 | Likely Pathogenic |
| 9:27548342:A:T | I447K | 0.990 | Likely Pathogenic |
| 9:27548348:A:G | L445P | 0.994 | Likely Pathogenic |
| 9:27548351:T:A | D444V | 0.980 | Likely Pathogenic |
| 9:27548352:C:G | D444H | 0.983 | Likely Pathogenic |
| 9:27548355:C:G | G443R | 0.987 | Likely Pathogenic |
| 9:27548366:A:G | L439S | 0.985 | Likely Pathogenic |
| 9:27548372:A:G | L437P | 0.995 | Likely Pathogenic |
| 9:27548255:T:C | D476G | 0.803 | Likely Pathogenic |
| 9:27548261:T:A | E474V | 0.850 | Likely Pathogenic |
| 9:27548324:T:A | E453V | 0.985 | Likely Pathogenic |
| 9:27548269:A:C | S471R | 0.946 | Likely Pathogenic |
| 9:27548296:A:T | F462L | 0.950 | Likely Pathogenic |
| 9:27548290:A:T | F464L | 0.942 | Likely Pathogenic |
| 9:27548249:A:G | L478P | 0.912 | Likely Pathogenic |
| 9:27548306:A:C | L459R | 0.962 | Likely Pathogenic |
| 9:27548309:C:T | G458D | 0.866 | Likely Pathogenic |
| 9:27548309:C:A | G458V | 0.918 | Likely Pathogenic |
| 9:27548304:G:C | H460D | 0.897 | Likely Pathogenic |
| 9:27548362:C:T | E453K | 0.964 | Likely Pathogenic |
| 9:27548315:T:A | K456I | 0.750 | Likely Pathogenic |
| 9:27548320:T:A | K454N | 0.905 | Likely Pathogenic |
| 9:27548297:A:G | F462S | 0.918 | Likely Pathogenic |

#### Splice Effect Predictions
| Dataset | Count |
|---------|-------|
| SpliceAI | 0 |

*No SpliceAI predictions available in biobtree for C9orf72.*
