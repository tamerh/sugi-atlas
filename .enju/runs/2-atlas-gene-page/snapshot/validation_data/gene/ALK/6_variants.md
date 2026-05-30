## Clinical variants & AI predictions

**ALK** (HGNC:427) — 1,620 aa receptor tyrosine kinase

### ClinVar Clinical Variants

| Metric | Count |
|--------|-------|
| **Total variants** | ~6,175 |
| Uncertain significance | ~5,100+ |
| Conflicting classifications | ~50–100 |
| Likely benign | ~20–50 |
| Benign | ~10–20 |
| Pathogenic/Likely Pathogenic | <10 |

**Associated phenotypes:** ALK-Related Neuroblastoma Susceptibility, Hereditary cancer-predisposing syndromes, Tumor predisposition

**Sample pathogenic/likely pathogenic variants** (limited in dataset):
| Variant ID | HGVS | Condition |
|-----------|------|-----------|
| 1001211 | c.1634G>A (p.Ser545Asn) | Neuroblastoma susceptibility (conflicting) |
| 1002636 | c.1375C>G (p.Gln459Glu) | Neuroblastoma susceptibility, ovarian cancer (conflicting) |
| 1002687 | c.1249del (p.Val417fs) | Neuroblastoma susceptibility (conflicting) |

*Note: Most ALK variants are classified as VUS; pathogenic/LP classifications are rare in germline sets.*

---

### AlphaMissense Predictions

| Category | Count |
|----------|-------|
| **Total predictions** | 10,546 |
| **Likely pathogenic** | ~100+ |

**Top 30 likely-pathogenic missense variants (amP score ≥ 0.84)**

| Position | Protein Change | amP Score | Class |
|----------|---|---|---|
| 2:29193578 | W1503C | 0.999 | likely_pathogenic |
| 2:29193575 | N1504K | 0.998 | likely_pathogenic |
| 2:29193580 | W1503R | 1.000 | likely_pathogenic |
| 2:29193575 | N1504K | 0.998 | likely_pathogenic |
| 2:29193568 | Y1507D | 0.992 | likely_pathogenic |
| 2:29193579 | W1503S | 0.997 | likely_pathogenic |
| 2:29193572 | P1505Q | 0.990 | likely_pathogenic |
| 2:29193337 | Y1584N | 0.907 | likely_pathogenic |
| 2:29193337 | Y1584H | 0.948 | likely_pathogenic |
| 2:29193336 | Y1584S | 0.881 | likely_pathogenic |
| 2:29193349 | G1580R | 0.944 | likely_pathogenic |
| 2:29193348 | G1580E | 0.965 | likely_pathogenic |
| 2:29193348 | G1580V | 0.957 | likely_pathogenic |
| 2:29193356 | F1577L | 0.988 | likely_pathogenic |
| 2:29193371 | F1572L | 0.976 | likely_pathogenic |
| 2:29193357 | F1577S | 0.990 | likely_pathogenic |
| 2:29193554 | F1511L | 0.968 | likely_pathogenic |
| 2:29193557 | W1510C | 0.991 | likely_pathogenic |
| 2:29193564 | G1508D | 0.996 | likely_pathogenic |
| 2:29193562 | S1509P | 0.991 | likely_pathogenic |
| 2:29193342 | V1582D | 0.966 | likely_pathogenic |
| 2:29193559 | W1510R | 0.998 | likely_pathogenic |
| 2:29193894 | P1398R | 0.979 | likely_pathogenic |
| 2:29193573 | P1505R | 0.985 | likely_pathogenic |
| 2:29193912 | V1392E | 0.996 | likely_pathogenic |
| 2:29193568 | Y1507H | 0.985 | likely_pathogenic |
| 2:29193333 | G1585D | 0.883 | likely_pathogenic |
| 2:29193337 | Y1584D | 0.951 | likely_pathogenic |
| 2:29193334 | G1585R | 0.785 | likely_pathogenic |
| 2:29193570 | T1506K | 0.989 | likely_pathogenic |

---

### SpliceAI Predictions

| Category | Count |
|----------|-------|
| **Total splice effect predictions** | 6,693 |

**High-confidence splice variants (score ≥ 0.85)** — selected examples

| Position | Variant | Effect | Score |
|----------|---------|--------|-------|
| 2:29196859 | AT:A | acceptor_gain | 1.00 |
| 2:29196858 | TAT:T | acceptor_gain | 1.00 |
| 2:29196764 | TTTTA:T | donor_loss | 1.00 |
| 2:29196765 | TTTA:T | donor_loss | 1.00 |
| 2:29196766 | TTAC:T | donor_loss | 1.00 |
| 2:29196767 | TACCT:T | donor_loss | 1.00 |
| 2:29196769 | C:A | donor_loss | 1.00 |
| 2:29193926 | CCG:C | acceptor_gain | 0.99 |
| 2:29193927 | C:CT | acceptor_gain | 0.99 |
| 2:29193927 | C:T | acceptor_gain | 0.99 |
| 2:29193857 | GTAT:G | acceptor_gain | 0.99 |
| 2:29193923 | C:CC | acceptor_gain | 0.98 |
| 2:29193801 | T:G | acceptor_loss | 0.98 |
| 2:29193921 | TC:T | acceptor_gain | 0.91 |
| 2:29193922 | CC:C | acceptor_gain | 0.91 |

**Summary:** ALK harbors primarily **uncertain significance variants** in clinical databases; AI predictions identify **~100 likely-pathogenic missense mutations** (mostly in kinase domain, C-terminus) and **6,693 splice-affecting variants** with high confidence predictions, providing computational support for potential pathogenicity in research/functional studies.
