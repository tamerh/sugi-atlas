## Clinical variants & AI predictions

### ClinVar Summary

**Total variants:** ~53 variants annotated for BAX gene

**Classification breakdown:**
| Classification | Count | Representative |
|---|---|---|
| Pathogenic | 3 | c.199G>A (p.Gly67Arg), c.115_121del (p.Gly39fs), c.121del (p.Glu41fs) |
| Likely Pathogenic | 0 | — |
| Uncertain significance | ~40 | Most recent submissions (2021-2024) by Ambry Genetics |
| Likely Benign | ~5 | Synonymous/intronic variants |
| Benign | ~5 | Intronic variants |

**Top pathogenic/likely pathogenic variants:**
| Variant ID | HGVS notation | Associated condition |
|---|---|---|
| 9513 | c.199G>A (p.Gly67Arg) | Acute T-cell leukemia, T-cell acute lymphoblastic leukemia (somatic) |
| 9514 | c.115_121del (p.Gly39fs) | T-cell acute lymphoblastic leukemia (somatic) |
| 9512 | c.121del (p.Glu41fs) | Colon carcinoma (somatic) |

*Note: Most remaining 50 variants are classified as Uncertain Significance (UVS) from clinical testing providers, with limited functional evidence.*

---

### AlphaMissense predictions

**Total variants:** 1,250 missense predictions for BAX protein

**Likely pathogenic (am_class=="likely_pathogenic"):** 100+ predictions with high confidence

**Top 30 likely-pathogenic variants:**
| Position | Variant | am_pathogenicity | Position | Variant | am_pathogenicity |
|---|---|---|---|---|---|
| S15R | A>C | 0.758 | G23R | G>C | 0.998 |
| S15R | C>A | 0.758 | G23W | G>T | 0.999 |
| S16F | C>T | 0.672 | G23E | G>A | 0.999 |
| E17K | G>A | 0.652 | G23A | G>C | 0.735 |
| E17V | A>T | 0.700 | G23V | G>T | 0.991 |
| I19L | A>C | 0.738 | A24P | G>C | 0.980 |
| I19V | A>G | 0.584 | A24D | C>A | 0.685 |
| I19F | A>T | 0.963 | L25P | T>C | 0.923 |
| I19N | T>A | 0.990 | L26M | T>A | 0.887 |
| I19T | T>C | 0.992 | L26V | T>G | 0.809 |
| I19S | T>G | 0.994 | L26S | T>C | 0.999 |
| I19M | C>G | 0.789 | L26W | T>G | 0.996 |
| M20K | T>A | 0.911 | L26F | G>C | 0.979 |
| M20T | T>C | 0.915 | L27I | C>A | 0.783 |
| M20R | T>G | 0.884 | L27V | C>G | 0.849 |

*Highest confidence: I19T (0.963), I19N (0.990), I19T (0.992), I19S (0.994), G23W (0.999), G23E (0.999), L26S (0.999), L26W (0.996), L27P (0.999)*

---

### SpliceAI predictions

**Total splice variants:** ~837 predictions (donor/acceptor gain/loss effects)

**Top 30 high-impact splice predictions:**
| Position | Variant | Effect | Score |
|---|---|---|---|
| 19:48955118 | G>C | donor_gain | 0.50 |
| 19:48955118 | G>GT | donor_gain | 1.00 |
| 19:48955118 | G>T | donor_gain | 0.92 |
| 19:48954932 | G>GT | donor_gain | 0.99 |
| 19:48955199 | G>GT | donor_gain | 0.99 |
| 19:48955199 | G>T | donor_gain | 0.93 |
| 19:48954972 | G>GT | donor_gain | 0.95 |
| 19:48954972 | G>T | donor_gain | 0.91 |
| 19:48955191 | G>T | donor_gain | 0.91 |
| 19:48954917 | G>GT | donor_gain | 0.95 |
| 19:48954960 | G>GT | donor_gain | 0.98 |
| 19:48954960 | G>GG | donor_gain | 0.98 |
| 19:48954961 | G>GG | donor_gain | 0.95 |
| 19:48954961 | G>GGG | donor_gain | 0.97 |
| 19:48954962 | G>GG | donor_gain | 0.95 |
| 19:48954959 | G>GGGG | donor_gain | 0.98 |
| 19:48955216 | C>T | donor_gain | 0.95 |
| 19:48955220 | C>T | donor_gain | 0.81 |
| 19:48955226 | C>T | donor_gain | 0.71 |
| 19:48955256 | C>T | donor_gain | 0.90 |

*Highest-risk variants concentrate in intronic regions near exon boundaries; many donor_gain predictions (>0.90) indicate potential cryptic splice site creation.*
