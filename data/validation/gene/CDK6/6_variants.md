## Clinical variants & AI predictions

### ClinVar Classification Summary

| Classification | Count |
|---|---|
| Pathogenic | 1 |
| Likely Pathogenic | 0 |
| Uncertain Significance | 12 |
| Likely Benign | 11 |
| Benign | 15 |
| Other (conflicting/unclassified) | 14 |
| **Total** | **53** |

### Top Pathogenic/Likely Pathogenic ClinVar Variants

| ClinVar ID | HGVS | Condition/Classification |
|---|---|---|
| 157508 | NM_001145306.2(CDK6):c.589G>A (p.Ala197Thr) | Pathogenic |
| 434660 | NM_001145306.2(CDK6):c.328G>A (p.Asp110Asn) | Conflicting classifications of pathogenicity |

*Only 2 variants with pathogenic/conflicting pathogenic classifications in ClinVar. Most disease-associated signals appear in AI predictions rather than clinical annotations.*

### AlphaMissense Missense Pathogenicity Predictions

| Metric | Count |
|---|---|
| Total Likely Pathogenic | 177 |

#### Top 30 Likely Pathogenic Variants (by am_pathogenicity score)

| Protein Variant | am_pathogenicity | Location |
|---|---|---|
| R288S | 1.000 | 7:92615257 |
| W243R | 1.000 | 7:92618179 |
| W243C | 0.999 | 7:92618177 |
| R288I | 0.998 | 7:92615258 |
| W243S | 0.998 | 7:92618178 |
| L277P | 0.998 | 7:92618076 |
| R288T | 0.999 | 7:92615258 |
| L276P | 0.999 | 7:92618079 |
| F254L/L | 0.990 | 7:92618146, 7:92618144 |
| G273D | 0.996 | 7:92618088 |
| L295P | 0.990 | 7:92615237 |
| A294D | 0.997 | 7:92615240 |
| L228P | 0.998 | 7:92623051 |
| Y299D | 0.996 | 7:92615226 |
| L276Q | 0.990 | 7:92618079 |
| I231S/T/N | 0.990-0.991 | 7:92623042 |
| P244T | 0.950 | 7:92618176 |
| L249P | 0.957 | 7:92618160 |
| F254C | 0.943 | 7:92618145 |
| F300S | 0.995 | 7:92615222 |
| P238Q | 0.989 | 7:92618193 |
| L277H | 0.996 | 7:92618076 |
| F283S | 0.993 | 7:92615273 |
| Q227P | 0.999 | 7:92623054 |
| V234E | 0.870 | 7:92618205 |
| H297R | 0.984 | 7:92615231 |
| H297Q | 0.973 | 7:92615230 |
| R288G | 0.995 | 7:92615259 |
| G236E | 0.992 | 7:92618199 |
| L277F | 0.782 | 7:92618077 |

### SpliceAI Splice Effect Predictions

| Metric | Count |
|---|---|
| Total SpliceAI Variants | 2,308 |

#### Top 30 Splice Effect Variants (by delta score)

| Location | Effect | Score |
|---|---|---|
| 7:92615282:CACTT:C | acceptor_gain | 1.0000 |
| 7:92615284:CTT:C | acceptor_gain | 1.0000 |
| 7:92615287:C:CC | acceptor_gain | 1.0000 |
| 7:92614796:G:C | donor_gain | 0.9900 |
| 7:92614991:T:A | donor_gain | 0.9900 |
| 7:92615285:TC:T | acceptor_loss | 0.9900 |
| 7:92615285:TT:T | acceptor_gain | 0.9900 |
| 7:92615288:T:A | acceptor_loss | 0.9900 |
| 7:92615287:C:CA | acceptor_loss | 0.9900 |
| 7:92614791:A:C | donor_gain | 0.9400 |
| 7:92614814:C:CC | donor_gain | 0.8300 |
| 7:92615292:A:C | acceptor_gain | 0.9400 |
| 7:92615292:A:AC | acceptor_gain | 0.9600 |
| 7:92615289:G:C | acceptor_gain | 0.8900 |
| 7:92615289:G:GC | acceptor_gain | 0.9300 |
| 7:92614807:A:AC | donor_gain | 0.9100 |
| 7:92614748:A:AC | donor_gain | 0.9300 |
| 7:92614785:A:C | donor_gain | 0.9300 |
| 7:92614780:A:AC | donor_gain | 0.9500 |
| 7:92615283:ACTTC:A | acceptor_gain | 0.6400 |
| 7:92614753:CAAA:C | donor_gain | 0.6200 |
| 7:92614904:T:TC | acceptor_gain | 0.4900 |
| 7:92614859:T:TA | donor_gain | 0.7500 |
| 7:92615091:C:CT | donor_gain | 0.7600 |
| 7:92615037:T:TA | donor_gain | 0.7500 |
| 7:92610181:T:A | donor_gain | 0.7300 |
| 7:92610174:G:C | donor_gain | 0.7200 |
| 7:92614904:T:TC | acceptor_gain | 0.4900 |
| 7:92615201:TTGC:T | donor_gain | 0.8000 |
| 7:92615245:G:GT | acceptor_gain | 0.6600 |

**Summary**: CDK6 shows limited pathogenic variants in clinical databases (1 confirmed pathogenic, 1 conflicting), but extensive AI predictions reveal 177 likely-pathogenic missense variants and 2,308 potential splice-affecting variants, indicating substantial predicted functional constraint.
