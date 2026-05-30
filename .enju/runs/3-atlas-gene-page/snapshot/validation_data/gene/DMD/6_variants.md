## Clinical variants & AI predictions

### Clinical Variants (ClinVar)

**Summary**
- **Total variants**: ~11,261
- **Associated conditions**: Duchenne muscular dystrophy (DMD), Becker muscular dystrophy (BMD), dilated cardiomyopathy

**Classification Breakdown** (sample-based estimate from reviewed entries)
| Classification | Status |
|---|---|
| Uncertain significance | Predominant (~50-60%) |
| Conflicting classifications | ~10-15% |
| Likely benign | ~10-15% |
| Pathogenic / Likely pathogenic | ~5-10% |
| Benign | <5% |

**Top Pathogenic/Likely Pathogenic Variants** (representative samples)

| Variant ID | HGVS Notation | Condition |
|---|---|---|
| 1001786 | NM_004006.3(DMD):c.2400C>A (p.Ser800Arg) | Duchenne/Becker muscular dystrophy |
| 1003542 | NM_004006.3(DMD):c.5684A>T (p.Asp1895Val) | Duchenne/Becker muscular dystrophy, dilated cardiomyopathy |
| 1003878 | NM_004006.3(DMD):c.961-4del | Duchenne/Becker muscular dystrophy |
| 1003362 | NM_004006.3(DMD):c.10951G>A (p.Asp3651Asn) | Duchenne/Becker muscular dystrophy |
| 1004650 | NM_004006.3(DMD):c.4705T>C (p.Cys1569Arg) | Duchenne/Becker muscular dystrophy |

---

### AI-Based Variant Effect Predictions

**Splice Effect Predictions (SpliceAI)**
- **Total predictions**: 7,041
- **Effect types**: donor_gain, donor_loss, acceptor_gain, acceptor_loss
- **Score range**: 0.20–0.95

**Top 30 SpliceAI Predictions (by score)**

| Position | Variant | Gene | Effect | Score |
|---|---|---|---|---|
| 1 | X:31121931:C:CC | DMD | acceptor_gain | 0.95 |
| 2 | X:31121927:TGTC:T | DMD | acceptor_gain | 0.90 |
| 3 | X:31121928:GTC:G | DMD | acceptor_gain | 0.86 |
| 4 | X:31121766:A:AC | DMD | donor_gain | 0.86 |
| 5 | X:31121767:C:CC | DMD | donor_gain | 0.86 |
| 6 | X:31121758:AAACT:A | DMD | donor_gain | 0.72 |
| 7 | X:31126709:A:T | DMD | acceptor_gain | 0.95 |
| 8 | X:31121929:TC:T | DMD | acceptor_gain | 0.81 |
| 9 | X:31121930:CC:C | DMD | acceptor_gain | 0.81 |
| 10 | X:31126708:C:CT | DMD | acceptor_gain | 0.82 |
| 11 | X:31121276:G:C | DMD | donor_gain | 0.54 |
| 12 | X:31121762:T:TA | DMD | donor_gain | 0.52 |
| 13 | X:31121768:T:C | DMD | donor_gain | 0.60 |
| 14 | X:31126639:A:T | DMD | acceptor_gain | 0.76 |
| 15 | X:31126676:T:TC | DMD | acceptor_gain | 0.22 |
| 16–30 | (additional variants with scores 0.20–0.75) | DMD | mixed | 0.20–0.75 |

---

**AlphaMissense Pathogenicity Predictions**
- **Total likely_pathogenic variants**: 100+ (pagination available)
- **Pathogenicity score range**: 0.565–1.000
- **Classification**: am_class="likely_pathogenic"

**Top 30 Likely Pathogenic Missense Variants** (by am_pathogenicity score)

| Position | Genomic Position | Protein Variant | am_pathogenicity | Interpretation |
|---|---|---|---|---|
| 1 | X:31147294:A:G | L3593P | 1.000 | Highest confidence |
| 2 | X:31147290:C:A | R3594S | 0.997 | Very high confidence |
| 3 | X:31147290:C:G | R3594S | 0.997 | Very high confidence |
| 4 | X:31147296:C:A | R3592S | 0.997 | Very high confidence |
| 5 | X:31147296:C:G | R3592S | 0.997 | Very high confidence |
| 6 | X:31147294:A:C | L3593R | 0.997 | Very high confidence |
| 7 | X:31147291:C:A | R3594M | 0.995 | Very high confidence |
| 8 | X:31147291:C:G | R3594T | 0.995 | Very high confidence |
| 9 | X:31147297:C:A | R3592M | 0.995 | Very high confidence |
| 10 | X:31147297:C:G | R3592T | 0.995 | Very high confidence |
| 11 | X:31147282:A:G | L3597P | 0.998 | Very high confidence |
| 12 | X:31147285:A:G | L3596P | 0.999 | Very high confidence |
| 13 | X:31147294:A:T | L3593Q | 0.999 | Very high confidence |
| 14 | X:31147282:A:T | L3597Q | 0.991 | Very high confidence |
| 15 | X:31147298:T:C | R3592G | 0.992 | Very high confidence |
| 16 | X:31147288:T:G | Q3595P | 0.982 | Very high confidence |
| 17 | X:31147282:A:C | L3597R | 0.980 | Very high confidence |
| 18 | X:31134149:A:G | L3656S | 0.977 | High confidence |
| 19 | X:31147295:G:C | L3593V | 0.976 | High confidence |
| 20 | X:31147285:A:C | L3596R | 0.973 | High confidence |
| 21 | X:31147292:T:A | R3594W | 0.965 | High confidence |
| 22 | X:31147300:T:G | H3591P | 0.955 | High confidence |
| 23 | X:31134116:A:G | F3667S | 0.952 | High confidence |
| 24 | X:31134116:A:C | F3667C | 0.941 | High confidence |
| 25 | X:31134128:A:G | L3663P | 0.978 | High confidence |
| 26 | X:31134128:A:T | L3663H | 0.970 | High confidence |
| 27 | X:31146310:A:C | S3634R | 0.919 | High confidence |
| 28 | X:31146310:A:T | S3634R | 0.919 | High confidence |
| 29 | X:31134128:A:C | L3663R | 0.926 | High confidence |
| 30 | X:31134115:G:C | F3667L | 0.984 | Very high confidence |
