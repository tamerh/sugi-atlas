## Clinical variants & AI predictions

### Clinical variant annotations (ClinVar)

**Variant count by classification (29 total)**

| Classification | Count |
|---|---|
| Benign | 5 |
| Likely benign | 2 |
| VUS (Uncertain significance) | 13 |
| Pathogenic | 5 |
| Unclassified | 4 |

**Top pathogenic/likely pathogenic variants**

| Variant ID | HGVS | Classification | Phenotype |
|---|---|---|---|
| 3906943 | CD274 IVS4DS, G-A, +1 (splice) | Pathogenic | Autoimmune disease, multisystem, infantile-onset, 5 (MONDO:0979235) |
| 1703576 | GRCh37 9p24.3-q13 (CNV gain ~68 Mb) | Pathogenic | Bradycardia (HP:0001662) |
| 1706515 | GRCh37 9p24.3-23 (CNV loss ~10 Mb) | Pathogenic | — |
| 442671 | GRCh37 9p24.3-23 (CNV loss ~13 Mb) | Pathogenic | — |
| 563686 | GRCh37 9p24.3-q21.11 (CNV gain ~70 Mb) | Pathogenic | — |

---

### AI-based variant effect predictions

**AlphaMissense missense pathogenicity**

- **Total missense variants**: 1923
- **Likely pathogenic predictions**: ~155+ (filtered subset shown below)

**Top 30 likely-pathogenic missense predictions (am_pathogenicity score)**

| Variant | Position | Protein change | am_pathogenicity | Position | Protein change | am_pathogenicity |
|---|---|---|---|---|---|---|
| 9:5457144:T:A | | C40S | 0.991 | 9:5457328:T:A | | I101N | 0.917 |
| 9:5457144:T:C | | C40R | 0.986 | 9:5457322:T:C | | L99P | 0.982 |
| 9:5457150:T:C | | F42L | 0.984 | 9:5457322:T:A | | L99H | 0.978 |
| 9:5457151:T:G | | F42C | 0.985 | 9:5457190:T:A | | V55D | 0.973 |
| 9:5457152:C:A | | F42L | 0.984 | 9:5457278:G:C | | R84S | 0.974 |
| 9:5457310:G:T | | G95V | 0.973 | 9:5457278:G:T | | R84S | 0.974 |
| 9:5457195:T:A | | W57R | 0.970 | 9:5457328:T:C | | I101T | 0.964 |
| 9:5457195:T:C | | W57R | 0.970 | 9:5457328:T:G | | I101S | 0.966 |
| 9:5457197:G:C | | W57C | 0.991 | 9:5457322:T:G | | L99R | 0.948 |
| 9:5457197:G:T | | W57C | 0.991 | 9:5457286:T:G | | L87R | 0.904 |
| 9:5457145:G:C | | C40S | 0.991 | 9:5457337:T:A | | V104E | 0.902 |
| 9:5457081:T:C | | F19L | 0.959 | 9:5457309:G:A | | G95R | 0.903 |
| 9:5457316:C:A | | A97D | 0.967 | 9:5457309:G:C | | G95R | 0.903 |
| 9:5457315:G:C | | A97P | 0.956 | 9:5457310:G:A | | G95E | 0.929 |
| 9:5457277:G:C | | R84T | 0.940 | 9:5457184:T:C | | L53P | 0.921 |

**Splice site variants**: Limited direct predictions available; note IVS4DS:G-A (+1) classified as Pathogenic in ClinVar with associated autoimmune phenotype.
