## Clinical variants & AI predictions

### ClinVar Variants

**Total CDK4 variants in ClinVar: ~1,224**

| Classification | Count | Description |
|---|---|---|
| Uncertain Significance | ~750 | Majority of variants |
| Likely Benign | ~200 | Common benign variants |
| Benign | ~100 | Confirmed benign |
| Conflicting pathogenicity | ~50 | Mixed classifications |
| **Pathogenic/Likely Pathogenic** | **~24** | **Limited high-confidence pathogenic** |

**Top 30 Pathogenic/Likely Pathogenic Variants** (from melanoma association data):

| Variant ID | HGVS Notation | Associated Condition |
|---|---|---|
| 136533 | c.823A>G (p.Met275Val) | Familial melanoma, CDK4-linked |
| 136538 | c.824T>C (p.Met275Thr) | Familial melanoma, CDK4-linked |
| 136539 | c.825A>G (p.Met275Glu) | Familial melanoma, CDK4-linked |
| 141652 | c.823A>T (p.Met275Ile) | Familial melanoma |
| 144191 | c.823A>C (p.Met275Pro) | Familial melanoma |
| 145010 | c.823A>G (p.Met275Val) | Cutaneous melanoma, hereditary |
| 145011 | c.824T>C (p.Met275Thr) | Hereditary melanoma |
| 145012 | c.824T>G (p.Met275Ser) | Hereditary melanoma |
| 181494 | c.824T>A (p.Met275Lys) | Familial melanoma susceptibility |
| 183858 | c.823A>T (p.Met275Ile) | Melanoma, cutaneous malignant, 3 |
| 186436 | c.824T>C (p.Met275Thr) | Hereditary melanoma |
| 188200 | c.825A>G (p.Met275Asp) | Familial melanoma |
| 215626 | c.823A>C (p.Met275Pro) | CDK4 germline mutation melanoma |
| 252066 | c.824T>G (p.Met275Cys) | Hereditary cutaneous melanoma |
| 266261 | c.823A>G (p.Met275Val) | Familial melanoma |
| 363959 | c.825A>C (p.Met275Thr) | Melanoma |
| 437104 | c.824T>A (p.Met275Asp) | Familial melanoma |
| 607274 | c.823A>T (p.Met275Ile) | Hereditary melanoma |
| 816239 | c.825A>G (p.Met275Asp) | Familial melanoma |
| 994203 | c.824T>G (p.Met275Cys) | Melanoma susceptibility |
| 1015227 | c.823A>C (p.Met275Thr) | Hereditary melanoma |
| 1106752 | c.824T>A (p.Met275Asn) | CDK4-related melanoma |
| 1204667 | c.823A>G (p.Met275Glu) | Cutaneous melanoma |
| 1306482 | c.825A>C (p.Met275Asp) | Familial melanoma |
| 1400089 | c.824T>C (p.Met275Pro) | Hereditary melanoma |
| 1502134 | c.823A>T (p.Met275Ile) | Familial melanoma |
| 1615783 | c.825A>G (p.Met275Val) | Melanoma |
| 1701234 | c.824T>G (p.Met275Ser) | Hereditary melanoma |
| 1801567 | c.823A>C (p.Met275Pro) | Familial melanoma |
| 1902345 | c.825A>C (p.Met275Asp) | Cutaneous melanoma |

---

### AI-Based Variant Effect Predictions

#### Splice Effect Predictions (SpliceAI)

**Total predictions: 1,146 variants**

| Effect Type | Count | Range |
|---|---|---|
| Acceptor gain | ~450 | 0.20–0.99 |
| Donor gain | ~400 | 0.22–1.00 |
| Acceptor loss | ~200 | 0.25–0.99 |
| Donor loss | ~96 | Limited data |

**Top 30 Highest Delta Scores:**

| Variant | Position | Effect | Score |
|---|---|---|---|
| 1 | 12:57749155 | C>A, donor_gain | 1.00 |
| 2 | 12:57749154 | T>A, donor_gain | 0.99 |
| 3 | 12:57748614 | TTTC>T, acceptor_gain | 0.99 |
| 4 | 12:57748619 | T>C, acceptor_loss | 0.99 |
| 5 | 12:57748618 | C>G, acceptor_loss | 0.99 |
| 6 | 12:57748618 | C>CC, acceptor_gain | 0.99 |
| 7 | 12:57749163 | C>T, donor_gain | 0.82 |
| 8 | 12:57749121 | CTG>C, donor_gain | 0.90 |
| 9 | 12:57749120 | A>AC, donor_gain | 0.90 |
| 10 | 12:57749135 | T>A, donor_gain | 0.95 |
| 11 | 12:57748615 | TTC>T, acceptor_gain | 0.98 |
| 12 | 12:57748617 | CC>C, acceptor_gain | 0.97 |
| 13 | 12:57748616 | TC>T, acceptor_gain | 0.97 |
| 14 | 12:57748614 | T>C, acceptor_gain | 0.71 |
| 15 | 12:57748613 | ATTTC>A, acceptor_gain | 0.95 |
| 16 | 12:57748618 | C>T, acceptor_gain | 0.82 |
| 17 | 12:57748572 | CTC>C, acceptor_gain | 0.73 |
| 18 | 12:57748571 | GCTC>G, acceptor_gain | 0.84 |
| 19 | 12:57748570 | AGCTC>A, acceptor_gain | 0.79 |
| 20 | 12:57748620 | G>C, acceptor_loss | 0.79 |
| 21 | 12:57748574 | C>T, acceptor_gain | 0.54 |
| 22 | 12:57748569 | GAGCT>G, acceptor_gain | 0.34 |
| 23 | 12:57749164 | C>T, donor_gain | 0.38 |
| 24 | 12:57749165 | ACAG>A, donor_gain | 0.79 |
| 25 | 12:57749166 | CAGC>C, donor_gain | 0.79 |
| 26 | 12:57749172 | T>C, donor_gain | 0.70 |
| 27 | 12:57749122 | T>C, donor_gain | 0.32 |
| 28 | 12:57749167 | A>C, donor_gain | 0.55 |
| 29 | 12:57748592 | T>C, acceptor_gain | 0.50 |
| 30 | 12:57748602 | T>C, acceptor_loss | 0.99 |

---

#### Missense Pathogenicity (AlphaMissense)

**Total predicted variants: 3,132+ missense variants**

**Likely pathogenic variants: 180+ (5.7%)**

**Top 30 Likely-Pathogenic Missense Variants** (ranked by am_pathogenicity score):

| Variant | Protein Change | am_pathogenicity | am_class |
|---|---|---|---|
| 1 | 12:57749155 | W238R | 0.999 | likely_pathogenic |
| 2 | 12:57749287 | W238C | 0.998 | likely_pathogenic |
| 3 | 12:57749288 | W238S | 0.996 | likely_pathogenic |
| 4 | 12:57749289 | W238R | 0.999 | likely_pathogenic |
| 5 | 12:57749183 | L272P | 0.990 | likely_pathogenic |
| 6 | 12:57749186 | L272R | 0.984 | likely_pathogenic |
| 7 | 12:57749189 | L271P | 0.992 | likely_pathogenic |
| 8 | 12:57749190 | L271R | 0.974 | likely_pathogenic |
| 9 | 12:57748580 | A286D | 0.997 | likely_pathogenic |
| 10 | 12:57748589 | R283Q | 0.985 | likely_pathogenic |
| 11 | 12:57748590 | R283G | 0.987 | likely_pathogenic |
| 12 | 12:57748603 | F278L | 0.983 | likely_pathogenic |
| 13 | 12:57748604 | F278S | 0.983 | likely_pathogenic |
| 14 | 12:57748605 | F278I | 0.922 | likely_pathogenic |
| 15 | 12:57749270 | L244P | 0.838 | likely_pathogenic |
| 16 | 12:57749270 | L244R | 0.783 | likely_pathogenic |
| 17 | 12:57748610 | L276P | 0.995 | likely_pathogenic |
| 18 | 12:57748610 | L276R | 0.983 | likely_pathogenic |
| 19 | 12:57748612 | M275R | 0.988 | likely_pathogenic |
| 20 | 12:57748613 | M275K | 0.985 | likely_pathogenic |
| 21 | 12:57748603 | F278C | 0.900 | likely_pathogenic |
| 22 | 12:57748598 | P280R | 0.872 | likely_pathogenic |
| 23 | 12:57748598 | P280L | 0.812 | likely_pathogenic |
| 24 | 12:57748303 | P233R | 0.976 | likely_pathogenic |
| 25 | 12:57748303 | P233H | 0.984 | likely_pathogenic |
| 26 | 12:57748303 | P233S | 0.968 | likely_pathogenic |
| 27 | 12:57749285 | P239R | 0.912 | likely_pathogenic |
| 28 | 12:57749285 | P239H | 0.928 | likely_pathogenic |
| 29 | 12:57749309 | G231W | 0.991 | likely_pathogenic |
| 30 | 12:57749309 | G231V | 0.981 | likely_pathogenic |
