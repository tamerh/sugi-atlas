## Clinical variants & AI predictions

**Clinical Variants (ClinVar)**

| Classification | Count |
|---|---|
| Pathogenic | ~1,800+ |
| Likely Pathogenic | ~200+ |
| Uncertain Significance (VUS) | ~12,000+ |
| Conflicting Classifications | ~6,000+ |
| Likely Benign | ~150+ |
| Benign | ~30+ |
| **Total** | **~21,181** |

**Top 30 Pathogenic/Likely Pathogenic Variants (BRCA2)**

| ClinVar ID | HGVS Notation | Variant Type | Classification |
|---|---|---|---|
| 1012157 | NM_000059.4:c.9163del (p.Leu3055fs) | Deletion | Pathogenic |
| 1012158 | NM_000059.4:c.5934del (p.Phe1978fs) | Deletion | Pathogenic |
| 1012159 | NM_000059.4:c.5566_5567inv (p.His1856Cys) | Inversion | Pathogenic |
| 1012160 | NM_000059.4:c.5362del (p.Ser1788fs) | Deletion | Pathogenic |
| 1012161 | NM_000059.4:c.5297del (p.Asn1766fs) | Deletion | Pathogenic |
| 1012162 | NM_000059.4:c.1561del (p.Ser521fs) | Deletion | Pathogenic |
| 1012163 | NM_000059.4:c.1053del (p.Lys351fs) | Deletion | Pathogenic/Likely Pathogenic |
| 1012164 | NM_000059.4:c.728del (p.Asn243fs) | Deletion | Pathogenic |
| 1012165 | NM_000059.4:c.691_692delinsGA (p.Ser231Asp) | Indel | Pathogenic |
| 1012166 | NM_000059.4:c.2588del (p.Asn863fs) | Deletion | Pathogenic |
| 1012167 | NM_000059.4:c.7177del (p.Lys2392_Met2393insTer) | Deletion | Pathogenic |
| 1012168 | NM_000059.4:c.10248del (p.Lys3416fs) | Deletion | Pathogenic |
| 1012202 | NM_000059.4:c.8423_8427delinsA (p.Leu2808fs) | Indel | Pathogenic |
| 1012203 | NM_000059.4:c.8487+2T>G | Splice Site | Pathogenic |
| 1012631 | NM_000059.4:c.1490_1493del (p.Ser497fs) | Deletion | Pathogenic |
| 1027606 | NM_000059.4:c.4057del (p.Glu1353fs) | Deletion | Likely Pathogenic |

**AI Predictions: AlphaMissense (Missense Pathogenicity)**

**Total Predictions:** 22,763 variants  
**Likely Pathogenic Predictions:** ~1,200+ variants

**Top 30 Likely-Pathogenic Missense Variants (by am_pathogenicity score)**

| Genomic Position | Protein Change | am_pathogenicity | Effect |
|---|---|---|---|
| 13:32319100 | W31R | 0.992 | Critical loss-of-function |
| 13:32319100 | W31R | 0.992 | Critical loss-of-function |
| 13:32326562 | W194R | 0.997 | Critical loss-of-function |
| 13:32326562 | W194R | 0.997 | Critical loss-of-function |
| 13:32326564 | W194C | 0.991 | Critical loss-of-function |
| 13:32319102 | W31C | 0.960 | Critical loss-of-function |
| 13:32326577 | A199P | 0.950 | Protein structural disruption |
| 13:32319105 | F32L | 0.955 | Hydrophobic core mutation |
| 13:32319104 | F32S | 0.938 | Loss of function |
| 13:32326568 | S196R | 0.990 | Charge reversal |
| 13:32319082 | G25R | 0.917 | Bulky residue insertion |
| 13:32319304 | F15C | 0.790 | Aromatic loss |
| 13:32319100 | W31G | 0.862 | Large aromatic loss |
| 13:32319113 | L35P | 0.981 | Rigidity introduction |
| 13:32326253 | S163R | 0.750 | Phosphorylation site disruption |
| 13:32326259 | F165L | 0.896 | Hydrophobic substitution |
| 13:32326535 | G185R | 0.935 | Charge insertion |
| 13:32326521 | I180S | 0.897 | Hydrophobic loss |
| 13:32319089 | I27K | 0.920 | Charge reversal |
| 13:32319091 | S28R | 0.969 | Phosphorylation disruption |
| 13:32326583 | P201Q | 0.723 | Proline flexibility loss |
| 13:32319198 | K63N | 0.904 | Charge loss |
| 13:32319113 | L35H | 0.922 | Hydrophobic loss |
| 13:32326572 | S197Y | 0.813 | Aromatic substitution |
| 13:32319082 | G25A | 0.569 | Backbone disruption |
| 13:32326549 | D189E | 0.826 | Conservative change |
| 13:32319493 | F15S | 0.918 | Aromatic loss |
| 13:32326529 | S183R | 0.977 | Charge reversal |
| 13:32319093 | S28R | 0.969 | Phosphorylation disruption |
| 13:32326578 | A199D | 0.953 | Charge introduction |

**AI Predictions: SpliceAI (Splice Effect Predictions)**

**Total Predictions:** 3,855 variants  
**Distribution:** Primarily donor gain and donor loss events

**Top 30 SpliceAI Predictions (by impact score)**

| Position | Variant | Effect | Score |
|---|---|---|---|
| 13:32315088 | G:C | donor_gain | 0.99 |
| 13:32315136 | C:A | donor_gain | 0.99 |
| 13:32315313 | A:AC | donor_gain | 0.93 |
| 13:32315314 | C:CC | donor_gain | 0.93 |
| 13:32315262 | T:TA | donor_gain | 0.90 |
| 13:32315167 | AGGT:A | donor_gain | 0.90 |
| 13:32315312 | CA:C | donor_gain | 0.90 |
| 13:32315316 | C:CA | donor_gain | 0.86 |
| 13:32315322 | C:CA | donor_gain | 0.86 |
| 13:32315288 | A:AC | donor_gain | 0.87 |
| 13:32315289 | A:C | donor_gain | 0.89 |
| 13:32315237 | T:C | donor_gain | 0.71 |
| 13:32315194 | C:T | donor_gain | 0.73 |
| 13:32315314 | CT:C | donor_gain | 0.73 |
| 13:32315567 | C:CA | donor_gain | 0.77 |
| 13:32315309 | A:T | donor_gain | 0.76 |
| 13:32315571 | G:A | donor_gain | 0.69 |
| 13:32315572 | G:A | donor_gain | 0.75 |
| 13:32315198 | C:CT | donor_gain | 0.69 |
| 13:32315167 | AGGTC:A | donor_gain | 0.86 |
| 13:32315312 | CACT:C | donor_gain | 0.61 |
| 13:32315312 | C:T | donor_gain | 0.55 |
| 13:32315310 | C:T | donor_gain | 0.55 |
| 13:32315313 | ACT:A | donor_gain | 0.69 |
| 13:32315314 | CTC:C | donor_gain | 0.69 |
| 13:32315315 | T:A | donor_gain | 0.54 |
| 13:32315315 | T:C | donor_gain | 0.54 |
| 13:32315312 | CACTC:C | donor_gain | 0.59 |
| 13:32315342 | A:AC | donor_gain | 0.65 |
| 13:32315343 | C:CC | donor_gain | 0.65 |
