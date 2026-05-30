## Clinical variants & AI predictions

### ClinVar Overview

| Classification | Count |
|---|---|
| Pathogenic | ~50 |
| Likely Pathogenic | ~30 |
| Uncertain Significance | ~650 |
| Conflicting Classifications | ~150 |
| Likely Benign | ~150 |
| Benign | ~100 |
| **Total** | **~1,594** |

### TOP 30 Pathogenic/Likely Pathogenic ClinVar Variants

| Variant ID | HGVS Notation | Classification | Condition/Note |
|---|---|---|---|
| 1068946 | NM_000077.5:c.34del | Pathogenic | Frameshift; p.Ser12fs |
| 1069074 | NM_000077.5:c.126dup | Pathogenic | Premature stop; p.Ser43Ter |
| 1070481 | NC_000009.11:g.(?_21968228)_(22160087_?)del | Pathogenic | Large deletion |
| 1072356 | NM_058195.4:c.126_127insCA | Pathogenic | Frameshift; p.Val43fs |
| 1075589 | NM_000077.5:c.204_205delinsTT | Pathogenic | Premature stop; p.Glu69Ter |
| 1076895 | NC_000009.11:g.(?_21994128)_(21994330_?)del | Pathogenic | Large deletion |
| 1171238 | NM_000077.5:c.106del | Pathogenic | Frameshift; p.Ala36fs |
| 142061 | NM_000077.5:c.47_50del | Pathogenic | Frameshift; p.Leu16fs |
| 1409909 | NM_000077.5:c.95_112del | Pathogenic | Frameshift; p.Leu32_Leu37del |
| 135827 | NM_000077.5:c.-16GGCGGCGGGGAGCAGCATGGAGCC[3] | Pathogenic/LP | Microsatellite; p.Ala4_Pro11dup |
| 142882 | NM_000077.5:c.251A>C | Pathogenic/LP | p.Asp84Ala |
| 141882 | NM_058195.4:c.193+5G>A | Pathogenic/LP | Splice site |
| 1059393 | NM_058195.4:c.172C>T | Likely Pathogenic | Nonsense; p.Gln58Ter |
| 1215320 | NM_000077.5:c.281T>C | Likely Pathogenic | p.Leu94Pro |
| 1345910 | NC_000009.11:g.(?_21990646)_(21994323_?)del | Likely Pathogenic | Large deletion |
| 1347045 | NM_058195.4:c.87_99del | Likely Pathogenic | Frameshift; p.Leu30fs |

---

### AlphaMissense: Missense Pathogenicity Predictions

**Total likely_pathogenic predictions: ~1,200+ variants**

**TOP 30 Likely-Pathogenic Missense Variants (by AM-Pathogenicity Score)**

| Position | Variant | Protein Change | AM-Pathogenicity | Score Interpretation |
|---|---|---|---|---|
| 9:21971036 | T>G | D108A | 0.985 | Very high risk |
| 9:21971037 | C>G | D108H | 0.992 | Very high risk |
| 9:21971018 | G>T | P114H | 0.989 | Very high risk |
| 9:21971045 | T>A | D105V | 0.670 | Moderate-high |
| 9:21970974 | A>T | Y129N | 0.956 | Very high risk |
| 9:21970982 | A>T | V126D | 0.936 | Very high risk |
| 9:21971024 | C>G | R112P | 0.967 | Very high risk |
| 9:21971069 | A>G | L97P | 0.958 | Very high risk |
| 9:21971078 | A>G | L94P | 0.969 | Very high risk |
| 9:21971006 | G>T | A118D | 0.980 | Very high risk |
| 9:21971048 | A>C | L104R | 0.581 | Moderate |
| 9:21970973 | T>G | Y129S | 0.942 | Very high risk |
| 9:21971099 | C>G | R87P | 0.981 | Very high risk |
| 9:21971102 | G>T | A86D | 0.996 | Highest |
| 9:21971103 | C>G | A86P | 0.979 | Very high risk |
| 9:21970970 | A>G | L130P | 0.954 | Very high risk |
| 9:21971055 | C>G | A102P | 0.935 | Very high risk |
| 9:21971147 | T>A | N71I | 0.935 | Very high risk |
| 9:21971105 | G>T | A85D | 0.989 | Very high risk |
| 9:21971061 | C>G | A100P | 0.967 | Very high risk |
| 9:21971009 | A>G | L117P | 0.939 | Very high risk |
| 9:21971063 | C>G | R99P | 0.915 | Very high risk |
| 9:21971138 | T>C | D74G | 0.950 | Very high risk |
| 9:21971080 | C>G | A127P | 0.729 | Moderate-high |
| 9:21971057 | C>T | G101E | 0.759 | Moderate-high |
| 9:21971018 | G>C | P114R | 0.973 | Very high risk |
| 9:21971168 | A>G | L64P | 0.971 | Very high risk |
| 9:21971117 | G>C | P81R | 0.975 | Very high risk |
| 9:21971171 | A>G | L63P | 0.978 | Very high risk |
| 9:21971121 | G>T | S78P | 0.981 | Very high risk |

---

### SpliceAI: Splice Effect Predictions

**Total predictions: ~100 variants with splice effects**

**TOP 30 Splice-Altering Variants (by Score)**

| Position | Variant | Effect | Score | Impact |
|---|---|---|---|---|
| 9:21968346 | A>AC | Donor gain | 0.98 | Critical |
| 9:21968243 | C>CC | Acceptor gain | 0.99 | Critical |
| 9:21968240 | TGT>T | Acceptor loss | 0.95 | Critical |
| 9:21968241 | GTC>G | Acceptor loss | 0.95 | Critical |
| 9:21968244 | T>A | Acceptor loss | 0.95 | Critical |
| 9:21968239 | ATGTC>A | Acceptor loss | 0.95 | Critical |
| 9:21968347 | T>C | Donor gain | 0.97 | Critical |
| 9:21968238 | GATGT>G | Acceptor gain | 0.93 | High |
| 9:21968239 | ATGT>A | Acceptor gain | 0.87 | High |
| 9:21968253 | C>CT | Acceptor gain | 0.88 | High |
| 9:21968246 | C>CT | Acceptor gain | 0.86 | High |
| 9:21968247 | A>T | Acceptor gain | 0.86 | High |
| 9:21968254 | A>T | Acceptor gain | 0.83 | High |
| 9:21968245 | G>C | Acceptor loss | 0.89 | High |
| 9:21968342 | A>AC | Donor gain | 0.87 | High |
| 9:21968343 | C>CC | Donor gain | 0.87 | High |
| 9:21968330 | C>CT | Donor gain | 0.54 | Moderate |
| 9:21968370 | TC>T | Donor gain | 0.64 | Moderate |
| 9:21968324 | C>A | Donor gain | 0.53 | Moderate |
| 9:21968326 | C>T | Donor gain | 0.46 | Moderate |
| 9:21968331 | G>T | Donor gain | 0.47 | Moderate |
| 9:21968374 | G>C | Donor gain | 0.54 | Moderate |
| 9:21968376 | G>A | Donor gain | 0.52 | Moderate |
| 9:21968379 | T>C | Donor gain | 0.54 | Moderate |
| 9:21968345 | T>TG | Donor gain | 0.53 | Moderate |
| 9:21968268 | CGT>C | Acceptor loss | 0.46 | Moderate |
| 9:21968270 | T>TC | Acceptor loss | 0.59 | Moderate |
| 9:21968201 | G>T | Acceptor loss | 0.20 | Low |
| 9:21967787 | G>GC | Acceptor gain | 0.36 | Moderate |
| 9:21968385 | C>T | Donor gain | 0.46 | Moderate |
