## Protein interactions & networks

### Protein-Protein Interactions (PPIs)

**Total Interaction Count (Approximate):**
- STRING: ~8,154 interactions
- BioGRID: 232 interactions
- IntAct: 143 interactions
- **Total: ~8,500+ interactions**

**TOP 30 Highest-Confidence Interacting Proteins (STRING Database):**

| Rank | Gene | Protein Name | STRING Score |
|------|------|--------------|--------------|
| 1 | XIAP | E3 ubiquitin-protein ligase XIAP | 998 |
| 2 | CASP3 | Caspase-3 (self-interaction) | 998 |
| 3 | APAF1 | Apoptotic protease-activating factor 1 | 986 |
| 4 | CASP8 | Caspase-8 | 984 |
| 5 | BCL2 | Apoptosis regulator Bcl-2 | 980 |
| 6 | BAX | Bcl-2-like protein 1 | 975 |
| 7 | CASP9 | Caspase-9 | 958 |
| 8 | PARP1 | Poly [ADP-ribose] polymerase 1 | 954 |
| 9 | MCL1 | Induced myeloid leukemia cell differentiation protein Mcl-1 | 940 |
| 10 | AKT1 | RAC-alpha serine/threonine-protein kinase | 923 |
| 11 | MAPK8 | Mitogen-activated protein kinase 8 | 923 |
| 12 | BIRC3 | Baculoviral IAP repeat-containing protein 3 | 922 |
| 13 | FADD | FAS-associated death domain protein | 920 |
| 14 | ANXA5 | Annexin A5 | 919 |
| 15 | PXDNL | Probable oxidoreductase PXDNL | 918 |
| 16 | PXDN | Peroxidasin homolog | 918 |
| 17 | DFF45 | DNA fragmentation factor subunit beta | 914 |
| 18 | HSPA4 | Heat shock 70 kDa protein 4 | 903 |
| 19 | IL1B | Interleukin-1 beta | 899 |
| 20 | BECN1 | Beclin-1 | 899 |
| 21 | ACTB | Actin, cytoplasmic 1 | 894 |
| 22 | GAPDH | Glyceraldehyde-3-phosphate dehydrogenase | 890 |
| 23 | FASLG | Tumor necrosis factor ligand superfamily member 6 | 888 |
| 24 | TDT | DNA nucleotidylexotransferase | 883 |
| 25 | MYC | Myc proto-oncogene protein | 881 |
| 26 | IL6 | Interleukin-6 | 874 |
| 27 | HSP90AA1 | Heat shock protein HSP 90-alpha | 874 |
| 28 | SMAC | Diablo IAP-binding mitochondrial protein | 864 |
| 29 | PTGS2 | Prostaglandin G/H synthase 2 | 863 |
| 30 | TNF | Tumor necrosis factor | 899 |

**Notable High-Confidence Interactions (IntAct & BioGRID):**
- **XIAP**: 0.870 confidence (multiple independent interactions) - inhibitor of caspase-3
- **TXN (Thioredoxin)**: 0.750 confidence - protein cleavage substrate
- **APAF1**: 0.650 confidence - upstream activator in apoptosome
- **CASP9**: 0.640 confidence - upstream initiator caspase
- **PARP1**: 0.570-0.620 confidence - key apoptotic substrate
- **CASP6**: 0.570 confidence - downstream effector caspase
- **BID**: 0.440 confidence - BH3-only substrate for cleavage
- **BIRC2/BIRC7**: 0.400 confidence - IAP inhibitors
- **MDM2**: 0.440-0.570 confidence - cleavage substrate

---

### Protein Similarity

**Structural/Embedding Similarity (ESM2 - Top 20):**
All 23 ESM2-similar proteins are caspase homologs from various species. Most are non-human orthologs. Human caspases showing high similarity:

| Rank | UniProt | Gene | Top Similarity Score |
|------|---------|------|----------------------|
| 1 | P55210 | CASP7 | 0.9973 |
| 2 | P55212 | CASP6 | 0.9957 |

Other homologs represent orthologs from mouse, rat, zebrafish, and other organisms (O08738, O35397, Q60431, etc.) with similarly high ESM2 scores (0.994-0.998).

**Sequence Homology (DIAMOND - Top 20 by Identity):**

| Rank | UniProt | Protein/Gene | Identity % | BitScore |
|------|---------|--------------|-----------|----------|
| 1 | O75601 | CASP (human variant) | 100.00 | 755.00 |
| 2 | Q5E9C1 | CASP (human variant) | 100.00 | 755.00 |
| 3 | P29594 | CASP (variant) | 97.60 | 897.00 |
| 4 | P55215 | CASP (variant) | 97.60 | 896.00 |
| 5 | G5ECW5 | CASP3-like | 94.20 | 230.00 |
| 6 | Q9TZP5 | CASP3-like | 94.20 | 232.00 |
| 7 | O08738 | Casp3 (mouse) | 94.60 | 542.00 |
| 8 | O35397 | Casp3 (mouse) | 94.60 | 543.00 |
| 9 | Q2PFV2 | Casp3 (ortholog) | 96.00 | 535.00 |
| 10 | Q5IS99 | Casp3 (ortholog) | 96.00 | 531.00 |
| 11 | Q60431 | Casp3 (mouse) | 92.10 | 523.00 |
| 12 | P42573 | CASP2 (human) | 84.50 | 843.00 |
| 13 | P55213 | CASP3 variant | 93.50 | 534.00 |
| 14 | P70677 | CASP3 variant | 93.50 | 537.00 |
| 15 | P55214 | CASP3 variant | 90.80 | 567.00 |
| 16 | P97864 | CASP3 variant | 90.80 | 570.00 |
| 17 | Q3T0P5 | CASP3-related | 88.30 | 517.00 |
| 18 | Q08DY9 | CASP3-related | 87.30 | 488.00 |
| 19 | O08736 | Casp3 (species ortholog) | 87.90 | 744.00 |
| 20 | O89110 | CASP-related | 77.40 | 737.00 |

**Key Interaction Networks:**

1. **Apoptotic Cascade**: CASP3 is activated by initiator caspases (CASP8, CASP9) via APAF1-mediated apoptosome formation and inhibited by XIAP/BIRC proteins
2. **Substrate Specificity**: Cleaves key apoptotic substrates including PARP1, CASP6, BID, DFFA, and GSDME
3. **Cross-caspase Activation**: Directly activates other executioner caspases (CASP6, CASP7)
4. **Apoptotic Signaling**: Interfaces with both extrinsic (FAS/FADD/CASP8) and intrinsic (mitochondrial/CASP9) pathways
5. **Regulatory Interactions**: TXN provides redox regulation; IAPs (XIAP, BIRC2/3/7) provide inhibitory control
