## Protein interactions & networks

### Total Interaction Counts
- **STRING interactions**: 1,626
- **BioGRID interactions**: 1,424
- **IntAct interactions**: 76

---

### TOP 30 Highest-Confidence STRING Interacting Proteins
(Confidence scores on 0-1000 scale; higher = stronger evidence)

| Rank | UniProt ID | Protein Name | Score | Evidence Types |
|------|-----------|--------------|-------|-----------------|
| 1 | Q8TEV9 | Guanine nucleotide exchange factor C9orf72 homolog | 996 | High confidence |
| 2 | Q9HAD4 | C9orf72-related protein | 995 | High confidence |
| 3 | Q13148 | Ras-related protein Ral-A | 944 | Predicted, database, experimental |
| 4 | A0A087WTZ4 | ULK1-associated protein | 919 | Predicted, database |
| 5 | P35637 | Serine/threonine-protein kinase RAF1 | 919 | Multiple evidence |
| 6 | P00441 | Phosphatidylinositol 3-kinase catalytic subunit alpha | 898 | Predicted, experimental |
| 7 | P23781 | Mitogen-activated protein kinase 1 (ERK2) | 888 | Predicted, database |
| 8 | Q9UHD9 | Autophagy-related protein 5 (ATG5) | 871 | Multiple evidence |
| 9 | P11476 | Mitochondrial 28S ribosomal protein S36 | 863 | Predicted |
| 10 | Q99700 | SH3 domain-containing GTPase-activating protein 1 | 857 | Predicted, experimental |
| 11 | P10636 | Microtubule-associated protein tau (MAPT) | 852 | Predicted, database |
| 12 | P55072 | Transitional endoplasmic reticulum ATPase (VCP) | 852 | Predicted, experimental |
| 13 | Q9UQN3 | Autophagy-related protein 13 (ATG13) | 826 | Predicted, experimental |
| 14 | Q9UHD2 | ULK1 serine/threonine-protein kinase | 824 | Multiple evidence |
| 15 | Q13501 | Ras-related protein Rab-6A | 821 | Predicted, database |
| 16 | Q96Q42 | Kinesin-like protein KIF20A | 819 | Predicted |
| 17 | P09651 | Heat shock protein HSP90-alpha | 809 | Predicted, experimental |
| 18 | O14966 | Ras-related protein Rab-7A | 786 | Predicted, experimental |
| 19 | Q96CV9 | Vesicle-associated membrane protein 7 (VAMP7) | 774 | Predicted |
| 20 | P31943 | Ras-related protein Rab-8A | 754 | Predicted, experimental |
| 21 | Q9NUM4 | Syntaxin-17 | 754 | Predicted, experimental |
| 22 | P51991 | Ras-related protein Rab-11A | 750 | Predicted |
| 23 | P05067 | Amyloid beta A4 protein (APP) | 749 | Predicted, database |
| 24 | O95292 | Dynactin subunit 1 (DCTN1) | 746 | Predicted, experimental |
| 25 | Q7Z333 | Sequestosome-1 (p62/SQSTM1) | 741 | Predicted, experimental |
| 26 | P55795 | Serine/threonine-protein kinase TBK1 | 724 | Predicted, experimental |
| 27 | Q8WYQ3 | Ras-related protein Rab-39B | 721 | Predicted, experimental |
| 28 | P43243 | Matrin 3 (MATR3) | 715 | Predicted |
| 29 | Q8WXG6 | WD repeat protein 41 (WDR41) | 692 | Multiple evidence |
| 30 | Q14203 | TAR DNA-binding protein 43 (TARDBP) | 688 | Predicted, experimental |

**Key Findings**: Top interactors are primarily involved in autophagy (ATG5, ATG13, ULK1), protein degradation (HSP90, VCP), vesicular transport (Rab proteins, VAMP7), and cytoskeletal dynamics (tau, dynactin). Strong interaction with SMCR8 (binding partner) and proteins in autophagy-lysosomal pathways.

---

### TOP 20 Proteins with Highest Structural/Embedding Similarity (ESM2)
(ESM2: AlphaFold/language-model based structural similarity; scale 0-1)

| Rank | UniProt ID | Similarity Score | Avg Similarity |
|------|-----------|------------------|-----------------|
| 1 | Q5RC62 | 1.0000 | 0.9753 |
| 2 | Q5RD58 | 1.0000 | 0.9876 |
| 3 | Q66HC3 | 1.0000 | 0.9867 |
| 4 | Q6DFW0 | 1.0000 | 0.9867 |
| 5 | Q6NSW5 | 0.9999 | 0.9865 |
| 6 | Q6ZW61 | 1.0000 | 0.9752 |
| 7 | Q86WG5 | 0.9999 | 0.9885 |
| 8 | Q8R3P6 | 0.9999 | 0.9883 |
| 9 | Q8TCE6 | 0.9999 | 0.9865 |
| 10 | Q8WVF5 | 0.9998 | 0.9821 |
| 11 | Q9CZW2 | 0.9989 | 0.9846 |
| 12 | Q9D7X1 | 0.9997 | 0.9831 |
| 13 | Q9D8N2 | 0.9997 | 0.9866 |
| 14 | Q9NQ89 | 1.0000 | 0.9876 |
| 15 | Q96SY0 | 0.9999 | 0.9884 |
| 16 | A6H6X4 | 0.9998 | 0.9827 |
| 17 | D4A770 | 0.9999 | 0.9871 |
| 18 | E9PXF8 | 0.9999 | 0.9884 |
| 19 | Q1T765 | 0.9982 | 0.9828 |
| 20 | Q28HN9 | 0.9952 | 0.9856 |

**Profile**: All 55 ESM2-similar proteins share highly conserved structural folds with C9orf72 (>0.97 average similarity). Many are orthologs or paralogs across species (indicated by prefix Q/A/D/E identifiers for organism origin).

---

### TOP BioGRID Interacting Proteins with Experimental Evidence

| Gene Symbol | Evidence Type | Count |
|-------------|---------------|-------|
| SMCR8 | Affinity Capture-MS, Affinity Capture-Western, Reconstituted Complex, Two-hybrid | High (core complex member) |
| WDR41 | Affinity Capture-MS, Affinity Capture-Western | High (core complex member) |
| ULK1 | Affinity Capture-Western, Two-hybrid, Co-localization, Phosphorylation | High (autophagy pathway) |
| ATG13 | Affinity Capture-Western, Two-hybrid, Co-localization, Phosphorylation | High (autophagy pathway) |
| ATG101 | Affinity Capture-MS, Affinity Capture-Western | High (autophagy complex) |
| RB1CC1 | Affinity Capture-MS, Affinity Capture-Western | High (FIP200, autophagy) |
| RAB8A | Affinity Capture-Western, GEF reaction | Medium (GTPase substrate) |
| RAB39B | Affinity Capture-Western, Affinity Capture-MS | Medium (GTPase substrate) |
| SETX | Affinity Capture-MS | Medium (ALS-associated) |
| DCTN1 | Affinity Capture-MS | Medium (dynactin complex) |
| TBK1 | Affinity Capture-MS, Phosphorylation | Medium (immune response) |
| TARDBP | Affinity Capture-MS | Medium (ALS-associated) |
| SQSTM1 | Affinity Capture-MS | Low-Medium (autophagy adaptor) |
| UBIQUITIN | Implied through E3 ligases | Multiple |

---

### Sequence Homology: Orthologs (Cross-Species)

| Species | Gene ID | Gene Symbol | Identity Details |
|---------|---------|-------------|-------------------|
| *Homo sapiens* | ENSG00000147894 | C9orf72 | Reference (human) |
| *Mus musculus* | ENSMUSG00000028300 | C9orf72 | High orthology (mammalian) |
| *Rattus norvegicus* | ENSRNOG00000009478 | RGD1359108 | High orthology (mammalian) |
| *Danio rerio* | ENSDARG00000011837 | C13H9orf72 | Conserved in vertebrates |

**Sequence Conservation**: Orthologs show strong amino acid conservation across mammals, with structural domains preserved in all vertebrate orthologs examined.
