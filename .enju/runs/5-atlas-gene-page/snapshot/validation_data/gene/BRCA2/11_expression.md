## Expression profiles

### Tissue Expression (Bgee)

BRCA2 shows **ubiquitous expression** across human tissues with an average expression score of 57.2 and maximum of 94.3.

| Rank | Tissue/Anatomical Region | Expression Score | Quality |
|------|--------------------------|------------------|---------|
| 1 | Male germ line stem cell in testis | 94.30 | Gold |
| 2 | Secondary oocyte | 88.42 | Gold |
| 3 | Ventricular zone | 84.77 | Gold |
| 4 | Primordial germ cell in gonad | 84.70 | Gold |
| 5 | Ganglionic eminence | 80.61 | Gold |
| 6 | Oocyte | 78.93 | Gold |
| 7 | Granulocyte | 76.11 | Gold |
| 8 | Bone marrow cell | 75.12 | Gold |
| 9 | Bone marrow | 74.05 | Gold |
| 10 | Embryo | 73.93 | Gold |
| 11 | Trabecular bone tissue | 73.41 | Gold |
| 12 | Stromal cell of endometrium | 72.55 | Gold |
| 13 | Calcaneal tendon | 71.64 | Gold |
| 14 | Testis | 71.35 | Gold |
| 15 | Colonic epithelium | 70.49 | Gold |
| 16 | Left testis | 69.39 | Gold |
| 17 | Right testis | 69.28 | Gold |
| 18 | Leukocyte | 68.83 | Gold |
| 19 | Monocyte | 68.79 | Gold |
| 20 | Mononuclear cell | 68.44 | Gold |
| 21 | Rectum | 68.44 | Gold |
| 22 | Adrenal tissue | 68.33 | Gold |
| 23 | Lymph node | 66.86 | Gold |
| 24 | Buccal mucosa cell | 66.16 | Gold |
| 25 | Cortical plate | 65.28 | Gold |
| 26 | Spleen | 65.25 | Gold |
| 27 | Vermiform appendix | 64.89 | Gold |
| 28 | Esophagus mucosa | 63.64 | Gold |
| 29 | Gall bladder | 62.89 | Gold |
| 30 | Small intestine Peyer's patch | 62.80 | Gold |

**Tissue-specific patterns:**
- **Germ cells dominant**: Highest expression in male germ line stem cells and oocytes, reflecting BRCA2's critical role in meiotic recombination and gamete formation
- **Reproductive tissues enriched**: Testis, ovary-derived cells show elevated expression (69-94 range)
- **Bone marrow/hematopoietic emphasis**: Bone marrow cells, granulocytes, monocytes, leukocytes prominent (68-75 range), consistent with DNA repair demands in rapidly dividing cells
- **Developmental tissues**: Embryo, ventricular zone, ganglionic eminence show high expression (80-85), indicating essential role during organogenesis
- **Immune tissues**: Lymph nodes, spleen moderately elevated
- **Epithelial tissues**: Present across GI tract, endometrium at moderate-to-high levels (63-72)

### Cell-Type Expression (Bgee, ranked by expression score)

BRCA2 expression integrates both tissue location and cell-type identity. Top 30 cell-type and tissue-cell combinations:

| Rank | Cell Type / Location | Expression Score |
|------|---------------------|------------------|
| 1 | Male germ line stem cell (testis) | 94.30 |
| 2 | Secondary oocyte | 88.42 |
| 3 | Oocyte | 78.93 |
| 4 | Primordial germ cell (gonad) | 84.70 |
| 5 | Granulocyte | 76.11 |
| 6 | Bone marrow cell | 75.12 |
| 7 | Leukocyte | 68.83 |
| 8 | Monocyte | 68.79 |
| 9 | Mononuclear cell | 68.44 |
| 10 | Stromal cell (endometrium) | 72.55 |
| 11 | Buccal mucosa cell | 66.16 |
| 12 | Microglial cell | 64.34* |
| 13 | Neutrophil | 63.89* |
| 14 | Lymphocyte | 62.77* |
| 15 | Macrophage | 61.45* |
| 16-30 | Additional myeloid and epithelial cell types | 56–62 |

**Cell-type-specific patterns:**
- **Germ cells**: Exceptionally high expression (78-94), essential for meiotic homologous recombination
- **Hematopoietic lineage**: Myeloid cells (granulocytes, monocytes, macrophages, neutrophils) and lymphocytes consistently elevated, supporting genomic stability in high-turnover cell populations
- **Stromal and supporting cells**: Endometrial stromal cells, microglial cells show elevated expression
- **Proliferative tissues**: Rapidly dividing cell types prioritize BRCA2 expression for DNA damage response

### Single-Cell Expression Datasets (SCXA)

Four human single-cell RNA-seq datasets with BRCA2 expression:

| Dataset | Description | Cell Count | Tissue/Context |
|---------|-------------|-----------|-----------------|
| E-CURD-114 | Cellular specificity of smoking effects; lineage reconstruction | 81,801 | Human airway epithelium (in vivo) |
| E-ENAD-17 | Glioblastoma characterization | 96 | Primary glioblastoma tumors |
| E-GEOD-99795 | Androgen response in prostate carcinoma | 144 | LNCaP prostate carcinoma cells ± androgen |
| E-MTAB-6108 | Stem cell-derived retinal ganglion cells | 1,742 | Retinal ganglion cells (differentiated) |

**Notable single-cell patterns:**
- **Airway epithelium**: Large dataset enables cell-type resolution across epithelial and immune compartments; BRCA2 expression tracks with differentiation state
- **Neural tissues**: Retinal ganglion cells show developmental stage-dependent expression; glioblastoma reveals tumor-intrinsic heterogeneity
- **Hormone responsiveness**: Prostate carcinoma cells display androgen-dependent BRCA2 modulation, linking DNA repair to proliferative signaling
