## Drug & pharmacology data

### EGFR as Drug Target
EGFR (Epidermal growth factor receptor) is a **well-established and extensively targeted** protein in drug development. Over **10,000 molecules** have been tested against EGFR in ChEMBL, with multiple approved drugs on the market.

### Approved EGFR-Targeting Molecules (Top 5, Phase 4)

| Molecule | ID | Mechanism | Status | Clinical Trials |
|----------|----|-----------|---------|----|
| **Osimertinib** | CHEMBL3353410 | EGFR TKI (irreversible, 3rd generation) | Phase 4 | 228 trials |
| **Erlotinib** | CHEMBL553 | EGFR TKI (reversible, 1st generation) | Phase 4 | 496 trials |
| **Gefitinib** | CHEMBL939 | EGFR TKI (reversible, 1st generation) | Phase 4 | 294 trials |
| **Afatinib** | CHEMBL1173655 | EGFR/HER2 TKI (irreversible, 2nd generation) | Phase 4 | 179 trials |
| **Lapatinib** | CHEMBL554 | EGFR/HER2 dual TKI | Phase 4 | 261 trials |

All five approved drugs are small-molecule tyrosine kinase inhibitors. Additional ~9,995 molecules in earlier development phases.

### Clinical Trials (Selected Top 20 by Drug)

**Erlotinib (CHEMBL553)** — Primary indications: NSCLC, pancreatic cancer, head/neck cancer
- NCT01287754: NSCLC with EGFR mutations | Phase 4 | COMPLETED
- NCT01609543: 1st-line lung adenocarcinoma with EGFR mutations | Phase 4 | COMPLETED  
- NCT00446225: NSCLC with EGFR TK domain mutations | Phase 3 | COMPLETED
- NCT01024413: Erlotinib vs Gefitinib in advanced NSCLC with EGFR exon 19/21 mutations | Phase 3 | COMPLETED
- NCT00349219 (TORCH): Erlotinib vs chemotherapy for advanced NSCLC | Phase 3 | COMPLETED
- NCT02296125: Osimertinib (AZD9291) vs Gefitinib/Erlotinib in NSCLC | Phase 3 | COMPLETED
- NCT02411448 (RELAY): Ramucirumab + Erlotinib for EGFR-mutant NSCLC | Phase 3 | ACTIVE

**Gefitinib (CHEMBL939)** — Primary indications: NSCLC (adenocarcinoma, especially never-smokers)
- NCT00076388 (IRESSA vs Docetaxel): Gefitinib vs chemotherapy | Phase 3 | COMPLETED
- NCT00322452 (IPASS): 1st-line Gefitinib vs carboplatin/paclitaxel in Asia | Phase 3 | COMPLETED
- NCT01774721 (ARCHER1050): Dacomitinib vs Gefitinib for 1st-line NSCLC | Phase 3 | COMPLETED
- NCT01404260: Gefitinib intercalating with chemotherapy | Phase 3 | COMPLETED
- NCT02296125: Osimertinib vs Gefitinib/Erlotinib | Phase 3 | COMPLETED
- NCT02588261: ASP8273 vs Erlotinib/Gefitinib with EGFR mutations | Phase 3 | TERMINATED

**Afatinib (CHEMBL1173655)** — Primary indications: NSCLC, head/neck cancer, squamous cell carcinoma
- NCT00949650 (LUX-Lung 3): Afatinib 1st-line vs chemotherapy in EGFR-mutant NSCLC | Phase 3 | COMPLETED
- NCT01121393 (LUX-Lung 4): Afatinib vs gemcitabine/cisplatin | Phase 3 | COMPLETED
- NCT01523587 (LUX-Lung 8): Afatinib vs Erlotinib in squamous NSCLC | Phase 3 | COMPLETED
- NCT01466660 (LUX-Lung 7): Afatinib vs Gefitinib for 1st-line EGFR-mutant adenocarcinoma | Phase 2 | COMPLETED

**Osimertinib (CHEMBL3353410)** — Primary indication: NSCLC (especially T790M resistance mutations)
- NCT02296125: Osimertinib vs Gefitinib/Erlotinib in NSCLC | Phase 3 | COMPLETED
- NCT02151981 (AURA): Osimertinib in EGFR-mutant NSCLC with acquired T790M | Phase 2 | COMPLETED

### Pharmacogenomics & Drug Response Predictors

**Key EGFR Mutations Affecting Drug Response:**

| Mutation | Clinical Relevance | Drug Sensitivity |
|----------|-------------------|------------------|
| **Exon 19 deletion** | ~45% of EGFR+ NSCLC | Sensitive to all 1st/2nd gen TKIs (erlotinib, gefitinib, afatinib) |
| **L858R (exon 21 point mutation)** | ~40% of EGFR+ NSCLC | Sensitive to 1st/2nd gen TKIs |
| **T790M** | Acquired resistance mechanism (~50% after progression on 1st/2nd gen TKI) | Resistant to 1st/2nd gen TKIs; **sensitive to osimertinib** (3rd gen) |
| **Exon 20 insertion** | ~5% of EGFR mutations | Generally resistant to 1st/2nd gen TKIs; variable osimertinib response |
| **G719X** | ~5% of mutations | Intermediate sensitivity |

**Dosing Considerations (from approved labels):**
- **Erlotinib**: 150 mg daily (150 mg daily oral); adjust for CYP3A4 inducers/inhibitors
- **Gefitinib**: 250 mg daily oral
- **Afatinib**: Dose reduced from 40-50 mg if diarrhea occurs (most common limiting toxicity)
- **Osimertinib**: 80 mg daily (adjusted to 40 mg if not tolerated); superior CNS penetration
- **Lapatinib**: 1250 mg daily (with capecitabine in HER2+ breast cancer)

**No major pharmacogenomic variant panels** (e.g., DPYD, NAT2) are standard for EGFR TKIs, but **EGFR genotyping is mandatory** to guide drug selection. T790M testing at progression determines osimertinib eligibility.
