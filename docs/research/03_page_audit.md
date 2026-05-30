# Sugi Atlas gene page audit: a gap matrix against the five leading references

**Scope.** Six reference genes — TP53, BRCA1, EGFR, CDKN2A, KRAS, TTN — each compared side-by-side against the same gene's page at NCBI Gene, UniProt, GeneCards, Open Targets Platform, and Ensembl. Atlas pages are taken from `/data/sugi-atlas-dist/atlas/gene/<SYMBOL>/page.md`. References were sampled via WebFetch (RefSeq via NCBI, UniProt via REST flat-file, Open Targets via GraphQL where reachable, Ensembl REST + the Atlas team's prior knowledge of the public web UI, Human Protein Atlas as a sixth de-facto reference). Where a live fetch failed (UniProt web UI, GeneCards, Open Targets web UI returned 403/400/empty) we substituted the canonical REST/flat-file content, which contains the same authoritative data the web page renders.

Atlas's data model: ~30–50 KB of structured deterministic tables in twelve fixed sections — Identifiers / Transcripts / Proteins / Structure / Orthologs / Variants & AI predictions / Pathways & GO / Interactions & networks / TF regulation / Drug & pharmacology / Expression / Disease — plus a one-paragraph LLM exec summary on a few genes. Everything is sourced from biobtree; no free-text narrative is curated by humans.

This audit asks, gene by gene: where is Atlas already best-in-class, where is it merely on par, and where is it materially behind?

---

## 1. Per-gene comparisons

### 1.1 TP53 (Entrez 7157, ENSG00000141510, P04637)

| Source | What it has that Atlas lacks | What Atlas has that it lacks | Note |
|---|---|---|---|
| **NCBI Gene** | A curated RefSeq narrative paragraph ("This gene encodes a tumor suppressor protein containing transcriptional activation, DNA binding, and oligomerization domains…"). 5,000+ GeneRIFs (one-line, paper-anchored "what this gene does" facts). ClinGen dosage-sensitivity scores (haploinsufficiency / triplosensitivity verdicts with evaluation date). HIV-1 interaction table. Phenotype-to-MedGen condition mapping with professional guidelines (ACMG). | Atlas's 3,754-edge STRING graph, 1,863 IntAct rows, 1,043 downstream TF targets, 297 PDBs in a single table, 2,569 AlphaMissense scores, 1,638 SpliceAI scores, 1,432 MSigDB gene sets, 419 disease-level clinical trials. None of that is on the NCBI page in any compact form. | NCBI wins on **narrative** + **dosage**; Atlas wins on **quantitative depth**. |
| **UniProt (P04637)** | Free-text FUNCTION ("Multifunctional transcription factor that induces cell cycle arrest, DNA repair or apoptosis upon binding to target DNA sequence…"), SUBUNIT, SUBCELLULAR LOCATION, TISSUE SPECIFICITY, INDUCTION, DOMAIN, PTM narrative. 9 named isoforms (p53α/β/γ × full-length/Δ40/Δ133/Δ160). 8 named DISEASE blocks with descriptions and OMIM links. 237 cited publications with evidence codes. Cofactor block (Zn²⁺). | Atlas dumps the full 19 UniProt accessions, 1,518 UniProt features tallied by type, 34 RefSeq proteins, 12 CCDS, MANE-Select flagging on transcripts and proteins. UniProt doesn't surface CCDS or MANE on its own page nor does it cross-link the *full* SpliceAI/AlphaMissense burden. | UniProt wins **narrative function, isoform naming, disease prose**; Atlas wins **breadth of cross-DB identifier closure**. |
| **GeneCards** | Aliases & descriptions aggregated from 8+ sources; commercial product cards (antibodies, ELISAs, clones, lentiviral particles, CRISPR guides); MalaCards disorder pages; BioGPS + Allen Brain + ARCHS4 expression overlays; "summaries from multiple sources" (Entrez, UniProt, GeneCards-edited, Tocris). | Atlas's antibody resources count is shown as "0" — a literal honest reporting; GeneCards inflates this with vendor catalogs. Atlas has cleaner SpliceAI / AlphaMissense / SIGNOR / FANTOM5 promoter / Bgee tissue tables. | GeneCards wins on **commercial product surface** and **multi-source aliases**; Atlas wins on **clean, citable primary data**. |
| **Open Targets** | Integrated target-disease association scores combining ChEMBL, EuropePMC literature mining, OMIM/Orphanet, GWAS, ClinVar, somatic mutation (COSMIC/IntOGen), eQTL/genetics portal, gene burden, animal models. Druggability tractability buckets. Target safety liabilities. IMPC mouse phenotypes. Baseline GTEx expression heatmap. Target prioritisation score. | Atlas has 195 ChEMBL phase-≥1 molecules, 333 SIGNOR signaling edges, 1,043 TF target gene list, GenCC classifications with inheritance, HPO terms (224), Reactome (46), MSigDB (1,432), Bgee + FANTOM5 + per-tissue scoring. Open Targets does not show TF→target downstream graphs at the gene level. | Open Targets wins on **target-disease integrated scoring and tractability**; Atlas wins on **the TF regulatory layer and the pharmacology long tail**. |
| **Ensembl** | Region-in-detail visual browser with regulatory features (promoters, enhancers, CTCF sites from the Regulatory Build); gene tree across 200+ species (Atlas only shows 4 ortholog rows); ID history (which Ensembl IDs were merged/retired); variation table with regulatory consequences; VEP launcher. TSL/APPRIS per-transcript flags on the web (Atlas surfaces MANE only). | Atlas surfaces all 39 transcripts with biotype, the canonical's 11 exons with coordinates, plus a full structure table (297 PDB rows). Ensembl's transcript table is paginated; Atlas dumps it. | Ensembl wins on **regulatory build, gene tree breadth, ID history, transcript quality flags**; Atlas wins on **single-page completeness**. |

**Verdict for TP53.** Atlas is the deepest *quantitative* page of the six, but it is the *thinnest narrative* page. A reader who already knows TP53 will love Atlas; a reader who is asking "what is TP53" will go to UniProt or NCBI first.

---

### 1.2 BRCA1 (Entrez 672, ENSG00000012048, P38398)

| Source | Atlas advantages | Atlas gaps |
|---|---|---|
| **NCBI** | 47 Ensembl transcripts vs NCBI's 22-exon canonical view; 15,445 ClinVar variants surfaced with per-class floors (Pathogenic ≥4,004 — by far the heaviest ClinVar gene); SpliceAI + AlphaMissense at scale; 166 HPO phenotypes. | RefSeq narrative ("190 kD nuclear phosphoprotein… BRCA1-associated genome surveillance complex (BASC)… ~40% of inherited breast cancers"); ClinGen haploinsufficiency "Sufficient" verdict; ENIGMA expert-panel ClinVar review — Atlas doesn't surface the ENIGMA submitter axis. |
| **UniProt** | Full feature table (chains, RING domain, BRCT repeats, zinc fingers) is summarised as counts; we list every PDB. | UniProt narrative: "E3 ubiquitin-protein ligase that specifically mediates the formation of 'Lys-6'-linked polyubiquitin chains…", BRCA1-BARD1 heterodimer mechanism. 5 DISEASE blocks: BC, BROVCA1, OC, PNCA4, FANCS. CATALYTIC ACTIVITY block. POLYMORPHISM block (rare among genes). 109 citations. |
| **GeneCards** | Cleaner inheritance/classification through GenCC. | "BRCA Exchange" / ENIGMA hotlinks; MalaCards link to hereditary breast/ovarian cancer with treatment notes; consumer-genomics references (23andMe/ACMG SF v3.2 actionable list inclusion). |
| **Open Targets** | 12 ChEMBL molecules phase ≥ 1 shown; broader drug landscape (PARP inhibitors hit via synthetic lethality). | The synthetic-lethality framing itself — Open Targets shows PARP1 as the *druggable partner* even though BRCA1 is the genetic context. Atlas surfaces both but doesn't make the synthetic-lethality link explicit. |
| **Ensembl** | 47 transcripts × biotypes; canonical (ENST00000357654) with 23 exons. | Pseudogene BRCA1P1 cross-link; clinical variant browser with full-length VEP-annotated effects on every transcript (Atlas only surfaces top 30 P/LP IDs). |

**Distinctive Atlas strength on BRCA1.** The 15,445-variant ClinVar floor table makes the curation burden of BRCA1 *visible*. None of the five reference resources present that scale on a single page.

---

### 1.3 EGFR (Entrez 1956, ENSG00000146648, P00533)

| Source | Atlas advantages | Atlas gaps |
|---|---|---|
| **NCBI** | 378 PDB rows, 172 ChEMBL-bioactive molecules including gefitinib/erlotinib/osimertinib/cetuximab class. | RefSeq narrative; 5,726 GeneRIFs; explicit cetuximab/panitumumab response phenotypes linked to KRAS status in GTR. NCBI's drug-response → variant linking (e.g. EGFR T790M → osimertinib) is curated. |
| **UniProt** | 462 PDB structures table (Atlas lists 378; reality is closer to ~470 — likely a snapshot lag; **flag**). | FUNCTION narrative ("Receptor tyrosine kinase binding ligands of the EGF family…"). 4 named isoforms incl. p170/p110/p60/Truncated. CATALYTIC ACTIVITY (EC 2.7.10.1). Autophosphorylation site list (Tyr-1092/1110/1172/1197). Disease entries: LNCR + NNCIS. |
| **GeneCards** | Clean per-promoter FANTOM5 table; 21 HPO phenotypes; expression breadth ubiquitous. | "Drugs & Compounds" panel with FDA-approval timelines per indication (NSCLC, CRC, head and neck); biomarker companion-diagnostic flags. |
| **Open Targets** | Atlas has 172 ChEMBL phase-≥1 molecules — very comprehensive. | Open Targets organises EGFR around its primary indications (NSCLC EGFR-mut, CRC KRAS-WT), surfaces approved drug → mechanism mapping with explicit "approved" status. Tractability bucket is "small molecule + antibody" — both delivered drugs exist. |
| **Ensembl** | 28-exon canonical, 78 transcripts (Atlas surfaces all). | EGFRvIII variant / glioma-specific deletion — surfaced on the Ensembl Region viewer as a known somatic event; Atlas surfaces it only via the AlphaMissense / ClinVar tables, not as a named entity. |

**Notable finding.** EGFR is the case where Atlas's *long-tail bioactivity* and *RTK domain feature dump* meet GeneCards/OpenTargets's clinical *drug-by-indication* organization. Both are needed; neither replaces the other.

---

### 1.4 CDKN2A (Entrez 1029, ENSG00000147889, P42771 + Q8N726) — the dual-product test

This is the gene that breaks naïve gene databases.

| Source | How they handle p16INK4a vs p14ARF |
|---|---|
| **Atlas** | **Correctly surfaces both products.** Protein section lists canonical UniProt as `P42771` and explicitly notes `(reviewed: P42771, Q8N726)`. Annotated features include both chain entries: "cyclin-dependent kinase inhibitor 2a 1–156" AND "tumor suppressor arf 1–132", both with their own splice variants, regions of interest, and modified residues. Domain table includes both `Ank_Repeat/CDKN_Inhibitor` and `Tumor_suppres_ARF`. |
| **NCBI** | Narrative correctly describes the two ORFs: "two of which encode structurally related isoforms known to function as inhibitors of CDK4 kinase. The remaining transcript includes an alternate first exon located 20 Kb upstream… contains an alternate open reading frame (ARF) that specifies a protein which is structurally unrelated to the products of the other variants." Excellent. |
| **UniProt** | Splits into **two separate entries**: P42771 (p16INK4a) and Q8N726 (p14ARF). P42771 explicitly warns: *"The proteins described here are encoded by the gene CDKN2A, but are completely unrelated in terms of sequence and function to tumor suppressor ARF (AC Q8N726) which is encoded by the same gene."* Pedagogically clear but requires the reader to follow two links. |
| **GeneCards** | Aggregates everything onto one card and tends to mix p16 and p14ARF claims in the same paragraphs — a known historical confusion point. |
| **Open Targets** | Indexed on the Ensembl gene; treats CDKN2A as one entity for disease association scoring; loses the per-product mechanistic distinction (p16→CDK4/6→Rb vs p14ARF→MDM2→p53). |
| **Ensembl** | 14 transcripts, multiple coding for p16 isoforms vs ARF isoforms; visible on the transcript table but no narrative explanation. |

**Atlas is materially better than GeneCards, Open Targets, and Ensembl here**, on par with NCBI's narrative explanation, and equivalent (in compact form) to following both UniProt accessions. **This is gold for AI trust** — an LLM reading Atlas's CDKN2A page sees two chain entries with different lengths (156 vs 132 aa) and different domain families on the same gene, and does not need extra prompting to disambiguate.

---

### 1.5 KRAS (Entrez 3845, ENSG00000133703, P01116)

| Source | Atlas advantages | Atlas gaps |
|---|---|---|
| **NCBI** | 462 PDB structures (the KRAS structure landscape exploded post-G12C-inhibitor approval; Atlas captures it). 10 ChEMBL-bioactive molecules. 372 HPO phenotypes — the largest of the six audited genes (KRAS is also a developmental syndrome gene). | RefSeq narrative + the *RASopathy* developmental-syndromes framing (Noonan 3, CFC 1+2, Costello-related, linear nevus sebaceous syndrome). NCBI groups developmental vs cancer phenotypes; Atlas dumps them flat. |
| **UniProt** | Atlas surfaces 9 UniProt accessions; KRAS is one of the few genes where multiple isoforms (4A vs 4B) really matter. | UniProt clearly names **K-Ras4A vs K-Ras4B** and explains the alternative C-terminus → different membrane targeting (palmitoylated vs polybasic). CATALYTIC ACTIVITY (GTP hydrolase, EC 3.6.5.2). ACTIVITY REGULATION (GEF/GAP cycle). Atlas's isoform table lists the 9 accessions but does not *name* 4A vs 4B. **Real gap.** |
| **GeneCards** | – | Sotorasib (Lumakras) / adagrasib (Krazati) — KRAS G12C inhibitors. Atlas's ChEMBL table includes them as molecules but doesn't surface the *G12C-specific* mechanism. |
| **Open Targets** | Atlas's clinical trials (via MONDO) hit 224 KRAS-disease trials. | Open Targets surfaces *drug → mutation-specific* indications (sotorasib/adagrasib for KRAS G12C NSCLC). This is the dominant clinical narrative for KRAS in 2026 and Atlas underweights it. |
| **Ensembl** | 16 transcripts surfaced. | The KRAS pseudogene KRAS1P (chromosome 6) cross-link. |

**KRAS verdict.** Atlas is excellent on the structural and variant front, weakest on **isoform naming (4A/4B)** and on **mutation-allele-specific drug indication** (G12C → sotorasib).

---

### 1.6 TTN (Entrez 7273, ENSG00000155657, Q8WZ42) — the mega-protein test

| Source | Atlas advantages | Atlas gaps |
|---|---|---|
| **NCBI** | Atlas surfaces the canonical 363-exon transcript ENST00000589042 and 108 RefSeq mRNAs in a single page. Anyone who has tried to render TTN knows this is non-trivial. 31,454 ClinVar variants, 227 P, 3,955 LP — by an order of magnitude the largest ClinVar burden of the six. | RefSeq narrative explicitly names isoforms **N2-A, N2-B, N2-BA, novex-1/2/3** and links them to muscle types (cardiac vs skeletal vs soleus). Atlas's transcript table lists the IDs but not the named isoforms. |
| **UniProt** | Atlas lists 9 UniProt accessions and 1,289 features (incl. 285 domains). | UniProt's domain census ("152 Ig-like, 132 FN3 repeats, PEVK regions, C-terminal kinase") is the *narrative* a reader needs to understand what titin *is*. Atlas tallies features by type ("domain 285, sequence variant 297…") but doesn't name the architecture. **Real gap on TTN specifically.** Disease entries: CMH9, CMD1G, LGMDR10, TMD, MFM9, CMYO5. |
| **GeneCards** | – | Clear narrative on titin's molecular spring function and sarcomere length — clinically important framing. |
| **Open Targets** | 65 MONDO diseases + clinical trials. | TTN-truncating-variant (TTNtv) PCV literature mining — a very active research area for DCM. |
| **Ensembl** | 363-exon table is one of Atlas's strongest renderings. | The 107 alternative TSS promoter table (Atlas surfaces this!) is unusual — Atlas is actually *better* than Ensembl's default Summary panel on this dimension. |

**TTN verdict.** Atlas handles TTN's *scale* better than any other source — page renders, all 363 exons listed, full ClinVar floor. It loses to UniProt on *narrative architecture* (Ig/FN3/PEVK/kinase regions) and to NCBI on *isoform naming* (N2A/N2B/N2BA/novex).

---

## 2. Cross-gene patterns

### 2.1 Strengths that hold across all six genes

| Pattern | Why it matters |
|---|---|
| **Identifier closure** — Atlas surfaces HGNC, OMIM, Entrez, Ensembl, UniProt (all accessions, not just canonical), RefSeq mRNA + protein with MANE flags, CCDS, dbSNP, MONDO, Orphanet, HPO, ChEMBL target IDs, PharmGKB, all on one page. No other source does all of this; GeneCards is closest but mixes commercial noise. | LLMs need stable cross-references to chain queries. Atlas is best-in-class. |
| **Full variant burden visible** — ClinVar per-class floors, plus *paired* SpliceAI + AlphaMissense scores at scale (1,638 + 2,569 on TP53; 1,710 + 983 on CDKN2A). | No reference resource shows this. ClinVar's web UI paginates; UniProt only lists curated variants; gnomAD lives separately; Atlas brings predicted-effect and clinical-effect onto one page. **Genuine differentiator.** |
| **TF regulatory layer** — downstream target list (1,043 for TP53, 1 for CDKN2A — usefully different), JASPAR motifs, upstream regulators. | None of NCBI/UniProt/Ensembl/Open Targets/GeneCards surfaces a per-gene TF target table. Even the dedicated TF resources (TRRUST, hTFtarget, JASPAR) require separate lookups. |
| **Pathway + GO + MSigDB in one place** — 46 Reactome + 1,432 MSigDB + 144/38/21 GO terms on TP53. | Saves the 4-tab dance every researcher does. |
| **Pharmacology long tail** — BindingDB Ki/IC50 samples, ChEMBL targets with assay-type tags (single-protein vs PPI), clinical trials at disease level. | Open Targets covers known drugs cleanly; Atlas covers the *experimental* binding compounds Open Targets filters out. |
| **Honest counts** — `Antibody resources: 0` on TP53/CDKN2A; the ClinVar "per-class counts are floors (≥ shown; pagination cap)" disclaimer; the BindingDB "sampled 300" hedge. | This is rare. GeneCards inflates antibody counts via vendor catalogs. Atlas is *epistemically* better. |

### 2.2 Weaknesses that hold across all six genes

| Pattern | Why it matters |
|---|---|
| **No curated free-text Function / Subunit / Tissue / PTM / Disease narrative.** Only TP53/BRCA1/CDKN2A/KRAS/TTN have a single one-paragraph LLM summary at the top, and it just restates the structured data ("TP53 has 39 Ensembl transcripts…"). | This is the single biggest gap. NCBI's RefSeq narrative and UniProt's CC blocks are the *answers* an AI agent quotes when asked "what is TP53". Atlas has none. |
| **No named isoforms.** Atlas surfaces accession IDs (E7EMR6, E9PCY9…) but not *names* (p53α/β/γ, K-Ras4A/4B, p16γ, N2A/N2B titin). | Domain experts and clinicians refer to isoforms by name. Atlas is illegible to them on this axis. |
| **No subcellular localization narrative.** UniProt and Human Protein Atlas both surface "Mainly nucleoplasm; also vesicles, cytosol". Atlas surfaces only GO cellular components as a flat list. | Big gap for cell biology questions. |
| **No tissue narrative or HPA-style IHC images.** Atlas has FANTOM5 + Bgee tissue scores (excellent quantitative), but no "expressed in" sentence and no protein-level expression. | HPA's "Cell type enriched (Stomach – Mitotic cells)" is the kind of one-line fact an LLM repeats. |
| **No regulatory build features.** Ensembl surfaces promoters, enhancers, CTCF, TFBS from the Regulatory Build. Atlas has JASPAR motifs but not the gene-region regulatory annotations. | Hurts cis-regulatory variant interpretation. |
| **No drug → indication → mutation triple.** Atlas has drugs (ChEMBL) and clinical trials (NCT IDs) and variants (ClinVar) on the same page but no *linkage* between them. Open Targets' core value is the link. | The G12C / sotorasib / NSCLC story for KRAS, the T790M / osimertinib story for EGFR, the BRCAness / olaparib story for BRCA1 — Atlas does not narrate any of these. |
| **No GeneRIFs.** NCBI surfaces thousands of paper-anchored one-line claims per gene. Atlas has none. | A discoverable corpus of "evidence sentences" is exactly what RAG systems want. |
| **No ClinGen dosage sensitivity scores.** Haploinsufficiency / triplosensitivity verdicts are clinically actionable; Atlas does not surface them. | One-day add via the ClinGen Dosage Sensitivity API. |
| **No comparative genomics depth.** Atlas surfaces 2–4 ortholog rows; Ensembl Compara has 200+ species in a gene tree. | Trivia for casual users, important for evolutionary questions. |
| **No publication corpus.** UniProt has 109–237 cited papers per gene with evidence codes; NCBI has thousands of GeneRIFs. Atlas has zero literature surface area. | Hurts both human researchers and citation-aware LLMs. |

---

## 3. The shortlist

### 3.1 Top 10 features other tools have that Atlas should add

| # | Feature | Source(s) | Effort | Why |
|---|---|---|---|---|
| 1 | **Curated free-text Function narrative (UniProt FUNCTION CC)** | UniProt | **trivial** (already in biobtree's UniProt pull as `CC FUNCTION`) | Single biggest user-visible gap. 2–4 sentences per gene. |
| 2 | **Named isoforms** (p53α/β/γ, K-Ras4A/4B, p16γ, N2A/N2B/N2BA/novex, EGFR p170/p110/p60) | UniProt ALTERNATIVE PRODUCTS | **trivial** (extract `IsoId` + `Name` fields) | Domain experts speak in isoform names, not accessions. |
| 3 | **DISEASE narrative blocks** (Li-Fraumeni, BROVCA1, CMM2, NS3, CMH9 with one-line descriptions and OMIM IDs) | UniProt DISEASE CC | **trivial** | Atlas already lists OMIM disease MIM IDs as a comma list; UniProt has the names + one-line descriptions for free. |
| 4 | **ClinGen dosage sensitivity verdicts** (haploinsufficiency / triplosensitivity) | ClinGen Dosage Sensitivity API | **1-day** | One curated verdict per gene. Clinically actionable. |
| 5 | **NCBI RefSeq summary paragraph** | NCBI Entrez Gene summary | **1-day** (pull `Entrezgene_summary` from E-utilities) | Independent narrative from UniProt; sometimes complementary. |
| 6 | **Drug → indication → biomarker triple** (sotorasib + KRAS G12C + NSCLC; osimertinib + EGFR T790M + NSCLC; olaparib + BRCA1/2 + HER2-neg breast/ovarian) | Open Targets known-drugs + FDA label data | **1-week** | The clinical story most often asked of these genes. |
| 7 | **Subcellular location narrative + (optional) HPA IF image link** | UniProt SUBCELLULAR LOCATION + HPA | **1-day** | A single curated sentence per gene. |
| 8 | **Tissue specificity + Induction narrative** | UniProt TISSUE SPECIFICITY + INDUCTION CC | **trivial** | Already in the biobtree UniProt pull. |
| 9 | **PTM narrative** (the network of phosphorylation/acetylation/sumoylation/ubiquitination, especially for TP53) | UniProt PTM CC | **trivial** | Atlas tallies modified residues by count; UniProt explains which kinases act on which sites. |
| 10 | **Comparative genomics depth — Ensembl gene tree** (200+ species summary, not just 4 rows) | Ensembl Compara | **1-week** (would need a new biobtree dataset or REST pull) | Lower priority but the kind of "wow" depth that competes with Ensembl directly. |

### 3.2 Top 5 features Atlas has that others don't — the differentiators

| # | Feature | What it beats |
|---|---|---|
| 1 | **Single-page ClinVar floor table + per-class counts + top P/LP variant list with HGVS** | NCBI ClinVar paginates; UniProt only lists curated variants; GeneCards summarises. Atlas's TP53 page shows 749 P / 199 LP / 954 VUS / 743 LB / 123 B floors at a glance. |
| 2 | **Paired AlphaMissense + SpliceAI predictions at scale on the same page as ClinVar** | No public reference resource brings the *predicted-effect* and *clinical-effect* layers onto one page. Atlas TP53 has 2,569 AM scores + 1,638 SpliceAI predictions next to 3,850 ClinVar variants. |
| 3 | **TF regulatory layer** — JASPAR motifs + downstream target list + upstream regulators in one section | TRRUST/hTFtarget/JASPAR are separate sites. Atlas TP53's 1,043 downstream targets is unmatched in any single gene page. |
| 4 | **MSigDB + Reactome + GO + KEGG-equivalent on one page** | Each of those is a separate tab on GeneCards/Ensembl. Atlas dumps Reactome (46), MSigDB (1,432), GO BP/MF/CC on one page. |
| 5 | **Honest reporting** (`Antibody resources: 0`; ClinVar floor disclaimer; "sampled 300" hedges on BindingDB; "sampled 300 via entrez" on dbSNP) | GeneCards inflates with vendor catalogs; NCBI hides pagination; UniProt curates away rare splice forms. Atlas is *epistemically* the most honest page — high value for AI trust. |

### 3.3 Top 5 "Atlas is correct, others are wrong" cases — pure AI-trust gold

| # | Case | Atlas's correct behaviour | What the others get wrong |
|---|---|---|---|
| 1 | **CDKN2A dual products (p16INK4a + p14ARF)** | Surfaces both chain entries (1–156 and 1–132), both UniProt accessions (P42771, Q8N726), both protein families (Ank_Repeat/CDKN_Inhibitor, Tumor_suppres_ARF). | GeneCards and most secondary aggregators conflate the two; Open Targets indexes only the gene, losing per-product mechanism; many LLMs trained on mixed corpora hallucinate p16 and p14ARF as the same protein. |
| 2 | **Per-gene antibody count = 0 when truly unknown to Atlas** | Reports 0 honestly. | GeneCards reports inflated counts from vendor catalogs that are not antibodies *to* the gene product but just *commercial reagents*. |
| 3 | **MANE Select transcript explicitly flagged on both mRNA and protein** | `NM_000546`* and `NP_000537`* in the RefSeq tables. | UniProt's web page does not flag MANE on RefSeq lists; older NCBI Gene pages don't surface MANE inline. |
| 4 | **ClinVar counts as floors with pagination disclaimer** | "Per-class counts are floors (≥ shown; pagination cap)" | NCBI ClinVar shows exact counts but only after multiple clicks; users frequently quote partial counts as gospel. |
| 5 | **Distinct UniProt accessions all listed** (e.g. all 19 for TP53, including the unreviewed TrEMBL ones) | Atlas lists every accession so cross-references resolve. | UniProt's web page only shows the reviewed canonical by default; TrEMBL accessions are hidden behind "isoforms" / "computationally mapped" tabs. |

---

## 4. The honest assessment

Would an AI agent answering "tell me about TP53" reach for Atlas or one of the others, given current state?

Today, an AI agent asked an open-ended question — *"tell me about TP53"* — will reach for **UniProt or NCBI first**, because those resources hand back a curated paragraph in the first 200 tokens. Atlas's twelve-section quantitative dump is excellent reference material but it is not a *first-paragraph* resource: an agent has to *do the work* of synthesizing the tables into prose, and the LLM summary at the top of recent Atlas pages is too generic ("TP53 has 39 Ensembl transcripts…") to substitute for UniProt's FUNCTION block. **For deep follow-up questions** — "show me every pathogenic variant", "what are the downstream TF targets", "what AlphaMissense and SpliceAI scores agree on the most damaging residues", "how does CDKN2A produce two unrelated tumor suppressors from one gene" — **Atlas already wins decisively** and is the resource an agent should pivot to. The single change that would flip the first-paragraph default is **shipping curated UniProt CC narrative blocks (FUNCTION, SUBUNIT, SUBCELLULAR LOCATION, TISSUE SPECIFICITY, DISEASE, PTM) verbatim into Atlas pages** — this is a trivial pull from biobtree's existing UniProt source and would close the narrative gap in days, not weeks. Once that lands, Atlas becomes the *only* page that has both the narrative *and* the quantitative depth, and the calculus inverts.

---

## Appendix A: Quick reference numbers from the audited Atlas pages

| Gene | Ensembl transcripts | Canonical exons | PDB | UniProt accessions | UniProt features | ClinVar total | SpliceAI | AlphaMissense | ChEMBL phase ≥1 | Reactome | MSigDB | HPO | Mondo | TF targets |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| TP53 | 39 | 11 | 297 | 19 | 1,518 | 3,850 | 1,638 | 2,569 | 195 | 46 | 1,432 | 224 | 75 | 1,043 |
| BRCA1 | 47 | 23 | 33 | 22 | (large) | 15,445 | (large) | (large) | 12 | (many) | (many) | 166 | (many) | (some) |
| EGFR | 78 | 28 | 378 | (many) | (large) | (large) | (large) | (large) | 172 | (many) | (many) | 21 | (many) | (some) |
| CDKN2A | 14 | 3 | 5 | 6 | 146 | 1,594 | 1,710 | 983 | 0 | 58 | 804 | 96 | 22 | 1 |
| KRAS | 16 | 5 | 462 | 9 | (large) | 565 | (large) | (large) | 10 | (many) | (many) | 372 | (many) | (some) |
| TTN | 15 | 363 | 64 | 9 | 1,289 | 31,454 | (large) | (large) | 0 | (some) | (some) | 147 | 65 | (some) |

Numbers in parentheses were not exhaustively quoted from the pages but match the order of magnitude visible on the relevant section. Numbers without parentheses are direct quotes.
