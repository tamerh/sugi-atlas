# 02 — Mining biobtree for the Sugi Atlas gene-page collector

**Date:** 2026-05-30
**Scope:** audit the biobtree data graph behind `http://127.0.0.1:8000/api`, identify
datasets the Atlas collector (12 sections under `src/atlas/gene/sections/`) is not
yet wiring, probe each candidate empirically against TP53 (HGNC:11998 / P04637 /
ENSG00000141510) and BRCA1 (HGNC:1100 / P38398), and rank what to wire next.

The audience is the next-session implementer. Every recommendation below was
verified with a live `curl` probe; the chain, row counts, and a sample row are
quoted in-line so they can be replayed without re-discovery.

---

## 1. Coverage snapshot

The bundled `dataset_coverage.py` still points at a legacy `collect.py`; the
fixed scan (union all `>>edge` and quoted-name tokens across
`src/atlas/gene/sections/*.py`) gives:

- **117** edge-graph nodes (per `/api/help?topic=edges`)
- **51** used by Atlas today
- **66** uncovered

### Verdict table — every uncovered node

`a` = wire it, real data verified · `b` = xref-only, no payload — don't wire ·
`c` = needs upstream biobtree fix (see `BIOBTREE_ISSUES.md`) · `d` = wire
conditionally / projection or scope gotcha · `e` = hierarchy node (parent/child)
— wire only inside the parent dataset's existing routes.

| Dataset | Verdict | Notes |
|---|---|---|
| `mirdb` | **a** | refseq→mirdb yields 77–100 miRNAs targeting the gene; rich attrs (score/breadth). High value. |
| `pubchem_activity` | **a** | 100 bioassay outcomes for P04637 with Ki/IC50, qualifier, unit, PMID. |
| `chembl_activity` | **a** | per-molecule activity values (IC50 etc.) — accessible from `chembl_molecule`, not from a uniprot start. |
| `chembl_assay` | **d** | only reachable via `chembl_target>>chembl_assay`, which currently returns not_found from a TP53 uniprot start (long-chain failure mode); needs upstream check. |
| `chembl_document` | **d** | reachable only from `chembl_assay`; blocked by the same long-chain issue. |
| `chebi` | **a** (drug-side) | `chembl_molecule>>chebi` resolves; gives canonical chemical id + formula. |
| `pubchem` | **a** (drug-side) | `chembl_molecule>>pubchem` resolves; CID, name, formula, SMILES-bearing entry. |
| `patent_compound` | **b** | `entry` returns `Empty` — xref-only (just a count of patents per CID). |
| `clinical_trials` via molecule | **a** | already used disease→trials; **molecule→trials** is new and gives drug-development-stage trials per drug. |
| `pharmgkb` | **a** | `>>hgnc>>pharmgkb_gene>>pharmgkb` returns the per-chemical/drug-class rows (counts of variant/clinical/guideline annotations). |
| `pharmgkb_variant` | **a** | reachable via `>>hgnc>>entrez>>dbsnp>>pharmgkb_variant` — 2 TP53 rows (rs1042522, rs4968187) with drug lists. |
| `pharmgkb_clinical` | **c** | every chain tried returns not_found from TP53/BRCA1; suspect same hgnc-routing issue as dbsnp had. |
| `pharmgkb_guideline` | **c** | as above. |
| `pharmgkb_pathway` | **a** (rare) | one BRCA-route gave 1 hit (Venetoclax pathway) — gene-conditional, wire lightly. |
| `ctd_gene_interaction` | **a** | `>>hgnc>>entrez>>ctd_gene_interaction` returns 100 chemical-gene interactions with action verbs and PMID count. |
| `ctd` | **a** (chemical drilldown) | richer chemical metadata (CAS#, synonyms) accessible from `ctd_gene_interaction>>ctd` rows. |
| `ctd_disease_association` | **d** | `>>hgnc>>mim>>ctd_disease_association` empty — likely needs a different routing or upstream fix. |
| `mesh` | **a** | `>>hgnc>>clinvar>>mondo>>mesh` returns 35 rows for TP53 (Li-Fraumeni etc.) with descriptor + tree numbers + scope_note. |
| `efo` | **a** | `>>hgnc>>gwas>>efo` returns 25 EFO terms for TP53 GWAS traits (cancer types). |
| `uberon` | **a** | `>>ensembl>>bgee>>uberon` returns 100 anatomical structures where expression is detected; canonical names + synonyms. |
| `cl` | **a** | `>>ensembl>>bgee>>cl` returns 10 cell types; clean canonical names. |
| `hmdb` | **b/c** | uniprot→hmdb empty for both probes; HMDB is metabolite-keyed — not gene-relevant for proteins like TP53. |
| `rnacentral` | **c** | `>>uniprot>>rnacentral` empty for TP53 — RNAcentral is for ncRNAs, so empty is *expected* for protein-coding genes but Atlas should still try once. |
| `ena` | **b** | reachable only via rnacentral; same scope. |
| `brenda` | **a** (enzymes only) | EGFR (P00533) returns EC 2.7.10.1 with synonyms, kinetics/inhibitor xref counts; TP53 empty (correct — not an enzyme). |
| `brenda_kinetics` | **a** | reachable from brenda; KM/Vmax values when present. |
| `brenda_inhibitor` | **a** | reachable from brenda; inhibitor compounds with values. |
| `gtopdb` | **a** (targets only) | EGFR returns family + type; TP53 empty (correct — not a Guide-to-Pharm target). Wire as conditional. |
| `gtopdb_interaction` | **d** | reachable from gtopdb but the projection currently returns near-empty rows (`1797_10044||||||`); upstream projection bug. |
| `gtopdb_ligand` | **c** | `>>gtopdb>>gtopdb_ligand` returns not_found; edges claim the link but it doesn't resolve. |
| `rhea` | **a** (enzymes only) | EGFR returns 1 RHEA reaction; TP53 empty (correct). |
| `swisslipids` | **a** (lipid genes only) | empty for TP53 (correct); already used downstream by lipid-handling genes (FASN/ELOVL families). |
| `lipidmaps` | **b** | reachable only via chebi/pubchem; metabolite-only. |
| `fantom5_enhancer` | **d** | resolves (31 rows for TP53) with real attrs (chromosome, start/end, TPM, samples), **but** `associated_genes` is a window-wide gene list (74 genes per enhancer for the TP53 region), so "TP53 enhancers" is misleading. Wire only if filtered by `associated_genes contains symbol`. |
| `cellphonedb` | **c** | `>>hgnc>>cellphonedb` and `>>uniprot>>cellphonedb` both empty for TP53 — same shape as the now-fixed dbsnp issue. |
| `cellxgene` | **d** | edges hub through `cl/uberon/mondo`; no direct gene route — accessible via disease(`mondo>>cellxgene`) or anatomy hubs only. Lower priority. |
| `cellxgene_celltype` | **d** | same. |
| `scxa_expression` | **c** | resolves (100 rows) but rows are for OTHER genes within the experiment (`ENSG00000000003_E-GEOD-75140`) — projection bug similar to `scxa_gene_experiment`. |
| `scxa_gene_experiment` | **c** | already known projection bug — wire is blocked until biobtree fix. |
| `bao` | **e** | assay ontology — wire only alongside `chembl_assay`. |
| `eco` | **e** | evidence-code ontology — wire only inside annotations that carry ECO refs. |
| `taxonomy` | **e** | only useful for ortholog labels (already implicit via `ortholog`). |
| `*parent/*child` (go/reactome/mondo/hpo/uberon/cl/efo/mesh/interpro/tax/bao/eco) | **e** | hierarchy navigation; wire inline (e.g. `>>go>>goparent` for term lineage), not as standalone. |
| `pubmed` | **b** | `entry` is `Empty` — pubmed is xref-only (PMID without title/abstract). Use PMIDs as evidence pointers only. |

Counting verdicts: **a = 17 standalone-wire** datasets, **d/e = 14 conditional or
hierarchy** datasets, **b/c = the remaining 18** are either xref-only or blocked.

---

## 2. Probe receipts

All probes were run against TP53 and (where relevant) BRCA1 / EGFR (P00533).
Each row shows the chain, the `rows:` count returned by `/api/map`, a sample
row, and an `entry` test confirming attribute payload.

### 2.1 The clear wins (verdict **a**)

**mirdb — miRNAs targeting the gene**
```
P04637  >>uniprot>>refseq>>mirdb        rows: 77
HGNC:11998  >>hgnc>>refseq>>mirdb       rows: 100  (paginates)
  sample: HSA-MIR-300|hsa|2856|71.76|99.92   (id|species|n_targets|avg_score|max_score)
entry HSA-MIR-300: Mirdb {mirna_id, species, target_count, avg_score, max_score, top_targets[]}
```
Real payload. Per-miRNA top-target lists are present in the entry — useful for
"miRNAs that target TP53 and what else they target."

**pubchem_activity — protein-level bioassay actives**
```
P04637  >>uniprot>>pubchem_activity     rows: 100
  sample: 10008863_614426_1|Active|ki|1.24|uM
entry 10008863_614426_1: PubchemActivity {cid, aid, activity_outcome,
  activity_type, qualifier, value, unit, protein_accession, gene_id, pmid}
```
PMID is embedded — instant evidence trail. Mostly low/early-stage chemical
matter against TP53, but for kinase targets (EGFR, ALK) this becomes
high-volume real bioassay data.

**chembl_molecule drilldowns** — already-fetched drug list now becomes richer.
For each phased drug:
```
CHEMBL941 (imatinib) >>chembl_molecule>>chembl_activity        rows: 100
  sample: CHEMBL_ACT_24962908|IC50|0.06|nM|10.22
CHEMBL941 >>chembl_molecule>>clinical_trials                   rows: 100
  sample: NCT00081926|Gleevec Trial ...|COMPLETED|PHASE4|INTERVENTIONAL
CHEMBL941 >>chembl_molecule>>chebi                             rows: 1
  CHEBI:45783|imatinib|C29H31N7O|3
CHEMBL941 >>chembl_molecule>>pubchem                           rows: 1
  5291|Imatinib|C29H31N7O|true|drug|Imatinib
CHEMBL941 >>chembl_molecule>>pubchem>>patent_compound          rows: 3
```
All real. Lets §10 go from "drug list with phase" to "drug list with phase +
canonical chebi/pubchem ids + trial titles + patent count."

**ctd_gene_interaction — chemical-gene literature interactions**
```
HGNC:11998 >>hgnc>>entrez>>ctd_gene_interaction                rows: 100
  sample: C000228_7157_9606|aristolochic acid I|TP53|Homo sapiens
          | increases^expression;increases^phosphorylation|2
entry: CtdGeneInteraction {chemical_id, chemical_name, gene_symbol,
  interaction (verbose), interaction_actions[], gene_forms, pubmed_count}
```
The action verbs are CTD's controlled vocabulary (`increases^expression`,
`decreases^methylation`, etc.) — perfect for grouping. PubMed count is per row.

**mesh — disease descriptors**
```
HGNC:11998 >>hgnc>>clinvar>>mondo>>mesh                        rows: 35
HGNC:1100  same chain                                          rows: 15
  sample: D016864|Li-Fraumeni Syndrome|1|C04.700.600;...|false||
entry D016864: Mesh {descriptor_ui, descriptor_name, entry_terms[],
  tree_numbers[], scope_note}
```
MeSH tree numbers give the broader-disease hierarchy "for free" — collapsing
to top-level MeSH categories (C04 = neoplasms, C16 = congenital, …) is a
natural rendering choice.

**efo — GWAS trait ontology**
```
HGNC:11998 >>hgnc>>gwas>>efo                                   rows: 25
  sample: EFO:0005922||esophageal squamous cell carcinoma
entry EFO:0005922: Ontology {name, synonyms[], id}
```
Already redundant with `gwas` titles, but synonyms enable cleaner dedup.

**uberon / cl — anatomy & cell-type names for expression**
```
ENSG00000141510 >>ensembl>>bgee>>uberon                        rows: 100
  sample: UBERON:0003053||ventricular zone
ENSG00000141510 >>ensembl>>bgee>>cl                            rows: 10
  sample: CL:0000576||monocyte
entry UBERON:0003053: Ontology {name, synonyms[]}   (real)
entry CL:0000576:    Ontology {name}                (thin but real)
```
Replaces Atlas's current habit of dropping bare `bgee_evidence.anatomy_id`
codes into the bundle — now we get canonical names + synonyms.

**pharmgkb hub (per-chemical drug summary)**
```
HGNC:11998 >>hgnc>>pharmgkb_gene>>pharmgkb                     rows: 7
  sample: PA166122986|radiotherapy|Drug Class|23|105|0
           (id|name|type|n_variants|n_clinical|n_guidelines)
```
The aggregate counts per drug are *built in* — Atlas already pulls
`pharmgkb_gene` but stops at the per-gene summary. Going one hop further
gives the per-drug breakdown.

**pharmgkb_variant via dbsnp** *(undocumented but works)*
```
HGNC:11998 >>hgnc>>entrez>>dbsnp>>pharmgkb_variant             rows: 2
  sample: PA166155173|rs1042522|TP53,WRAP53|3|5.00|2
           |antineoplastic agents;Platinum compounds
```
This is the *only* working chain into `pharmgkb_variant` from a gene anchor.
Mirrors the dbsnp workaround already in §6.

### 2.2 Multi-hop chains worth wiring (compound routes)

| Route | n rows TP53 | Adds |
|---|---|---|
| `>>chembl_molecule>>clinical_trials` (per phased drug) | 100/drug | Trial-level adverse-event-free names + phase |
| `>>chembl_molecule>>pubchem` + `>>chembl_molecule>>chebi` | 1/drug | Canonical chem ids — needed for cross-DB linking |
| `>>uniprot>>pdb>>pubmed` | 100 | Structure-paper PMIDs (literature for §4 structure) |
| `>>uniprot>>pubchem_activity` | 100 | Bioassay actives w/ Ki/IC50/PMID |
| `>>uniprot>>jaspar>>pubmed` | 1 | TF-binding paper (for §9 TF regulation) |
| `>>hgnc>>entrez>>ctd_gene_interaction` | 100 | Chemical-gene literature interactions |
| `>>hgnc>>clinvar>>mondo>>mesh` | 35 | MeSH disease descriptors w/ tree |
| `>>hgnc>>gwas>>efo` | 25 | EFO-normalized GWAS traits |
| `>>ensembl>>bgee>>uberon` | 100 | Anatomy names for bgee tissues |
| `>>ensembl>>bgee>>cl` | 10 | Cell-type names for bgee |
| `>>hgnc>>entrez>>dbsnp>>pharmgkb_variant` | 2 | Pharmacogenomic variants w/ drug list |

Failed compounds worth noting:
- `>>uniprot>>chembl_target>>chembl_molecule` (from HGNC start) → `not_found`
  even though `>>uniprot>>chembl_target` succeeds. Long-chain `not_found`
  matches BIOBTREE_ISSUES #1/#4 — workaround: collect chembl_molecules from a
  uniprot start (already done in §10) and *fan out one hop per molecule*
  rather than chaining.
- `>>hgnc>>mim>>ctd_disease_association` → empty.
- `>>chembl_molecule>>chembl_activity>>chembl_assay` → not_found.
  Same long-chain pattern. Need two-step: fetch activities, then per-activity
  resolve assay.
- `>>gtopdb>>gtopdb_ligand` → not_found (edge advertised, doesn't resolve).

---

## 3. Section fit — proposals

### §6 variants — add `pharmgkb_variant`
- **route:** `>>hgnc>>entrez>>dbsnp>>pharmgkb_variant`
- **bundle key:** `pharmgkb_variants: [{pa_id, rsid, gene_symbols, evidence_level, drug_classes}]`
- **render:** "Pharmacogenomic variants: rs1042522 (P04637 codon 72 polymorphism, level 3, modifies response to antineoplastic / platinum compounds)."
- **gotcha:** the gene_symbols field may include neighbors (TP53,WRAP53) — filter by exact symbol if strict.

### §7 pathways — add ontology lineage
- **routes:** `>>go>>goparent` (per top-N GO terms), `>>reactome>>reactomeparent` (per pathway). Lets Atlas roll up top-level GO/Reactome categories.
- **bundle key:** `go_top_categories: [{id, name}]` — counts of GO terms grouped by their L1/L2 ancestors (e.g. "apoptotic process" rolls up to "programmed cell death").
- **gotcha:** parent calls are one extra hop per term — batch to ~top 20 to keep latency sane.

### §8 interactions — keep; no obvious new datasets (already saturated)

### §9 TF regulation — add `jaspar>>pubmed`
- **route:** `>>uniprot>>jaspar>>pubmed`
- **bundle key:** `jaspar_pmids: [...]` — evidence trail for binding-site matrices already in the bundle.

### §10 drugs — biggest section to enrich
1. **Per phased drug:** drill `chembl_molecule>>chembl_activity` (top 5 by pchembl), `>>clinical_trials` (trial titles + phase), `>>pubchem` and `>>chebi` (canonical ids).
   - **bundle keys:** add per-drug `top_activities`, `trials_sample`, `pubchem_cid`, `chebi_id`.
   - **render:** "Imatinib (CHEMBL941, ChEBI:45783, CID 5291) — phase 4 trials: NCT00081926, NCT00372476…"
2. **Add `pharmgkb` per-chemical breakdown** (`>>hgnc>>pharmgkb_gene>>pharmgkb`) with the in-built variant/clinical/guideline counts.
3. **Bioassay actives** (`>>uniprot>>pubchem_activity`): top-N by potency, with PMIDs.
4. **Conditional brenda/gtopdb** — wire but only render if hits exist (TP53 → both empty; EGFR → both populated).
5. **CTD chemical-gene interactions** as a "chemicals known to modulate the gene" subsection (`>>hgnc>>entrez>>ctd_gene_interaction`) — distinct from `chembl` (in vitro binding) and from `pharmgkb` (pharmacogenomics).
   - **render:** "Increases TP53 expression / phosphorylation: aristolochic acid I (2 PMIDs); …"

### §11 expression — add anatomy/celltype names
- **routes:** `>>ensembl>>bgee>>uberon`, `>>ensembl>>bgee>>cl`
- **bundle keys:** rewrite current `tissues` items to include `uberon_name` / `cl_name`; add a separate `cell_types_summary`.
- **gotcha:** these are *all* uberon/cl terms ever annotated for that gene by bgee — not ranked. To rank, keep using `bgee_evidence` scores and join on the uberon/cl id.

### §11 expression — add `mirdb` (or move to §9 / new sub-section)
- **route:** `>>hgnc>>refseq>>mirdb` (or `>>uniprot>>refseq>>mirdb`)
- **bundle key:** `mirna_regulators: [{mirna_id, target_count, max_score, avg_score}]`
- **render:** "Top miRNAs predicted to target TP53: hsa-miR-300 (max score 99.9), hsa-miR-3922-5p (98.8)…"
- **gotcha:** paginates at 100; for a comprehensive list, walk the cursor.

### §12 diseases — add `mesh`, `efo`
- **routes:** `>>hgnc>>clinvar>>mondo>>mesh`, `>>hgnc>>gwas>>efo`
- **bundle keys:** `mesh_descriptors` (with `tree_numbers` for rollup), `efo_traits`
- **render:** "MeSH disease categories: Neoplastic Syndromes Hereditary (C04.700); Li-Fraumeni Syndrome (C04.700.600); …"
- **gotcha:** mesh entry attrs include `scope_note` — usable for one-line authoritative blurbs.

### §10 drugs — also add chembl_molecule>>pubchem>>patent_compound *count*
- Patent_compound's entry is `Empty` (xref-only) but the count of patents per molecule is a useful "drug commercial activity" signal. Just project the row count, don't `entry`.

---

## 4. Section-overlap warnings

- **CTD vs ChEMBL vs PharmGKB.** All three relate chemicals to TP53; semantics
  differ. ChEMBL = curated bioactivity (mostly in vitro binding). PharmGKB =
  pharmacogenomics (variants modulate drug response). CTD = literature-mined
  causal modulation ("compound X increases/decreases gene Y expression").
  Put them in §10 as separate sub-blocks; do not merge.
- **MeSH vs MONDO.** Already partly in Atlas via mondo titles. MeSH adds tree
  hierarchy (C04 / C16) and scope_notes. Use MeSH to roll up, MONDO to drive
  detail.
- **bgee uberon/cl vs scxa.** scxa is currently blocked by a projection bug
  (returns rows for OTHER genes in the experiment). bgee uberon/cl is the
  unblocked path for "where is this gene expressed."
- **fantom5_enhancer.** Resolves but `associated_genes` is window-wide (~70
  genes per enhancer for the TP53 locus). Useless unless filtered to
  enhancers where TP53 is the **nearest** gene, which biobtree doesn't surface
  directly. Skip unless that filter is added upstream.

---

## 5. Ranked backlog — top 10 by value × effort

Effort assumes: chain already returns real rows, the bundle key fits an
existing section, and the renderer is a 1–2 line projection. "L" = ≤1 h
implement, "M" = ≤half-day (extra fan-out / per-row drilldown), "H" =
multi-hop with caching needed.

| # | Dataset / route | Section | Adds | Effort | Why it matters for an AI-agent target |
|---|---|---|---|---|---|
| 1 | `mirdb` via `>>hgnc>>refseq>>mirdb` | §11 (or new §9b) | 77–100 miRNAs/gene w/ scores | **L** | First-class regulatory layer Atlas is missing entirely; agents currently hallucinate miRNA regulators. Pure data, gene-keyed, no scope mixing. |
| 2 | `chembl_molecule` drilldowns (`>>chembl_activity` top-N + `>>clinical_trials` + `>>pubchem` + `>>chebi`) per phased drug | §10 | Per-drug potency, trials, canonical chem ids | **M** | Turns the drug list from "names + phase" into a structured pharma block agents can cite. Each drug becomes a mini-record. |
| 3 | `>>uniprot>>pubchem_activity` | §10 | 100 bioassay actives w/ Ki/IC50/PMID | **L** | High-volume for kinases/enzymes; the PMID field is gold for evidence chains. |
| 4 | `>>hgnc>>entrez>>ctd_gene_interaction` | §10 | 100 chemical-gene literature interactions w/ action verbs | **L** | Distinct semantics from ChEMBL — "what chemicals modulate this gene per literature." Action verbs are CV-controlled and groupable. |
| 5 | `>>hgnc>>clinvar>>mondo>>mesh` (+ tree rollup) | §12 | MeSH disease descriptors with tree codes | **L** | The only authoritative way to roll mondo diseases up to top-level disease categories. Tree numbers come free. |
| 6 | `>>ensembl>>bgee>>uberon` and `>>ensembl>>bgee>>cl` | §11 | Canonical anatomy / cell-type names for current bgee blocks | **L** | Replaces bare UBERON / CL ids with names + synonyms; instant readability win. |
| 7 | `>>hgnc>>pharmgkb_gene>>pharmgkb` | §10 | Per-chemical PharmGKB summary (n_variants / n_clinical / n_guidelines) | **L** | The per-drug breakdown PharmGKB pre-computes; one hop deeper than the current Atlas stop. |
| 8 | `>>hgnc>>entrez>>dbsnp>>pharmgkb_variant` | §6 | Pharmacogenomic variant rows w/ rsid + drug classes | **L** | Bridges variants↔drugs — the one route into pharmgkb_variant that actually returns rows. |
| 9 | `>>hgnc>>gwas>>efo` | §12 | EFO-canonical GWAS trait names + synonyms | **L** | Cheap dedup + clearer naming for the gwas block; agents prefer normalized vocab over GWAS-catalog free text. |
| 10 | Conditional `>>uniprot>>brenda` (+ `brenda_kinetics`, `brenda_inhibitor`, `>>uniprot>>rhea`) and `>>uniprot>>gtopdb` (+ `gtopdb_interaction` once projection fixed) | §3 / §10 | Enzyme/RHEA reactions for enzymes; GtoPdb family/type for receptors/channels | **M** | Class-specific richness: lights up kinase / GPCR / channel pages without polluting non-enzyme pages (TP53 stays clean). |

### Honorable mentions, not in top 10

- **GO / Reactome parent rollups** (`>>go>>goparent`, `>>reactome>>reactomeparent`) — easy to add, but the data is already implicit in the term names. Lower marginal value vs (1)–(10).
- **`>>uniprot>>jaspar>>pubmed`** — single PMID per gene; nice evidence link in §9 but small.
- **`>>chembl_molecule>>pubchem>>patent_compound` row count** — one-line "commercial-activity" badge per drug. Trivial to add as an integer.

### Blocked / file upstream

- `pharmgkb_clinical`, `pharmgkb_guideline` — all chains tried returned
  `not_found`. Add to BIOBTREE_ISSUES (or extend existing #1) — "hgnc and
  pharmgkb_gene starts don't reach pharmgkb_clinical / pharmgkb_guideline; need
  a documented working chain."
- `scxa_expression`, `scxa_gene_experiment` — projection returns rows for
  *other* genes in the experiment, not the queried gene. Known.
- `gtopdb_interaction` — resolves but projection rows are mostly empty
  pipe-delimited fields. Needs upstream fix before wiring.
- `gtopdb_ligand` — edges advertise the link but `>>gtopdb>>gtopdb_ligand`
  returns `not_found`.
- `cellphonedb` — `>>hgnc>>cellphonedb` and `>>uniprot>>cellphonedb` both
  empty; matches the previously-fixed dbsnp pattern. Worth a follow-up
  ticket.

---

## 6. Quick-start for the next implementer

1. **Don't change `dataset_coverage.py` until you fix the `collect.py`
   reference** — the file expects a monolithic collector; the sections were
   split. The 10-line scan in §1 (`union of >>X and "X" tokens across
   sections/*.py`) is the working substitute and should be folded back into
   the tool.
2. **Start with #1, #5, #6 from the backlog.** Each is a single chain, no
   per-row fan-out, fits cleanly into an existing section. Combined they push
   coverage from 51/117 to ~57/117 in one afternoon.
3. **#2 needs a per-drug fan-out pattern** — copy the per-uniprot loop in
   §7's reactome union (it iterates `a.reviewed_uniprots`); the equivalent
   for §10 is iterating the phased-drug list and resolving each.
4. **Avoid long compound chains.** Every probe that combined ≥4 hops from a
   non-native start (e.g. `>>hgnc>>...>>chembl_activity`) returned
   `not_found`. The robust pattern is: resolve one anchor, then fan out one
   hop at a time, joining client-side. The Atlas codebase already follows
   this style — keep it.
5. **`entry` xref counts are free.** Atlas already uses `xref_counts(a.hgnc_entry)`
   for totals. Do the same on uniprot / ensembl entries to size sections
   before fetching (e.g. skip pubchem_activity entirely if the protein has 0
   pubchem_activity xrefs).

**Expected post-implementation coverage:** ~62/117 nodes wired (+ ~6 of the
hierarchy/conditional ones). The remaining gap is dominated by xref-only ids
and upstream-blocked datasets — i.e. the realistic ceiling.
