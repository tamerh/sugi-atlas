# Drug → Indication → Biomarker triple — design spec

**Status:** unimplemented; decision pending
**Owner:** open
**Related:** audit item #6 in [`docs/research/03_page_audit.md`](docs/research/03_page_audit.md),
parked entry in [`docs/research/NEXT.md`](docs/research/NEXT.md) Path C.

## 1. What the triple is

A row that joins three columns about a drug's clinical use:

| Drug | Variant / biomarker | Indication (disease) |
|---|---|---|
| sotorasib (LUMAKRAS) | KRAS G12C | non-small cell lung cancer |
| osimertinib (TAGRISSO) | EGFR T790M | non-small cell lung cancer |
| olaparib (LYNPARZA) | BRCA1/BRCA2 deleterious mutation (BRCAness) | HER2-negative breast cancer; advanced ovarian cancer |

These are precision-medicine narratives — *"drug X is used in disease Y when the patient carries variant Z"* — that are the canonical answer for the clinically-loaded questions an LLM is most often asked of a gene page.

## 2. What Atlas has today

Atlas's gene §10 and disease §13 already carry every column **independently:**

- ChEMBL phased molecules per gene (drug names, max phase)
- Disease-level clinical trials per Mondo id (NCT IDs, trial phases)
- ClinVar P/LP variants per gene (HGVS, classifications)
- BRENDA EC + UniProt features for the variant residues
- CIViC gene-level + intOGen somatic-driver flags for cancer cohorts

**What's missing is the join.** No row currently says "sotorasib + KRAS G12C + NSCLC" as a single fact.

## 3. Why the join is hard

Three independent sources, four normalization problems:

| Source | Carries | Missing for our purposes |
|---|---|---|
| **ChEMBL** | drug → indication (Mondo / MeSH / EFO) via clinical_trials | no variant column |
| **FDA approved-drug labels** | drug → variant requirement in free-text English at the top of the label | not in biobtree; would require label parsing |
| **ClinVar** | variant → gene + pathogenicity classification | no drug or indication context |

### Normalization issues

**(a) Drug name fuzziness.** Same drug appears as `sotorasib` / `AMG-510` / `LUMAKRAS` / `CHEMBL4439653` / `CHEMBL5314110` (salt forms / brand / generic / development code). Same problem we already addressed for `chembl_molecule.parent` (BIOBTREE_ISSUES #15 partial-fix), multiplied across sources.

**(b) FDA labels are free text.** Approved-use biomarker statements read like *"approved for adult patients with KRAS G12C-mutated locally advanced or metastatic NSCLC"*. Extracting the structured `{drug, variant, indication}` triple requires parsing English with rules or an LLM, and the format varies per drug.

**(c) Variant nomenclature is fragmented.** The same KRAS G12C is spelled:
- `p.G12C` (HGVS protein)
- `c.34G>T` (HGVS cDNA)
- `chr12:25245350G>T` (genomic)
- `rs121913530` (dbSNP)
- `VCV000045132` (ClinVar)

Each source uses a different one. Cross-linking requires variant normalization.

**(d) Disease nomenclature is fragmented.** NSCLC appears as:
- "NSCLC" / "non-small cell lung cancer" / "non-small cell lung carcinoma" (free text)
- `MONDO:0005233`
- `EFO:0003060`
- MeSH `D002289`
- ICD-10 `C34.x`

ChEMBL indication labels use a mix of MeSH + EFO; biobtree's mondo>>mesh+mondo>>efo edges resolve most cross-walks but not all.

## 4. Candidate solutions

### Option A — Consume Open Targets

Open Targets does this exact aggregation as their core product. Their `targets/{ensembl_id}/knownDrugs` GraphQL endpoint returns rows like:

```json
{
  "drug": {"name": "sotorasib", "chemblId": "CHEMBL4439653"},
  "disease": {"name": "non-small cell lung carcinoma", "id": "EFO_0003060"},
  "phase": 4,
  "status": "Completed",
  "mechanismOfAction": "KRAS GTPase inhibitor",
  "ctIds": ["NCT04303780", "NCT03600883"]
}
```

…but with **biomarker context attached via the precision-medicine pipeline**: their `target → mechanismOfAction` + `target → variantEffect` joins surface the G12C requirement.

**Pros**
- Already does the four normalizations (drug / variant / indication / source-of-evidence).
- Refreshed quarterly with FDA + ChEMBL + ClinVar updates.
- One GraphQL call per gene; cheap.
- Open Targets ID space (`EFO:...`, `CHEMBL...`, `ENSG...`) already cross-walks with what Atlas carries.

**Cons**
- Network dependency on a third-party API in the collect step (same reliability concern that ruled out direct UniProt REST earlier — see [`docs/research/NEXT.md`](docs/research/NEXT.md) "Important — needs design decision").
- Atlas would inherit Open Targets' curation decisions (good and bad) rather than running its own.
- Single-source dependency: if Open Targets' API shape changes or their service goes down, Atlas's §10 loses content.
- Their licensing is CC-BY 4.0 — attribution required on every page; fine but worth noting.

**Mitigation paths for the network concern**
- File as a biobtree feature request — biobtree ingests Open Targets nightly, exposes as a new dataset (`open_targets_drug` / `open_targets_indication` edges). Same pattern that resolved UniProt CC (BIOBTREE_ISSUES #9). **This is the cleanest fit with Atlas's "biobtree is single source" architecture.**
- Atlas-side on-disk cache (`<dist>/cache/open_targets/<symbol>.json`) with weekly refresh.

### Option B — Direct DIY ingestion

Atlas builds the three-way join from raw sources:
1. ChEMBL drug ↔ indication via existing biobtree edges.
2. FDA approved-drug labels parsed locally from `https://dailymed.nlm.nih.gov/` or `https://api.fda.gov/drug/label.json`. One label per approved drug; the biomarker requirement lives in the "Indications and Usage" section.
3. ClinVar variant ↔ gene via existing biobtree edges; normalize HGVS p.G12C ↔ rs121913530 via biobtree's dbsnp dataset.
4. Drug-name normalization across ChEMBL + FDA via ChEMBL `parent` + `altNames` (BIOBTREE_ISSUES #15) + RxNorm if needed.

**Pros**
- Atlas owns the curation; no third-party API dependency.
- Each layer is independently versionable.

**Cons**
- FDA label parsing is the hard part — free-text English with no structured biomarker field. Would either need rule-based parsers (brittle, drug-by-drug) or an LLM extraction step (slow, expensive, hallucination surface).
- Variant normalization between FDA's prose ("KRAS G12C-mutated") and ClinVar's HGVS is non-trivial; biobtree's dbsnp edge helps but doesn't cover protein-level shorthand consistently.
- Coverage: only ~30-50 drugs have FDA biomarker-restricted approvals today (compared with thousands of phased ChEMBL molecules), so the ingestion infrastructure is heavy relative to the row count produced.
- Maintenance: FDA labels are updated as new approvals roll in; Atlas would need a refresh job.

**Estimated scope:** weeks of focused work, not days. The audit's "1-week" estimate (item #6) assumes consuming a precomputed source.

### Option C — Hybrid: Open Targets via new biobtree dataset

The recommended path:

1. **File BIOBTREE_ISSUES #18:** request biobtree ingest Open Targets' `knownDrugs` table as a new dataset, with edges:
   - `>>chembl_target>>open_targets_drug` (drug per target)
   - `>>open_targets_drug>>mondo` (indication)
   - `>>open_targets_drug>>clinvar` (biomarker variant, where Open Targets has linked one)
2. **Atlas waits for the ingest** (same path that worked for UniProt CC #9, intOGen, CIViC, etc.). When it lands, Atlas's existing §10 drug-targets collector picks it up via cohort fan-out with a small renderer addition.
3. **Until then:** Atlas's §10 surfaces a one-line stub *"For drug-indication-biomarker linkage see Open Targets: `https://platform.opentargets.org/target/{ensembl_id}`"* — explicit attribution, no false claim that Atlas does the join.

**Pros**
- Matches existing Atlas architecture: biobtree is single source, no network in collect.
- Curation responsibility stays with Open Targets (who do it full-time).
- File-and-forget once the biobtree request lands.

**Cons**
- Wait time for upstream biobtree work; depends on biobtree dev priorities.
- Atlas stays *missing* this view until then (current state).

### Option D — Stay missing; link out

Surface the explicit link-out stub in §10 without trying to build the join at all. Lowest-effort; matches "honest reporting" principle (audit's differentiator #5).

**Pros**
- Zero engineering cost.
- No false-precision risk.

**Cons**
- Audit item #6 stays open; the question *"what drugs treat KRAS-G12C NSCLC"* sends the AI agent to Open Targets, not Atlas.

## 5. Recommendation (for when the decision is made)

**Lean:** Option C (file as biobtree feature request → wait → consume via existing collector). Matches the pattern that already worked for the audit's #1 gap (UniProt CC #9). Until then, Option D's link-out stub keeps Atlas honest about the gap.

**If Option C blocks for too long** and the audit-#6 gap becomes a launch blocker, fall back to Option A with an on-disk cache.

**Option B (direct ingest) is not recommended** unless biobtree explicitly cannot ingest Open Targets — the FDA-label parsing cost is high and the row count is small.

## 6. Open questions before any work starts

- Confirm Open Targets' licensing (CC-BY 4.0) is compatible with Atlas's intended distribution.
- Confirm biobtree can ingest a CC-BY 4.0 source; check existing datasets for licensing precedent.
- Survey whether `Open Targets > knownDrugs` actually carries the variant-biomarker column for the canonical examples (sotorasib/G12C, osimertinib/T790M, olaparib/BRCAness) or whether their variant linkage is shallower than the audit assumed. **One probe call** confirms or refutes.
- Define renderer placement: new §10 subblock, new §11 subblock, or insert into the cancer-overview lead.

## 7. Out of scope for this spec

- Generalizing the same triple-join pattern to non-cancer precision medicine (CYP enzyme variants + drug dosing — already partially covered by PharmGKB which biobtree already ingests).
- Variant-level CIViC ingestion — separate item (BIOBTREE_ISSUES.md "civic_variant / civic_evidence / civic_assertion" — routes return n=0 today; tracked separately).
