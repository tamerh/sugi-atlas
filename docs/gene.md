# How a gene page is built

This walks through exactly how Atlas turns a gene symbol into a published page —
the biobtree chains it queries, the sources behind each section, and the
curation rules that decide what is shown and what is suppressed. For the shared
machinery (the biobtree transport, the collect→render split, the page contract,
JSON-LD), see [how-it-works.md](how-it-works.md); this doc is gene-specific.

The build is two passes: **collect** resolves the gene once and runs twelve
section collectors into a structured bundle (no model, no prose); **render**
arranges that bundle into eight canonical, frozen page zones.

```
HGNC symbol
   │  resolve()            one set of anchor IDs, shared by all sections
   ▼
Anchors ─►  12 section collectors  (biobtree chains → structured JSON bundle)
   │
   ▼
render_all ─►  8 canonical H2 zones  (sections demoted to H3, stable {#anchors})
   │
   ▼
assemble_page ─►  frontmatter · ## Summary · body · ## Related Atlas pages
```

---

## 1. Resolving the gene

Everything starts from one resolution step (`gene/anchors.py`). A gene symbol
(or HGNC id) is resolved to an immutable `Anchors` record that every section
then reads — so the ~24 redundant HGNC/Ensembl/UniProt lookups the old
per-section code paid collapse to one.

**Symbol → HGNC** (`resolve_hgnc`): a `search(symbol, source="hgnc")` keeps
candidates matching `HGNC:\d+`. Two cases need care:

- **Ambiguous symbols** — `AR` is both *amphiregulin* and *androgen receptor*.
  When more than one HGNC row matches, each candidate's `entry` is fetched and
  its approved-symbol list is checked against the query; only the true match
  wins.
- **Failure is typed** — an unresolvable symbol raises `ValueError`, not
  `sys.exit`, so the batch driver skips and logs that one gene rather than
  crashing the whole corpus.

From the HGNC id, the anchor set is assembled in one place:

| Anchor | Source | Chain |
|---|---|---|
| `ensembl_id` | Ensembl | `>>hgnc>>ensembl` |
| `reviewed_uniprots` | UniProtKB | `>>hgnc>>uniprot` — **all** reviewed Swiss-Prot products |
| `canonical_uniprot` | UniProtKB | `reviewed_uniprots[0]` |
| `canonical_transcript` | MANE | `>>ensembl>>refseq[is_mane_select==true]` → `>>refseq>>transcript` |
| `entrez_id`, `ncbi_summary` | NCBI | `>>hgnc>>entrez` + entry summary |
| `clingen_dosage`, `depmap` | ClinGen, DepMap | `>>hgnc>>clingen_dosage`, `>>hgnc>>depmap` |

Two of these choices matter downstream:

- **`reviewed_uniprots` is a *tuple*, not a single accession.** A "gene" can
  encode more than one reviewed protein — *CDKN2A* makes both p16^INK4a
  (P42771) and p14^ARF (Q8N726) from alternative reading frames. Carrying all
  reviewed products lets the protein, structure, pathway and function sections
  union across them instead of silently dropping one.
- **The canonical transcript is MANE-Select**, picked with an explicit
  `[is_mane_select==true]` filter rather than by scanning RefSeq, so it is
  reproducible and never lost on genes with hundreds of RefSeq rows. It is the
  coordinate axis for AlphaMissense.

---

## 2. The eight zones

The published page is a frozen sequence of eight `## H2` zones (the contract is
in [PAGE_CONTRACT.md](PAGE_CONTRACT.md)). Each zone always appears — with an
informative placeholder if empty — so the table of contents is identical on
every gene page and deep links never break. The twelve collectors map onto the
zones like this:

| Zone `{#anchor}` | Built from | What it answers |
|---|---|---|
| **Summary** `{#summary}` | declarative lead + NCBI summary + JSON-LD | the one-line "what is this gene" |
| **Identifiers** `{#identifiers}` | §1 gene IDs | core IDs + an xref-count census |
| **Gene structure** `{#gene-structure}` | §2 transcripts, §11 expression, §9 regulation, §5 orthologs, derived genomics/GeneRIFs | the locus: isoforms, where it's expressed, how it's regulated |
| **Protein** `{#protein}` | §3 protein IDs, §4 structure, residue map | the product(s): domains, features, structures |
| **Function** `{#function}` | §7 pathways, §8 interactions | pathways, GO, the interaction network |
| **Disease & clinical** `{#disease}` | §6 variants, §12 associations, cancer significance | variants, disease links, cancer-driver role |
| **Drugs & pharmacology** `{#drugs}` | §10 drugs | targets, molecules, pharmacology, trials |
| **Related Atlas pages** `{#related}` | cross-entity [mesh](mesh.md) | links to associated diseases & drugs |

Sub-sections are demoted to `### H3` under their zone and carry **backend-owned
`{#anchor}` ids** — never Hugo's prose-derived ids, which would break a deep
link the moment a count in the heading changed.

### Summary
The lead is a *deterministic sentence* assembled from the bundle (no model): HGNC
name, location, gene class, the canonical protein and its UniProt FUNCTION first
sentence, the top CIViC verdict, and a DepMap/ClinGen dependency clause when
notable. Then the NCBI RefSeq curated summary, an "at a glance" digest, and the
inline schema.org JSON-LD (a `Gene` node encoding one typed `Protein` per
reviewed accession — see [how-it-works.md](how-it-works.md#schemaorg-json-ld)).

### Identifiers (§1)
Core HGNC fields plus an **xref-count census** read from the HGNC entry's xref
table — see the "count trick" in §4 below. RNAcentral is added only for
non-coding genes. Chains: `>>hgnc>>mim`, `>>hgnc>>entrez`, `>>hgnc>>rnacentral`.

### Gene structure
- **§2 transcripts** — Ensembl/RefSeq/CCDS transcript IDs and the MANE-Select
  exon structure (`>>ensembl>>transcript`, `>>transcript>>exon`). Transcript
  counts are marked *shrinkable*: NCBI periodically re-applies REVIEWED-only
  filtering (TP53 went 46→25), so a shrink is treated as drift, not regression.
- **§11 expression** — Bgee per-tissue scores, FANTOM5 CAGE promoters, and
  Single-Cell Expression Atlas (`>>ensembl>>bgee`, `>>ensembl>>fantom5_gene`,
  `>>ensembl>>scxa`). Anatomy IDs are resolved to names; only the top tissues
  are tabled.
- **§9 regulation** — the transcription-factor layer from CollecTRI, with the
  symbol interpolated into the filter: `>>hgnc>>collectri[tf_gene=="TP53"]` for
  downstream targets and `[target_gene=="TP53"]` for upstream regulators, plus
  JASPAR motifs and miRDB miRNAs. A gene is flagged a TF when it has downstream
  targets or JASPAR motifs.
- **§5 orthologs** — Ensembl Compara orthologs and paralogs
  (`>>ensembl>>ortholog`, `>>ensembl>>paralog`), emitted as returned (no species
  filter — see the correction in §4).
- **Derived**: gene-level ClinGen dosage + DepMap fitness (shown only when
  notable), and GeneRIF literature claims with PMID cites.

### Protein
Elides entirely for non-coding genes (placeholder: *"Non-coding RNA — no protein
product"*).
- **§3 protein IDs** — UniProt accessions, curated CC narratives (FUNCTION,
  SUBUNIT, SUBCELLULAR LOCATION, …), InterPro/Pfam domains, isoforms (this is
  what surfaces p53α/β/γ and p16/p14ARF), antibody resources, and BRENDA enzyme
  kinetics. Chains include `>>ensembl>>uniprot`, `>>uniprot>>interpro`,
  `>>uniprot>>pfam`, `>>uniprot>>ufeature`, `>>uniprot>>brenda`.
- **§4 structure** — PDB experimental structures and AlphaFold models
  (`>>uniprot>>pdb`, `>>uniprot>>alphafold`), unioned over all reviewed products.
  AlphaFold rows carry a `present` flag with a footnote for proteins too large
  to model (>~2700 aa: ATM, BRCA2, DMD, TTN, MUC16) — honest absence, not a
  broken link.
- **Residue map** (derived) — a render-only regrouping of UniProt sequence
  features into a drug-discovery view (active/ligand-binding sites, PTMs,
  disulfides, mutagenesis-validated residues), per product, each anchored
  `#protein-<accession>` to match the JSON-LD `@id`.

### Function
- **§7 pathways** — Reactome and GO, **unioned across every reviewed UniProt and
  the Ensembl gene route**, because per-product annotation diverges (CDKN2A's
  p14ARF mitophagy GO terms are absent from the p16/Ensembl route), plus MSigDB
  gene sets and GO/Reactome parent rollups.
- **§8 interactions** — STRING, IntAct, BioGRID, SIGNOR, and ESM2/Diamond
  structural similarity. The chains query the *interaction record*
  (`>>uniprot>>string_interaction`), not `>>…>>uniprot` (which collapses to bare
  partner IDs and loses the scores). True totals come from the entry xref count
  (TP53's 14,764 STRING partners: ~3 calls, not ~150 pages).

### Disease & clinical
Elides for non-coding genes.
- **Cancer significance** (derived) — folds the two best-grounded cancer signals
  (the CIViC gene narrative and the intOGen driver role) into one block; empty
  for non-cancer genes.
- **§6 variants** — ClinVar breakdown by germline classification
  (`>>hgnc>>clinvar[germline_classification=="Pathogenic"]`, one chain per
  class), SpliceAI, AlphaMissense **on the canonical transcript**, and a dbSNP
  sample (routed via Entrez because the direct `hgnc>>dbsnp` edge is unbacked).
  Per-class counts are labelled "(floor)" since they are pagination floors.
- **§12 associations** — the Mendelian/phenotype/complex-disease layer: OMIM,
  GenCC, MONDO, Orphanet, HPO, GWAS, EFO, MeSH, intOGen, CIViC, ClinGen
  Gene-Disease Validity. Disease-phenotype MIMs are reached through MONDO
  (`>>hgnc>>clinvar>>mondo>>mim`), not the gene's own `>>mim`. GenCC is collapsed
  to one row per disease, keeping the strongest classification.

### Drugs & pharmacology (§10)
The largest section — targets, molecules, pharmacology, pharmacogenomics, and
clinical evidence. Key curation:
- Molecules are filtered to **phased drugs only**
  (`>>chembl_target>>chembl_molecule[highestDevelopmentPhase>=1]`) — the raw edge
  returns thousands of screening compounds — and sorted by phase.
- **Curated affinities (Tier 1)** come from GtoPdb/IUPHAR; **BindingDB (Tier 2)**
  is human-filtered, normalised to nM, and carries a heterogeneity note.
- **Clinical trials are reached through the disease route**
  (`>>hgnc>>gencc>>mondo>>clinical_trials`), never `chembl_molecule>>
  clinical_trials`, which would pull in off-target drugs.

### Related Atlas pages
The human-facing cross-entity links (associated diseases, drugs that target the
gene). Built entirely from already-collected curated edges — see [mesh.md](mesh.md).

---

## 3. Non-coding genes

For a gene whose biotype is not `protein_coding`, the pipeline actively
**scrubs** the positionally-inherited blocks. biobtree links many datasets by
genomic overlap, so without this an antisense gene like *TTN-AS1* would inherit
*TTN*'s 6,528 ClinVar variants — confidently wrong, which is worse than empty.
The variant, disease and trial keys are cleared, and the Protein / Function /
Disease / Drugs zones collapse to their placeholders. A non-coding page is
honestly thin, not falsely full.

---

## 4. Why the data is trustworthy

These are the decisions that separate "scraped a lot of fields" from "got it
right":

- **Two-tier data rule.** Consolidated, comparable data leads plainly (UniProt
  features, GtoPdb affinities, PDB/AlphaFold). Heterogeneous data follows *with
  a provenance note doing real work* — BRENDA Km renders as a range with an
  "aggregated across organisms" caveat; BindingDB shows measured/human/total
  counts with a "heterogeneous assays — not directly comparable" note. The hard
  line: **wrong ≠ heterogeneous.** Contaminated edges are fixed or omitted (the
  disease-route trials above), not surfaced with a disclaimer.
- **Dual-product genes** union across all reviewed proteins (§3/§4/§7/§8), so a
  gene like *CDKN2A* shows both products' domains, structures and pathways
  instead of silently keeping only the first.
- **The xref-count trick.** Exact per-class totals (ClinVar, SpliceAI, MSigDB,
  HPO, GWAS, STRING, IntAct…) are read from an entry's xref table in *one call*
  rather than paginating thousands of rows — both faster and immune to the
  100-row pagination cap.
- **Column-shift guards + universal unescaping.** Every biobtree response is
  HTML-unescaped once at the single point all responses flow through (so prose
  reads `3'-end`, not `3&apos;-end`), and any row whose field count ≠ the schema
  column count is dropped — an unescaped pipe shifting columns is the bug class
  that silently corrupts every later value.
- **Notability gates.** DepMap is only shown when the gene is actually a
  dependency (≥10% dependent, or strongly-selective/common-essential); a 0.2%
  line would be a low-value orphan.

> **Correction worth recording for the preprint:** there is *no* human-only
> filter on orthologs — §5 emits Ensembl Compara orthologs/paralogs as returned.
> The human-only organism filter applies to **BindingDB** and **CTD**
> interactions, not orthologs.
