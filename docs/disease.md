# How a disease page is built

This walks through how Atlas turns a disease into a published page — how the
disease corpus is built, how a disease resolves to a gene cohort, and how each
section is mined. For the shared machinery (transport, the collect→render split,
the page contract, JSON-LD), see [how-it-works.md](how-it-works.md); this doc is
disease-specific.

The disease page is distinctive in one way: most of its biology comes from
**reusing the gene collectors over a cohort of associated genes**, then
aggregating.

```
Mondo OBO ──► corpus build (rank + gate) ──► seed list (by MONDO id)
   │
   ▼  per disease:
resolve(MONDO id) ──► DiseaseAnchors  (IDs + a ranked 4-route gene cohort)
   │
   ▼
14 collectors  (disease-level chains + cohort fan-out over the gene collectors)
   │
   ▼
render_all ──► 6 canonical H2 zones + derived views
   │
   ▼
assemble_page ──► frontmatter · ## Summary · body · ## Related Atlas pages
```

---

## 1. Building the disease corpus

There is no "list all diseases" call, so the corpus is built from the Mondo
ontology (`disease/corpus.py`):

1. **Parse** the Mondo OBO (cached locally), keeping `{id, name, parents}` and
   dropping obsolete terms.
2. **Admission gate.** Remove the disease-*characteristic* qualifier subtree
   (rooted at `MONDO:0021125`: *inherited / acquired / sporadic / congenital /
   X-linked* — these are qualities, not diseases, but they word-match a lot of
   evidence and otherwise float to the top of the ranked corpus). The gate is a
   pure graph walk and runs *before* probing, so it both cleans the corpus and
   saves the probes. It is deliberately narrow — it does **not** gate the broad
   `disease_grouping` subset, which contains real hub diseases (cardiomyopathy,
   AML, renal cell carcinoma).
3. **Score + rank.** One `entry` call per term reads its xref counts; a weighted
   `signal_score` favours curated over literature-mined over raw evidence
   (`gencc` 3.0, `civic_evidence`/`intogen` 2.0, `gwas`/`clinvar`/`hpo` 1.0), with
   `clinical_trials` damped (÷100) so trial-heavy diseases don't drown out
   mechanistically rich ones. Terms are ranked into tiers; zero-signal terms are
   dropped.

The seed list is emitted **by MONDO id**, so the per-disease build passes the id
to `resolve()` — biobtree's name search ranks by xref count, so `cardiomyopathy`
would otherwise resolve to the high-xref subtype *dilated cardiomyopathy 1G*
rather than the umbrella term.

---

## 2. Resolving the disease and its cohort

`resolve(MONDO id)` runs once and pre-computes everything the 14 sections need:

- **Federated IDs** — sibling ontology ids via one-hop maps (`>>mondo>>efo`,
  `>>mondo>>mesh`, `>>mondo>>mim`, `>>mondo>>orphanet`) and cross-ontology xrefs
  (DOID, UMLS, NCIT, ICD, GARD, MedDRA…) fetched *only when* the entry's xref
  count for that dataset is non-zero — saving about half the round-trips. UBERON
  anatomy (gated the same way) drives the JSON-LD `associatedAnatomy`.
- **Orphanet entry** — the first Orphanet record of type *Disease* (not a
  clinical subtype or group), which carries HPO phenotypes and prevalences.
- **The gene cohort** — the heart of the page.

**The gene cohort.** Four chains, each tagged with an evidence flag, are unioned:

| Route | Chain | Evidence |
|---|---|---|
| GWAS | `>>mondo>>gwas>>hgnc` | common-variant |
| GenCC | `>>mondo>>gencc>>hgnc` | curated gene-disease |
| ClinVar | `>>mondo>>clinvar>>hgnc` | germline variant |
| CIViC | `>>mondo>>civic_evidence>>hgnc` | somatic |

Genes are ranked by **how many evidence routes they hit** (ties by id), so
dual/triple-evidence genes survive the cap of 50. Each cohort gene is pre-resolved
to a gene anchor record, so the cohort sections can reuse the existing gene
collectors without re-paying the anchor cost. (The cap bounds wall-clock: ~50
genes × ~12 collectors × ~5 chains ≈ 3,000 biobtree calls per disease.)

---

## 3. The six zones

Fourteen collectors plus three render-only derived views fold into six frozen
`## H2` zones plus Summary and Related (see [PAGE_CONTRACT.md](PAGE_CONTRACT.md)).

| Zone `{#anchor}` | Built from | What it answers |
|---|---|---|
| **Summary** `{#summary}` | declarative lead + at-a-glance + JSON-LD | the one-line "what is this disease" |
| **Identifiers** `{#identifiers}` | §1 IDs | federated ids, synonyms, epidemiology, phenotypes |
| **Genetics & variants** `{#genetics}` | §2 GWAS, §3 variant tiers, §4 Mendelian/somatic overlap | the variant landscape |
| **Genes & proteins** `{#genes}` | §5–§9 cohort fan-out | the associated-gene cohort, profiled |
| **Function** `{#function}` | §14 pathways | pathways the cohort touches |
| **Therapeutics** `{#drugs}` | §10–§12 + derived tractability views | drugs against the cohort, repurposing signals |
| **Clinical trials & evidence** `{#trials}` | §13 trials + CIViC | trials for the disease, precision-oncology map |
| **Related Atlas pages** `{#related}` | cross-entity [mesh](mesh.md) | linked genes & drugs |

### Summary
A deterministic lead — *"{Name} (MONDO:id) is a {cancer|disease} with N cohort
genes (… GWAS associations; … CIViC somatic drivers; … ClinVar records) and T
clinical trials…"* — then the at-a-glance digest and inline JSON-LD. Two numbers
here are deliberately sourced carefully: the trial count is the **validated** §13
count (never the contaminated raw xref count), and the ClinVar count is the
accurate xref total (not §3's paginated floor).

### Identifiers (§1)
The federated identifier table, synonyms, a data-availability line (excluding
identity/ontology xrefs and the contaminated raw trial count), epidemiology
(prevalences, validated-first), and HPO clinical features. No new calls — pure
shaping of the anchor record.

### Genetics & variants
- **§2 GWAS** — totals, top hits by p-value, top studies (`>>mondo>>gwas`,
  `>>mondo>>gwas_study`). Odds ratios round through `fnum`; **p-values are not
  rounded** (rounding `3e-67` would collapse it to 0).
- **§3 variant tiers** — GWAS variants tiered (coding/splice/regulatory/intronic)
  via dbSNP, plus ClinVar germline (`>>mondo>>gwas>>dbsnp`, `>>mondo>>clinvar`).
  The fetched `clinvar_total` is presented as a *floor* so it never contradicts
  the accurate xref total.
- **§4 Mendelian & somatic overlap** — which cohort genes also cause Mendelian
  forms (GenCC/Orphanet/OMIM) or carry somatic-driver evidence (intOGen/CIViC, for
  cancers). GWAS + Mendelian = highest confidence. GenCC is deduped to one row per
  gene (a cohort fan pulls every submission — Lynch syndrome had MLH1 ×19), with
  the two-tier rule preferring the record *for this disease* before keeping the
  strongest classification (so BRCA2's Fanconi-D1 record can't outrank its
  on-disease record on a different cancer's page).

### Genes & proteins — the cohort fan-out
The defining mechanism: each of §5–§9 runs an existing **gene** section collector
over every cohort gene and aggregates the results — the disease sections add no
new chains of their own, they reuse the gene mining:

- **§5 genes→proteins** — gene IDs + protein IDs; cohort partitioned by evidence
  bucket (GWAS-only, GWAS+GenCC, multi-evidence…).
- **§6 protein families** — InterPro classification → druggable / difficult /
  unknown split; the family distribution is **ORA-ranked** (enriched families such
  as Kinase rise; the catch-all Other/Unknown sinks despite the largest raw count).
- **§7 expression context** — Bgee/single-cell, bucketed by expression breadth.
- **§8 interactions** — the intra-cohort interaction graph and hub genes.
- **§9 structural data** — PDB / AlphaFold coverage across the cohort.

### Function (§14)
Reactome pathways and GO biological processes the cohort touches, ranked by
**over-representation analysis** (ORA): a hypergeometric test of the cohort's
overlap against a genome-wide background, Benjamini-Hochberg FDR, sorted by
enrichment — so size-biased umbrella categories (Signal Transduction, Immune
System) no longer float to the top by raw count, and disease-specific pathways
surface. The raw cohort-gene count and gene members are kept (consumer
ground-truth); a fold-enrichment + FDR column is added. The background sizes are
counted via the SAME chain/classifier the cohort uses over all protein-coding
genes (`atlas.build_background` → `data/background/{reactome,go,family}.json`,
refreshed per biobtree data release); the statistics are pure + deterministic
(`atlas.ora`, no scipy).

### Therapeutics
- **Disease-direct indicated drugs (`#indicated`)** — drugs with a registered
  ChEMBL indication *for this disease*, independent of the gene cohort (so
  cohort-less / autoimmune conditions still surface real drugs). Tiered so it
  never overstates: **approved** drugs (phase 4, or an anticancer drug at phase 3
  vs this cancer — see [drug §4](drug.md)) tabled as indicated; phase 2–3 listed
  separately as **in clinical trials** (investigational — a trial record, not an
  indication). The `approved` decision is the same flag the drug §4 collector
  computes, so the two pages can't disagree.
- **§10 drug-target analysis** — fans the gene drug collector over the cohort:
  per-gene max development phase, approved/phased/undrugged buckets, and a
  cross-gene drug aggregation (deduped by molecule).
- **§11 bioactivity & enzyme** and **§12 pharmacogenomics** add screening depth
  and PharmGKB/CPIC/DPWG coverage (the upstream `is_vip` flag is dropped — it is
  degenerate, true for every record).
- **Derived views** — repurposing candidates, a druggability pyramid (tiers A–E
  from family + structure + phase), and undrugged-target profiles. The repurposing
  framing is gated on disease class: for cancers it's "repurposing candidates";
  for non-cancer/Mendelian diseases it's explicitly "chemical tractability — a
  research signal, **not** a therapeutic recommendation", because a cohort-gene
  bioactivity hit is often off-target screening binding.

### Clinical trials & evidence (§13)
Disease-level trials (not a cohort fan-out), the drugs tested, and the CIViC
precision-subtype map. The key correctness step is **title-validation**:
biobtree's `mondo→clinical_trials` edge is contaminated (*Vici syndrome* →
1,156 glaucoma/cataract trials), so a trial is kept only if its title names the
disease's canonical name or a Mondo synonym (tokens ≥4 chars, to drop ambiguous
short acronyms). The validated count is authoritative; the raw count is shown only
as a disclosure. Trial drugs fold salt-forms into parents; only drugs in ≥1
validated trial are kept. The CIViC map keeps predictive evidence, deduped and
ranked by level.

---

## 4. Why the data is trustworthy

- **Curated provenance over inference.** The cohort is built from curated/typed
  evidence routes; disease→drug edges come from *trials + CIViC*, never bioactivity
  inference — the same discipline as the gene and drug pages.
- **Title-validated trials.** The contaminated trials edge is filtered to trials
  that actually name the disease; the validated count is the one used everywhere
  (lead, glance, distribution).
- **Accurate counts over paginated floors.** ClinVar and trial counts come from
  the authoritative xref total, with paginated fetch counts labelled "floor".
- **Two-tier GenCC** (on-disease record preferred, then strongest class) and
  one-row-per-gene dedup keep a polygenic cohort honest.
- **Class-gated framing.** Repurposing is only called "repurposing" for cancers;
  elsewhere it's labelled a research signal, not a recommendation.
- **Suppressing noise.** Raw ontology accessions leaking where a label belongs
  (including `MP:` mouse-phenotype ids) are dropped, not shown as diseases; the
  degenerate PharmGKB VIP flag is dropped entirely.
