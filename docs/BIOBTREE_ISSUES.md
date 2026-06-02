# biobtree MCP — Issues & Improvement Requests

**Initial filing:** 2026-05-28
**Last updated:** 2026-06-01 (biobtree refresh resolved #12, #13, #14, #15, #18 + Mondo OBO xrefs; speculative asks #21/#22/#23/#24 removed; filed patent gaps #25/#26/#27; filed trial-edge contamination #28)

**Context:** Found while building a deterministic gene/disease reference-page
collector (Sugi Atlas) on top of the local biobtree REST API
(`http://127.0.0.1:8000/api`, with `/search`, `/entry`, `/map`) and the MCP
server (`biobtree_search` / `biobtree_map` / `biobtree_entry`).

Numbers are stable — when an issue resolves it stays here as a one-line entry
under "Resolved" so commit/PR references against the original number stay
unambiguous. Bodies are removed once resolved to keep the doc compact.

---

## Resolved upstream

| # | Title | Resolved |
|---|---|---|
| #1 | `not_found` overloaded — same response for "unknown input" vs "input recognized, chain dead-ends" | 2026-05-30 — now returns `message: "Recognized but no mapping for this chain: <input>. Try a shorter chain."` |
| #2 | Bare 400 on unknown dataset name | 2026-05-30 — now returns `"unknown dataset: 'omim'. Did you mean 'mim'?"` |
| #3 | Statless empty blob on certain bad chains | 2026-05-30 — consistent shape with explicit `message` field on every outcome |
| #5 | `map` source-dataset semantics non-obvious (e.g. `map(HGNC:11998, ">>ensembl")` silently empty) | 2026-05-30 — same signaling work as #1/#3 |
| #7 | `map` mixes species, no human filter | 2026-05-30 — `[genome=="homo_sapiens"]` chain filter now supported |
| #8 | `map` pagination "broken" — was a wrong param name (`p` not `page`); FastAPI silently dropped the unknown one | resolved 2026-05-28 (caller side) |
| #11 | `>>uniprot>>ufeature` leaked ortholog features (P04637 query returned mouse Tp53 rows) | 2026-05-30 — single-accession queries now return only that accession's features. Atlas's `id.startswith(u + "_")` post-filter removed |
| #9 | UniProt entry payload too thin (no CC, no `reviewed`, no isoforms) | **2026-05-31 — RESOLVED.** Entry now carries `comments` (function, subunit, subcellular_location, tissue_specificity, disease, ptm, cofactor, domain, induction, miscellaneous, similarity, caution) + `isoforms` list with named-isoform IDs (e.g. `P04637-1..9` for p53α/β/γ, Del40, Del133) + `is_canonical` flag on the principal isoform. Unreviewed accessions return empty attrs (implicit Swiss-Prot filter). |
| #10 | `>>uniprot>>alphafold` empty for very large proteins | **2026-05-31 — RESOLVED biobtree-side.** Coverage extended; MTOR (2549 aa) and similar now have data. The remaining empties (ATM/BRCA2/DMD/TTN/MUC16 — all > ~3000 aa) reflect that **AlphaFold DB upstream genuinely has no model** for these. Atlas's responsibility: render a graceful "AlphaFold DB does not provide a model for proteins > ~3000 aa" footnote when alphafold is empty for an extra-large protein. |
| #12 | `pubchem_activity` empty for KRAS (under upstream investigation) | **2026-06-01 — RESOLVED.** KRAS `>>uniprot>>pubchem_activity` returns 4927 rows. Atlas's chembl_activity fallback comment removed in commit 66ed71b. |
| #14 | Reactome pathway entries with empty `name` field | **2026-06-01 — RESOLVED.** Every probed pathway carries a name. Atlas's "Unnamed pathway (R-HSA-N…)" graceful fallback removed in commit 66ed71b. |
| #15 | `chembl_molecule` parent/child salt-form linkage | **2026-06-01 — RESOLVED.** Both directions work: `>>chembl_molecule>>chembl_moleculeparent` (child→parent, e.g. CHEMBL1642 Imatinib mesylate → CHEMBL941 Imatinib) and `>>chembl_molecule>>chembl_moleculechild` (parent→children). Unblocks drug entity at scale. |
| #18 | GtoPdb drug→target: interaction id substring contamination | **2026-06-01 — RESOLVED.** Accessible via `>>chembl_molecule>>gtopdb_ligand>>gtopdb_interaction` (3 hops). |
| — | Mondo OBO cross-ontology xrefs + UBERON anatomy | **2026-06-01 — RESOLVED.** `>>mondo>>{doid,sctid,umls,ncit,medgen,icd10cm,icd11,gard,meddra,nord,uberon}` all work. §1 federated identifier table extended; JSON-LD `sameAs` + `code` + `associatedAnatomy` populated (commit d911cb9). |
| #13 | `pharmgkb_guideline` / `_clinical` / `_variant` edges empty | **2026-06-01 — RESOLVED.** All three populated for pharmacogenes; the earlier "still empty" probe was a false negative against a non-pharmacogene (KRAS). CYP2C19 returns 37 guidelines, CYP2D6 69, TPMT/DPYD/VKORC1/SLCO1B1 6-14 each. Atlas §10 now renders all three blocks. |

---

## Open

## Issue #25 — Patent attributes documented but not populated (assignee / CPC / IPC)

`docs/datasets/patent.md` lists `asignee`, `cpc`, `ipcr`, `ipc` as stored
`PatentAttr` fields, and three of the doc's seven headline use cases
(competitive intelligence via assignee, freedom-to-operate + technology
landscaping via CPC/IPC) depend on them. In the current build none are
populated.

**Repro (2026-06-01):**
```
entry("EP-2914622-B1", "patent")  # EP, not just CN
→ Attributes.Patent keys: [title, country, publication_date, family_id, id]
   (no asignee, no cpc, no ipc, no ipcr — checked across CN/EP/US/WO samples)
```

**Atlas impact:** blocks the assignee breakdown ("who holds the IP") and the
technology-classification landscape (CPC/IPC) we'd surface on drug pages.
Only `title / country / publication_date / family_id` are available, so the
drug Patent section stays a coarse count + per-compound split.

---

## Issue #26 — `patent_compound` exposes a `patent` xref count but no `patent_family` rollup

The meaningful dedup of a patent footprint is *distinct families* (one
invention, many jurisdictions), not raw mention count. `patent` entries carry
`family_id`, and `patent>>patent_family` exists per-patent, but `patent_compound`
only exposes a `patent` xref — there's no `patent_compound>>patent_family` edge
or family count.

**Repro (2026-06-01):**
```
entry("3827", "patent_compound").xrefs
→ {pubchem:1, chebi:1, chembl_molecule:1, hmdb:1, patent:107733}   # no patent_family
```

So a drug's distinct-family count is only reachable by `entry()`-ing all N
patents for their `family_id` (107k+ for Imatinib's primary compound) —
infeasible. A `patent_family` xref on `patent_compound` (or a family-count
attribute) would make "N inventions across M jurisdictions" a one-call signal.

**Atlas impact:** can't show distinct patent families (the honest dedup metric);
forced to report raw SureChEMBL mention counts, which the doc itself warns are
inflated by promiscuous compounds.

---

## Issue #27 — `>>patent_compound>>patent` is ID-ordered, no date sort / no aggregate facets

`map(... >>patent_compound>>patent)` returns patents in patent-ID order, not by
date. A bounded sample is therefore unrepresentative — e.g. the first pages of
Imatinib's primary compound (107k patents) are all `CN-100…` from 2007–2015,
while a smaller compound's first pages are EP/US 2022–2025. Country/year/
"recent patents" computed from any sample is a sampling artifact, not a real
landscape.

**Repro (2026-06-01):**
```
map("3827", ">>patent_compound>>patent", cap=2)   # 300-row sample
→ 100% CN, years 2007–2015  (id-ordered; not the true jurisdiction/time mix)
```

**Suggested fix:** a date-sorted option (`order=publication_date desc`) for
"recent patents", and/or server-side facet counts (by country / year / CPC)
so a drug-patent landscape can be shown without enumerating 100k+ rows.

**Atlas impact:** can't honestly surface jurisdiction breakdown, filing
timeline, or "most recent patents" — only the accurate total + per-compound
split (which need no enumeration).

---

## Issue #28 — `mondo>>clinical_trials` edge is contaminated (unrelated trials)

`map(<mondo>, ">>mondo>>clinical_trials")` returns trials whose actual
conditions don't match the disease, with absurd counts on rare diseases.

**Repro (2026-06-02):**
```
map("MONDO:0009452" /* Vici syndrome, prevalence <1/1,000,000 */, ">>mondo>>clinical_trials")
→ 1,156 trials; sampled conditions are Glaucoma / Cataract
  (NCT00273221 "Phacotube vs Phacotrabeculectomy", NCT00312299 "Posterior Capsule
   Opacification Study") — none about Vici syndrome.
```
1,156 > cardiomyopathy's 317 for an ultra-rare disorder is the tell.

**Root cause (confirmed by biobtree dev, 2026-06-02):** the `clinical_trials→mondo`
edge is built by `collectOntologyIDs()→d.lookup(condition)`, which greedily adds
every ontology id a token/text lookup returns. Vici's Mondo synonyms include
"absent corpus callosum **cataract** immunodeficiency", so a trial with condition
"Cataract" resolves to Vici (and every rare syndrome whose multi-symptom synonym
mentions a common term). **Scope is bigger than trials:** `collectOntologyIDs`
is shared by **clinical_trials, intogen, AND civic** → the same over-linking
contaminates `mondo→intogen`, `mondo→civic_evidence`, and `mondo→clinical_trials`.
**What Atlas actually traverses:** only `mondo→clinical_trials` (§13) and
`mondo→civic_evidence` (§13 subtype map + disease→drug mesh + cohort CIViC route).
Atlas reaches intogen and civic *drivers* gene-first (`hgnc→intogen`, `hgnc→civic`
in §4) — the gene-side edges are keyed by gene and are NOT affected. So the
**mondo→intogen edge is never traversed by Atlas (blast radius 0)**; the only
contaminated edge in use is `mondo→civic_evidence`.

**Upstream fix (dev, in progress — needs a re-index):** require an EXACT
name/synonym match in `collectOntologyIDs` (condition == a full MONDO name or
synonym), not a token/text hit. One change fixes trials + intogen + civic; only
takes effect on a re-index.

**Atlas mitigation (commit 43c5736):** s13 title-validates trials (keep only those
whose brief_title names the disease/synonym) → trial_count/drugs/lead/At-a-glance
use the validated set. Interim only; brief_title is a too-strict proxy
(under-counts: cardiomyopathy 317→92). **No mitigation added for civic** — see
blast-radius result below.

**Blast-radius check on `mondo→civic_evidence` (curated head, 268 dense-test
diseases, 2026-06-02):** only 18/268 pages carry §13 CIViC rows at all (303 rows);
**97% (294/303) match the page's canonical disease name**, 16/18 pages 100% clean.
The 2 pages with non-canonical rows are both legitimate, NOT contamination:
`glioblastoma` (27 Glioblastoma + 3 *Glioma*, the parent class) and
`tumor-predisposition-syndrome-3` (6 *Glioma* rows — its own synonyms are "glioma
susceptibility 9" / "malignant glioma caused by mutation in POT1", so it genuinely
predisposes to glioma). **Conclusion: ~0 genuine contamination on the head.** The
Vici-style over-link is a long-tail rare-disease phenomenon that doesn't reach the
curated head, so **no interim civic guard is warranted — wait for the upstream
re-index.** (Re-run this check if the corpus expands deep into the rare-disease
tail before #28 lands.)

**⚠ REVERT-ON-RESOLVE:** when the re-index lands —
1. swap the s13 title-match → **exact condition-match** (self-deactivating: once
   the edge only returns condition-matching trials it becomes a no-op) — requires
   biobtree adding `conditions` to the clinical_trials map `compact_fields` (dev
   offered); OR remove the guard entirely and trust the edge;
2. re-check disease §4 (intogen) / §13 (civic) / the cohort CIViC route for residual
   contamination and drop any guards added there.

---

## Issue #4 — Long multi-hop chains fail silently / intersectively

Every individual hop can succeed while the combined chain returns empty,
with no signal indicating *which* hop broke the chain.

**Repro (still relevant 2026-05-30):**
```
map("TP53", ">>hgnc>>ensembl>>entrez>>mim>>refseq")
→ message: "Recognized but no mapping for this chain: TP53. Try a shorter chain."
```

After the #1/#3 signaling work, the caller now learns the chain produced
nothing, but still can't tell *which hop* dropped to zero (each individual
hop succeeds; some pair in the middle has no edge between them).

**Suggested fix:** return per-hop diagnostics — either a partial result up
to the break with a note of where it stopped, or a `chain_diagnostics`
field that says "hop 4 (mim→refseq) returned 0 rows". Either makes the
debugging path orders of magnitude faster.

---

## Issue #6 — `entry` xrefs return counts, not values

`entry` shows how many links exist to each dataset, but not the linked IDs.

**Repro:**
```
entry("HGNC:11998", "hgnc")
→ xrefs: ["ensembl|1", "entrez|1", "mim|1", "refseq|93", ...]   (counts only)
```

So you know an OMIM/Ensembl/Entrez link exists, but must issue a separate
`map` per dataset to get the actual value. There is no single
"give me all the cross-references" call, which multiplies round-trips and
the chance a weak model drops a field.

**Suggested fix:** an option to resolve xrefs to their target IDs
(`entry(..., expand_xrefs=true)`), or a convenience that returns all
standard gene IDs in one call.

**Atlas impact:** non-blocking — Atlas enumerates IDs via `map_all` per
chain. Quality-of-life win for weak callers.

---

## Retracted

| # | Title | Reason |
|---|---|---|
| #16 | No `list-ids` endpoint for a dataset (corpus enumeration) | Retracted 2026-05-31. Corpus enumeration is naturally upstream of biobtree — HGNC ships `hgnc_complete_set.txt` (all human genes), Mondo ships its OBO at obofoundry.org (all disease classes), ChEMBL ships SQLite/flat-file releases. biobtree itself consumes these source files; Atlas can parse the same files for "discover the corpus" without asking biobtree to duplicate that role. Filing was a "biobtree is single source" aesthetic, not a real engineering need. |
| #17 | No bulk xref-count check | Retracted 2026-05-31. The mondo entry already exposes xref counts; calling `/entry` per id is N round-trips but only once per release, easily cached locally. At ~25k Mondo nodes this is ~30 min serial / a few min parallel — not prohibitive. The motivation was that #16 made enumeration cheap and #17 made filtering cheap; once #16 is retracted, #17 loses its main case. A future "batch entry resolve" endpoint might help at scale, but no current bottleneck justifies filing it now. |
