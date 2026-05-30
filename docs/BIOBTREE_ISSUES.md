# biobtree MCP — Issues & Improvement Requests

**Initial filing:** 2026-05-28
**Last updated:** 2026-05-30 (now includes disease-page findings: #14, #15 + scale-out requests #16, #17)

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

---

## Open

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

## Issue #9 — UniProt entry payload is too thin (no CC, no `reviewed`, no isoforms)

Today `GET /api/entry?i=<acc>&s=uniprot` returns only:
`names / alternative_names / sequence / id / name`.

The full UniProt flat-file (e.g. `https://rest.uniprot.org/uniprotkb/P04637.txt`)
carries far richer content that's missing here. Three categories of asks,
in priority order:

**(1) CC ("Comments") narrative blocks** — the curated free-text
descriptions UniProt is famous for: `FUNCTION`, `SUBUNIT`,
`SUBCELLULAR LOCATION`, `TISSUE SPECIFICITY`, `DEVELOPMENTAL STAGE`, `PTM`,
`DISEASE`, `DISRUPTION PHENOTYPE`, `INDUCTION`, `DOMAIN`, `MISCELLANEOUS`,
etc. This is *the* highest-leverage piece of UniProt for any AI/agent
consumer — the paragraph an AI quotes when answering "tell me about gene X"
almost always paraphrases UniProt's `FUNCTION`.

**(2) `reviewed` flag** — whether the accession is reviewed (Swiss-Prot)
vs unreviewed (TrEMBL). Currently the only way to tell is to use
`>>hgnc>>uniprot` as a workaround (HGNC xrefs only the curated set).

**(3) Named alternative products** (`ALTERNATIVE PRODUCTS` section) —
isoform names like p53α/β/γ for TP53, K-Ras4A/K-Ras4B for KRAS, etc.

**Suggested response shape:**
```json
{
  "Attributes": {
    "Uniprot": {
      "id": "P04637",
      "name": "P53_HUMAN",
      "reviewed": true,
      "names": [...],
      "alternative_names": [...],
      "sequence": "...",
      "comments": {
        "FUNCTION": "Multifunctional transcription factor that induces cell cycle arrest, DNA repair or apoptosis...",
        "SUBUNIT": "Forms homodimers and homotetramers. Interacts with MDM2...",
        "SUBCELLULAR_LOCATION": "Cytoplasm. Nucleus. Nucleus, PML body.",
        "TISSUE_SPECIFICITY": "Ubiquitous; isoforms are expressed in a wide range of normal tissues...",
        "PTM": "Phosphorylated on Ser-15 by ATM, ATR and DNA-PK...",
        "DISEASE": "Li-Fraumeni syndrome 1 (LFS1) [MIM:151623]: An autosomal dominant familial cancer syndrome..."
      },
      "isoforms": [
        {"id": "P04637-1", "name": "p53alpha", "synonyms": ["p53"], "is_canonical": true},
        {"id": "P04637-2", "name": "p53beta", "synonyms": [], "is_canonical": false}
      ]
    }
  }
}
```

Notes for the implementer:
- Headers normalized to underscore-snake-case (`SUBCELLULAR_LOCATION`).
- Evidence codes (`{ECO:0000269|PubMed:9774970}`) — strip on the way out
  for clean consumer-side JSON.
- biobtree already parses UniProt for the existing thin payload; surfacing
  the additional fields is incremental work on the same parser.

**Atlas impact:** **biggest single content gap.** Audit (`docs/research/
03_page_audit.md`) flagged "no curated narrative" as the #1 reason an AI
agent picks UniProt/NCBI over Atlas today. Path C item #1 in
`docs/research/NEXT.md` is paused on this. Partial mitigation now available
via CIViC gene-level paragraphs (recent biobtree addition), but only for
cancer-relevant genes — UniProt CC is still required for the long tail.

---

## Issue #10 — `>>uniprot>>alphafold` is empty for very large proteins

> **FIXED in code 2026-05-30 (commit `fb3f917`) — ships in prod next release
> (~2026-05-31).**
>
> *Correction to the report's original premise:* current AlphaFold DB (v6)
> does **not** ship fragmented `AF-<acc>-F1..Fn` models for the headline
> examples. ATM (Q13315), BRCA2 (P51587) and the DMD canonical sequence
> return **0 models** from the AlphaFold prediction API — AlphaFold dropped
> large-protein fragments from the SwissProt tar and most are absent
> entirely now. Of the ~1,243 reviewed proteins >2700 aa, only ~1/3 have
> any model at all (DMD only has models for shorter isoforms).
>
> **Fix:** an opt-in backfill (`largeProteinBackfill=yes` on the alphafold
> dataset) enumerates reviewed proteins >2700 aa via UniProt REST, queries
> the AlphaFold prediction API per accession, and ingests whatever models
> exist (canonical fragments and/or isoform models) as one aggregate
> `AlphaFoldAttr` per protein (length-weighted mean pLDDT, pooled
> confidence fractions, total residues, fragment count). Off by default
> (~1,243 live API calls); enable + re-run alphafold update to populate.

For proteins >~2700 aa (ATM Q13315, BRCA2 P51587, DMD P11532) the alphafold
map returns 0 rows.

**Atlas workaround (still in place):** construct `AF-<acc>-F1` for every
reviewed protein and attach pLDDT only when the map provides it. Will lift
to consume the new aggregate attribute when the next release lands.

---

## Issue #12 — `pubchem_activity` index is empty for KRAS (and similar targets)

> **FIXED in code 2026-05-30 (commit `fb3f917`) — ships in prod next release
> (~2026-05-31).**
>
> *Actual root cause (differs from the snapshot-subset hypothesis in the
> original filing):* the ingested `bioactivities.tsv.gz` leaves its
> **Protein Accession and Gene ID columns empty** for the vast majority of
> rows (verified: all of the first 200k rows), so the activity→protein edge
> was never built. The authoritative per-assay target mapping lives in a
> **separate file biobtree did not join** — `Aid2GeneidAccessionUniProt.gz`.
> KRAS proves it: gene 3845 has 256 assays, 216 mapped to P01116 there, but
> none in the activity rows themselves.
>
> **Fix:** join `Aid2GeneidAccessionUniProt.gz` keyed by AID (panel-limit 25
> to drop a handful of unresolvable mega-panels) to emit activity→uniprot;
> plus derive uniprot from the curated Entrez entry for rows that do carry
> a gene_id. Restores `>>uniprot>>pubchem_activity` for KRAS and ~369k
> assays overall.

KRAS via every Atlas-tried route returned n=0:
```
>>uniprot>>pubchem_activity                          n=0   (direct)
>>hgnc>>uniprot>>pubchem_activity                    n=0   (via hgnc)
>>hgnc>>entrez>>pubchem_activity                     n=0   (via entrez)
>>uniprot>>chembl_target>>pubchem_activity           n=0   (via chembl_target)
```

While `>>uniprot>>bindingdb` n=100 and `>>uniprot>>chembl_target>>chembl_activity`
n=100 both work fine for KRAS — confirming KRAS has curated binding data,
just not in biobtree's pubchem_activity slice.

**Atlas mitigation (kept in place):** §10 now wires `chembl_activity`
alongside `pubchem_activity` (commit `189f98a`). Both blocks render where
data exists; KRAS now has 5,239 ChEMBL activities (4,825 potent at
pChembl≥5). When the pubchem_activity fix lands the PubChem block will
backfill automatically.

---

## Issue #13 — `pharmgkb_guideline` / `pharmgkb_clinical` / `pharmgkb_variant` edges present but empty

These edges exist in `/api/help?topic=edges` (and were present in the
2026-05-30 refresh batch alongside the new datasets), but every query route
returns n=0, including for the canonical pharmacogenes where these
annotations are well-known to exist in PharmGKB.

**Repro (CYP2C19 = HGNC:2621):**
```
map(i="HGNC:2621", m=">>hgnc>>pharmgkb_guideline")  → n=0
map(i="HGNC:2621", m=">>hgnc>>pharmgkb_clinical")   → n=0
map(i="HGNC:2621", m=">>hgnc>>pharmgkb_variant")    → n=0
```

CYP2C19 has at least 18 CPIC dosing guidelines and dozens of clinical
annotations on pharmgkb.org. Other canonical pharmacogenes (CYP2D6, TPMT,
DPYD, VKORC1, SLCO1B1) are equally well-represented in upstream PharmGKB
and equally empty here.

`pharmgkb_gene` is populated and works fine (Atlas's §10 wires
`>>hgnc>>pharmgkb_gene`); `pharmgkb_pathway` is also populated. The
empties cluster on the **clinical / guideline / variant** trio.

**Hypothesis:** these three datasets' edges are declared but the ingest
either:
1. hasn't pulled the actual PharmGKB data files for these record types
   yet (the "annotated but not loaded" case), or
2. is keyed off a join that isn't matching — e.g. variant ↔ rsID where
   the rsID linkage from gene didn't index the variant side.

**Suggested fix:** confirm whether the PharmGKB
annotation/guideline/variant tables are in the ingest pipeline. If they
are, check the join key (CYP2C19 should match by `gene_id`, `symbol`, or
via uniprot P33261). If they aren't, schedule the ingest.

**Atlas impact:** PharmGKB clinical guidelines (CPIC dosing tables) are
the single most-cited PGx-decision source in clinical pharmacology.
Without them, Atlas's §10 PharmGKB block can only state "there's a
pharmgkb_gene entry"; the *contents* are unreachable. Will be more
pronounced once drug-pages are tackled.

## Issue #14 — Reactome pathway entries with empty `name` field

Some pathway records are indexed in biobtree (they're returned from
`>>uniprot>>reactome` / `>>ensembl>>reactome` map chains and have working
`/entry` endpoints) but their `name` attribute is empty.

**Repro:**
```
search(R-HSA-549132)  → R-HSA-549132|reactome||54    # name column blank
search(R-HSA-425366)  → R-HSA-425366|reactome||284   # name column blank
entry(R-HSA-549132, "reactome")  → Attributes.Reactome is empty
entry(R-HSA-425366, "reactome")  → Attributes.Reactome is empty
```

Both pathways are real Reactome entries
(`reactome.org/PathwayBrowser/#/R-HSA-549132` = "Synthesis of bile acids
and bile salts via 24-hydroxycholesterol"; `R-HSA-425366` = "Transport of
bile salts and organic acids, metal ions and amine compounds"). The
biobtree index has the id but is missing the human-readable label.

**Atlas impact:** Sugi Atlas's §14 cohort-pathway aggregation surfaces
these as `[Unnamed pathway (R-HSA-NNNN)]` with a working Reactome link,
which is a graceful fallback but reads as a data hole on the page. We've
seen 1–2 unnamed pathways per disease in the gout / kidney / hepatic
cohorts (transport-heavy gene sets).

**Suggested fix:** the Reactome ingest's `name` resolution likely keys
off a separate file (`ReactomePathways.txt`) that wasn't fully joined for
these older pathway ids — refresh against the upstream pathway-name table.

## Issue #15 — `chembl_molecule` parent/child salt-form linkage exposed only via `childs` on the parent

ChEMBL treats salt forms and anhydrous forms as separate molecule IDs
(e.g. CHEMBL92 = "DOCETAXEL ANHYDROUS" parent, CHEMBL3545252 = "DOCETAXEL"
child; CHEMBL1542 = "AZATHIOPRINE" parent with multiple children). The
parent's `/entry` exposes the child list:
```
entry(CHEMBL92, "chembl_molecule").Attributes.Chembl.childs
  → ['CHEMBL3545252']
entry(CHEMBL1542, "chembl_molecule").Attributes.Chembl.childs
  → ['CHEMBL1200400', 'CHEMBL3785814', 'CHEMBL3785780']
```

But there is **no forward edge from child → parent**. Building a
de-duplication index requires either:
1. Walking every candidate molecule's `/entry` to read `childs`, then
   inverting to a child→parent map (Atlas's §13 disease collector does
   this — adds ~30 entry calls per disease just to dedupe drug names), or
2. A separate `>>chembl_molecule>>chembl_molecule` self-loop map which
   doesn't currently exist.

**Suggested fix:** surface `parent_molecule_chembl_id` (the field is
already in ChEMBL's source data) either as a target schema column on
relevant edges or as a queryable attribute, so a single map call can
return parent-deduplicated lists.

**Atlas impact:** the same drug shows up as 2-3 rows in §10 / §13 drug
tables (TAMOXIFEN + TAMOXIFEN CITRATE + TAMOXIFEN HEMICITRATE;
DOCETAXEL + DOCETAXEL ANHYDROUS) without this. Atlas works around it
today by paying the per-molecule `/entry` calls — fine at our scale
(top-30 drugs per disease) but won't scale to drug pages where the entire
ChEMBL graph is walked.

## Issue #16 — No `list-ids` endpoint for a dataset (corpus enumeration)

There's no API to enumerate every id in a given dataset (every Mondo node,
every HGNC gene, every ChEMBL molecule). `/api/search?s=mondo` returns
relevance-ranked search results for a *query*, not "all ids in dataset
mondo".

**Repro:**
```
search(i="", s="mondo")           → empty (requires a search term)
search(i=" ", s="mondo")          → empty / no useful matches
search(i="disease", s="mondo")    → 50 ranked-by-text-match results,
                                    not a full dataset traversal
```

**Why Atlas needs it.** Scaling beyond a curated reference list (we're at
6 genes + 18 diseases, parked-PROD curated 61 diseases) requires:
- enumerating Mondo to find every disease with material cohort signal
  (current Atlas estimate ~5–10k of Mondo's ~25k classes carry useful
  GWAS / CIViC / ClinVar / GenCC xrefs)
- enumerating HGNC similarly for the full ~43k human-gene corpus.

Today's workaround: drive enumeration off an external Mondo `.obo` /
HGNC TSV dump, parse, then loop. Brittle (version-skew between biobtree
and external dump), defeats the "biobtree is the single source" model.

**Suggested fix:** a dataset list endpoint along the lines of:
```
GET /api/list?s=mondo&limit=1000&page_token=<cursor>
  → {schema: "id|name|xref_count", data: ["MONDO:0001|...|N", ...],
     pagination: {next_token: "..."}}
```
Same schema as `search`, same pagination as `map`. Cheap, would unlock
full-corpus runs across every entity type biobtree indexes.

**Atlas impact:** blocks the "all diseases" / "all genes" scale-out.
Without it, every batch run depends on a hand-maintained list.

## Issue #17 — No bulk xref-count check (need one /entry per id to filter by signal)

Even with #16 resolved, we'd next want to filter `~25k Mondo nodes` down to
"only those with non-empty gwas/civic_evidence/clinvar/gencc xref counts".
Today this requires one `/entry` per node:

**Repro pattern:**
```
for mondo_id in [...]:
    en = entry(mondo_id, "mondo")
    # parse en['xrefs']['data'] for "gwas|N", "civic_evidence|N", ...
```
That's 25k HTTP calls just to skip the 60-80% of Mondo nodes with no
Atlas-relevant signal (deprecated terms, leaf-level rare-rare classes).

**Suggested fix:** add xref counts (or a subset — at least gwas /
gwas_study / clinical_trials / civic_evidence / clinvar / gencc /
intogen) to the search and (proposed) list response schemas. So:
```
list(s="mondo") schema: "id|name|xref_count|gwas|gwas_study|civic_evidence|clinvar|gencc|intogen|clinical_trials"
```
Then a single paginated walk produces both the corpus AND the filter
metadata. We'd skip nodes with all-zero columns.

Alternative: a `/api/xref_counts?ids=ID1,ID2,...&s=DS` batch endpoint
that returns one row per id.

**Atlas impact:** combined with #16, unlocks "rank Mondo nodes by Atlas
signal density, take top-N, run pipeline" — the actual workflow for
all-diseases coverage. Without these two, Atlas either ships a curated
seed-list forever, or pays prohibitive per-node entry-call cost just to
filter.
