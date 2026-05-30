# biobtree MCP — Issues & Improvement Requests

**Initial filing:** 2026-05-28
**Last updated:** 2026-05-30

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
