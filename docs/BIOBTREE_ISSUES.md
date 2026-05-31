# biobtree MCP — Issues & Improvement Requests

**Initial filing:** 2026-05-28
**Last updated:** 2026-05-31 (biobtree refresh resolved #9 + #10 fully + #15 partially; #12 mostly resolved; corpus-enumeration asks #16/#17 retracted on reflection; #18 filed during drug-entity build — GtoPdb drug→target routing)

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

## Issue #9 — RESOLVED 2026-05-31 (see top-of-doc Resolved table)

## Issue #10 — RESOLVED 2026-05-31 (see top-of-doc Resolved table)
## Issue #12 — pubchem_activity broadly restored; KRAS still empty (under upstream investigation)

**Status (2026-05-31 refresh):** the original "KRAS and similar targets"
filing is now narrowed. Sampled 12 reference / canonical drug-targetable
proteins:

```
✓ TP53     P04637   n=882
✓ EGFR     P00533   n=2000 (cap)
✓ AKT1     P31749   n=2000
✓ NR3C1    P04150   n=2000
✓ AR       P10275   n=2000
✓ ESR1     P03372   n=2000
✓ CAMK2A   Q9UQM7   n=1168
✓ BRCA1    P38398   n=28
✓ NFKBIA   P25963   n=130
— KRAS     P01116   n=0    (under upstream investigation)
— CDKN2A   P42771   n=0    (likely target-specific; CDK inhibitor)
— TTN      Q8WZ42   n=0    (likely target-specific; extreme size)
```

So the broad-index issue is fixed — pubchem_activity now returns
populated activity sets for the canonical drug-target proteome. **KRAS
specifically remains empty pending upstream investigation.** TTN and
CDKN2A may have target-specific explanations (no surprise that titin
isn't in a drug-target bioactivity index; CDKN2A is a CDK inhibitor with
limited small-molecule chemistry).

**Atlas mitigation (kept in place):** §10 still wires `chembl_activity`
alongside `pubchem_activity` (commit `189f98a`). Both blocks render where
data exists; for KRAS the PubChem block is empty but ChEMBL surfaces
5,239 activities (4,825 potent at pChembl≥5). When the KRAS-specific
investigation completes the PubChem block will backfill automatically.

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

## Issue #15 — `chembl_molecule` parent/child salt-form linkage — partial 2026-05-31

**Partial update 2026-05-31:** child entries now expose a `parent` field
(e.g. `entry('CHEMBL3545252', 'chembl_molecule').Attributes.Chembl.molecule.parent`
returns `CHEMBL92`). Atlas can dedupe with one entry call per *child*
instead of walking every candidate's `childs`. Remaining ask: a forward
edge / map column would still beat the per-entry pattern at drug-page scale.

ChEMBL treats salt forms and anhydrous forms as separate molecule IDs
(e.g. CHEMBL92 = "DOCETAXEL ANHYDROUS" parent, CHEMBL3545252 = "DOCETAXEL"
child; CHEMBL1542 = "AZATHIOPRINE" parent with multiple children). Today's
options:
1. ✓ child → parent: 1 entry call per child via the new `parent` field.
2. parent → children: 1 entry call per parent via `childs` list (unchanged).
3. ✗ no forward map edge (e.g. `>>chembl_molecule>>parent`) — would
   collapse the cost to a single map call.

**Suggested next step:** add `parent_chembl_id` to the target schema of
ChEMBL-molecule-emitting edges (e.g. `>>chembl_target>>chembl_molecule`,
`>>mondo>>clinical_trials>>chembl_molecule`) so the dedupe key is in the
map response itself.

**Atlas impact:** §13's per-child `entry()` dedupe path (atlas/disease/
sections/s13_clinical_trials.py) still pays ~30 calls per disease.
Acceptable at disease scale; will be a hot path on drug pages.

## Issue #18 — GtoPdb drug→target: interaction id substring contamination (FIXED upstream, pending re-update)

GtoPdb (Guide to Pharmacology) is the cleanest source of **curated mechanism
targets** for a drug. Two things were investigated:

**(b) `>>gtopdb_ligand>>gtopdb_interaction` leaked other ligands' interactions
by id-substring — ✅ FIXED upstream (2026-05-31), pending a gtopdb re-update to
take effect (~next release).** Querying ligand `7519` (olaparib) returned
interactions belonging to ligand `5662` (**AT-7519**, a CDK inhibitor), so
olaparib (a PARP inhibitor) appeared to target CDK1–9:
```
map(7519, ">>gtopdb_ligand>>gtopdb_interaction")  →
  1961_5662 | cyclin dependent kinase 1 | AT-7519 | ...   ← WRONG (ligand 5662)
  2771_7519 | poly(ADP-ribose) polymerase 1 | olaparib    ← correct (ligand 7519)
```
**Root cause (dev):** ligand 5662's "AT-7519" synonym was tokenized into the
bare number `7519` and indexed as a gtopdb_ligand keyword, colliding with the
numeric ligand-ID namespace. Fix: `extractSignificantWords` now skips all-digit
tokens (+ a guard on bare-numeric full synonyms). Builds clean; needs a gtopdb
re-update, after which `7519` resolves only to olaparib.

**(a) Antibody coverage — NOT a bug; a documented biologics data gap.** Small
molecules traverse the forward edge fine: `>>chembl_molecule>>gtopdb_ligand`
gives Imatinib→5687, Olaparib→7519. Trastuzumab returns 0 because gtopdb's
**antibody** ligands carry no ChEMBL id — there's no shared key to join on, and
forcing it would need fragile name matching. (Earlier framing of "(a) forward
edge unwired" was wrong — I'd only tested the antibody; the edge works for
small molecules.)

**Atlas handling (in place):** `atlas/drug/anchors.py` resolves the ligand by
the **ID-join first** (`>>chembl_molecule>>gtopdb_ligand`), falling back to
name search only when there's no xref (antibodies + a few small molecules like
Sotorasib). It also keeps an interaction-row guard
(`id.split("_")[-1] == ligand_id`) as the interim fix for (b) — a no-op once
the re-update ships. Trastuzumab→ERBB2, Imatinib→ABL1/DDR1/DDR2,
Olaparib→PARP1/2/3, Sotorasib→KRAS all resolve correctly (action + pAffinity).

> **Not filed (checked against biobtree's edges doc, `/data/biobtree/docs`):**
> - `opentargets` returning n=0 from every route is **not a bug** — it's a
>   *derived/identifier* namespace (appears in `/api/meta` but isn't one of the
>   76 edge-bearing datasets). Routing is defined only over the actual
>   datasets; refer to the edges doc, not `/meta`, for what's traversable.
> - `chembl_target>>hgnc` returning 0 is **by design** — `chembl_target`'s
>   cross-reference is to UniProt (`components.acc`); the gene is reached via
>   `chembl_target>>uniprot>>hgnc`. That's the documented pattern, not a gap.

## Retracted

| # | Title | Reason |
|---|---|---|
| #16 | No `list-ids` endpoint for a dataset (corpus enumeration) | Retracted 2026-05-31. Corpus enumeration is naturally upstream of biobtree — HGNC ships `hgnc_complete_set.txt` (all human genes), Mondo ships its OBO at obofoundry.org (all disease classes), ChEMBL ships SQLite/flat-file releases. biobtree itself consumes these source files; Atlas can parse the same files for "discover the corpus" without asking biobtree to duplicate that role. Filing was a "biobtree is single source" aesthetic, not a real engineering need. |
| #17 | No bulk xref-count check | Retracted 2026-05-31. The mondo entry already exposes xref counts; calling `/entry` per id is N round-trips but only once per release, easily cached locally. At ~25k Mondo nodes this is ~30 min serial / a few min parallel — not prohibitive. The motivation was that #16 made enumeration cheap and #17 made filtering cheap; once #16 is retracted, #17 loses its main case. A future "batch entry resolve" endpoint might help at scale, but no current bottleneck justifies filing it now. |

## Issue #21 — Batch-map / batch-entry endpoint to amortize HTTP overhead at scale

Filed 2026-05-31 during the all-diseases scale-out (Atlas's Phase 2 of
the disease corpus build, ~19.7k pages).

For each disease page, Atlas calls biobtree thousands of times — most of
the cost is HTTP round-trips, not biobtree compute. Each cohort gene's
§7 (pathways), §8 (interactions, now capped), §10 (drugs), §11
(bioactivity) involves N separate `/map` calls with different chains.
Across the full corpus this is hundreds of millions of HTTP calls.

biobtree's per-call response is fast (~6ms). The waste is HTTP-layer:
TCP setup, header serialization, JSON parsing per response.

**Proposed shape:**

```
POST /api/map_batch
{
  "ids":   ["HGNC:11998", "HGNC:1100", "HGNC:6407", ...],   // ≤500 ids per request
  "chains": [">>hgnc>>reactome", ">>hgnc>>uniprot>>interpro"]  // ≤10 chains
}
→ {
  "results": {
    "HGNC:11998": {
      ">>hgnc>>reactome":              {schema: "...", targets: [...]},
      ">>hgnc>>uniprot>>interpro":     {schema: "...", targets: [...]}
    },
    "HGNC:1100":   {...},
    ...
  },
  "pagination": {has_next: ..., next_token: ...}
}
```

Similar shape for `/api/entry_batch` (one POST with N ids + one source).

**Atlas impact:** the disease pipeline fans gene collectors over a 50-gene
cohort × ~9 sections × ~5 chains/section = ~2,250 HTTP calls per disease.
A batched call (50 ids × 1 chain per request) cuts that 50× — same total
data, ~45 calls instead of 2,250. Estimated per-disease drop from ~48s
to ~5-10s. Corpus build: ~12h parallel → ~2h.

**Why this is the right shape (vs alternatives):**

- *list-ids endpoint* — was filed as #16, retracted: corpus enumeration
  belongs upstream of biobtree.
- *bulk xref-count check* — was filed as #17, retracted: one-time per
  release is fine.
- *per-chain HTTP caching headers (ETag/Cache-Control)* — would help repeat
  callers but doesn't address the first-call HTTP volume.
- *batch-map/entry* — directly attacks the HTTP-call-count problem the
  hot path actually has.

**Not blocking — Atlas can ship at the current 48s/page rate** (12h corpus
build with 8-way parallel). This is a "ship the request, swap when it
lands" pattern.
