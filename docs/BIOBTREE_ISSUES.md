# biobtree MCP — Issues & Improvement Requests

**Initial filing:** 2026-05-28
**Last updated:** 2026-06-01 (biobtree refresh resolved #12, #14, #15, #18 + Mondo OBO xrefs; #13 partially resolved — clinical + variant work, guideline still empty; speculative asks #21/#22/#23/#24 removed)

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

## Issue #13 — `pharmgkb_guideline` still empty (clinical + variant resolved)

**Status (2026-06-01):** `pharmgkb_clinical` and `pharmgkb_variant`
populated and wired into Atlas §10 (commit 66ed71b). Only
`pharmgkb_guideline` remains empty.

**Repro:**
```
map(i="HGNC:6407", m=">>hgnc>>pharmgkb_guideline")  → n=0   (KRAS)
map(i="HGNC:2621", m=">>hgnc>>pharmgkb_guideline")  → expected to be non-empty (CYP2C19)
```

CYP2C19 has at least 18 CPIC dosing guidelines on pharmgkb.org. Other
canonical pharmacogenes (CYP2D6, TPMT, DPYD, VKORC1, SLCO1B1) equally
empty.

**Suggested fix:** confirm CPIC/PharmGKB clinical guideline tables are
in the ingest pipeline. Same likely root cause as the original #13
filing — ingest scheduled but data file not yet loaded.

**Atlas impact:** PharmGKB clinical guidelines (CPIC dosing tables) are
the single most-cited PGx-decision source. Without them, Atlas's §10
PharmGKB block can show the clinical-annotation triples but not the
clinical-grade dosing decisions. More pronounced once drug pages land.

## Retracted

| # | Title | Reason |
|---|---|---|
| #16 | No `list-ids` endpoint for a dataset (corpus enumeration) | Retracted 2026-05-31. Corpus enumeration is naturally upstream of biobtree — HGNC ships `hgnc_complete_set.txt` (all human genes), Mondo ships its OBO at obofoundry.org (all disease classes), ChEMBL ships SQLite/flat-file releases. biobtree itself consumes these source files; Atlas can parse the same files for "discover the corpus" without asking biobtree to duplicate that role. Filing was a "biobtree is single source" aesthetic, not a real engineering need. |
| #17 | No bulk xref-count check | Retracted 2026-05-31. The mondo entry already exposes xref counts; calling `/entry` per id is N round-trips but only once per release, easily cached locally. At ~25k Mondo nodes this is ~30 min serial / a few min parallel — not prohibitive. The motivation was that #16 made enumeration cheap and #17 made filtering cheap; once #16 is retracted, #17 loses its main case. A future "batch entry resolve" endpoint might help at scale, but no current bottleneck justifies filing it now. |
