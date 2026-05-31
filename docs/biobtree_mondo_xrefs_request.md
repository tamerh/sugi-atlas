# Expose Mondo's full OBO cross-references (and `disease_has_location` UBERON)

## Summary

biobtree exposes some Mondo cross-references (`MESH`, `OMIM`, `Orphanet`,
`EFO`) via mapping edges but drops the rest of what Mondo's OBO file
carries. The missing xref types are valuable for downstream consumers and
all come from the same upstream OBO biobtree already ingests.

The ask: surface the rest of the OBO xref types as edges /
xref-counts under the existing `mondo` dataset, plus the
`disease_has_location` UBERON link for `associatedAnatomy`.

## What's already exposed (good)

For a MONDO id, `bbmap(mondo_id, ">>mondo>>{efo,mesh,mim,orphanet}")`
returns the curated cross-references. Verified working today on
`MONDO:0016419` (hereditary breast carcinoma):

| edge | result |
| --- | --- |
| `>>mondo>>mesh` | `C562840` |
| `>>mondo>>mim` | `114480` |
| `>>mondo>>orphanet` | `227535` |

## What's missing

The Mondo OBO carries many more xref prefixes for the same term. From
`mondo.obo` (release 2026-05-28), the full distribution across non-obsolete
terms:

| prefix | terms covered | currently in biobtree? |
| --- | ---: | :---: |
| UMLS | 21,503 | no |
| MEDGEN | 21,503 | no |
| GARD | 15,947 | no |
| DOID | 12,091 | no |
| Orphanet | 10,492 | yes |
| OMIM | 10,177 | yes |
| SCTID (SNOMED) | 9,153 | no |
| MESH | 8,211 | yes |
| NCIT | 7,435 | no |
| ICD9 | 5,658 | no |
| icd11.foundation | 4,637 | no |
| OMIA | 3,196 | no |
| EFO | 2,400 | yes |
| NANDO | 2,345 | no |
| ICD10CM | 2,141 | no |
| MedDRA | 1,488 | no |
| NORD | 911 | no |
| ICDO | 769 | no |
| OMIMPS | 621 | no |
| HP | 579 | no |
| ICD10WHO | 209 | no |
| Wikipedia | 86 | no |
| GTR | 73 | no |

Concretely: `MONDO:0007254` (umbrella "breast cancer") has DOID/ICD10CM/
ICD11/NCIT/SCTID/UMLS/MEDGEN in the OBO but biobtree exposes none of these,
so Atlas's federated-identifier table renders empty for that term despite
the upstream data being present.

## Bonus ask — UBERON anatomy

Mondo's OBO encodes anatomical location via `intersection_of` clauses:

```
intersection_of: disease_has_location UBERON:0000310 ! breast
```

2,539 MONDO terms (mostly cancers and tissue-specific diseases) carry
this. Exposing it as a `>>mondo>>uberon` edge — or as a `disease_has_location`
field on the mondo entry's Attributes — would let downstream consumers
populate schema.org `associatedAnatomy` for ~half the cancer corpus
without re-parsing the OBO.

## Suggested shape

Two ways biobtree could surface these:

**Option A — extend existing edges.** Add the missing prefixes as
`>>mondo>>doid`, `>>mondo>>sctid`, etc., parallel to the existing
`>>mondo>>mesh`. Adds many dataset names but matches the existing
edge pattern.

**Option B — single `>>mondo>>xref` edge with `prefix` field.** Returns
one row per OBO xref with `{prefix, id}` columns. Compact; one edge
covers all 20+ prefixes; no schema explosion. Atlas would filter
client-side by prefix.

Either works for our use case. Option B is probably easier on biobtree's
schema and indexing side, but you know the cost trade-offs better.

For UBERON anatomy, a dedicated `>>mondo>>uberon` edge (with the
`disease_has_location` relation captured in metadata) would be cleanest.

## Atlas-side commitment

We *won't* parse the OBO for xrefs in Atlas — keeping the OBO usage
strictly limited to MONDO-id enumeration during corpus discovery
(retrieving the full term list once per Mondo release). Atlas's
principle is biobtree-as-single-source-of-truth for entity data;
duplicating xref parsing in Atlas would split that source.

## Verification

A simple check after implementation: pick five terms covering both
common (cancer umbrella) and rare (Mendelian) diseases, compare the
xref set biobtree returns to the OBO term block. They should match.

---

*Filed by the sugi-atlas project (sugi.bio/atlas). Discovered while
investigating empty federated-identifier rows on disease pages — turned
out to be Mondo curation sparsity for some terms, but the OBO carries
strictly more xref types than biobtree exposes for the same term.*
