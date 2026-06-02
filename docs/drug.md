# How a drug page is built

This walks through how Atlas turns a drug into a published page — the biobtree
chains it queries, the sources behind each section, and the curation rules that
decide what is shown. For the shared machinery (transport, the collect→render
split, the page contract, JSON-LD), see [how-it-works.md](how-it-works.md); this
doc is drug-specific.

The build resolves the drug once to an anchor record, runs twelve section
collectors, and renders them into five canonical, frozen page zones.

```
ChEMBL id (or name)
   │  resolve()            one DrugAnchors record, shared by all sections
   ▼
DrugAnchors ─►  12 section collectors  (biobtree chains → structured bundle)
   │
   ▼
render_all ─►  5 canonical H2 zones  (sections demoted to H3, stable {#anchors})
   │
   ▼
assemble_page ─►  frontmatter · ## Summary · body · ## Related Atlas pages
```

---

## 1. Resolving the drug

Resolution (`drug/anchors.py`) anchors at a **ChEMBL molecule id**. A ChEMBL id
input is a direct `entry` fetch (unambiguous); a name input is searched and
ranked by a three-key rule — *exact name match > named-over-nameless >
xref count* — because a nameless, phaseless screening entry can otherwise outrank
the real drug on raw xref count (the documented *Vemurafenib* case). The corpus
seeds and the batch driver therefore pass the **id, not the name**, so this
decision is pinned once.

The resolver does ~15–25 biobtree calls to assemble the `DrugAnchors` record:

| Field | Source | How |
|---|---|---|
| `chembl_id`, `canonical_name`, `molecule_type`, `max_phase`, `atc_codes`, `alt_names` | ChEMBL | the molecule entry |
| `parent_chembl`, `child_chembls` | ChEMBL | salt/parent linkage (see below) |
| chemistry (CID, InChIKey, SMILES, IUPAC, formula, MW, ChEBI definition + roles, FDA flag) | PubChem + ChEBI | `>>chembl_molecule>>pubchem`, `>>chembl_molecule>>chebi` |
| `targets`, `bioactivity_targets` | GtoPdb / ChEMBL | curated vs bioactivity (see §2) |
| `indications` | ChEMBL + EFO/MeSH→MONDO | `_resolve_indications` |

Two resolution choices shape the page:

- **Chemistry comes from PubChem + ChEBI, not ChEMBL.** The descriptors and the
  drug-class semantics aren't on the ChEMBL molecule entry; they're fetched from
  the open-licensed PubChem and ChEBI records. Only the raw ATC *code* is kept
  (the WHO ATC name table is licence-restricted) — the class label is derived
  from ChEBI roles instead.
- **Salt → parent collapse.** ChEMBL treats salt and anhydrous forms as distinct
  molecule ids (*imatinib* `CHEMBL941` ↔ *imatinib mesylate* `CHEMBL1642`). The
  anchor reads both directions, which is used to render navigation lines, to
  exclude a drug's own salts from its "related molecules", and to register the
  drug in the cross-entity mesh under its own id, its parent, and all children.

**Drug class from ChEBI roles.** A ChEBI role list is arbitrarily ordered and
mixes pharmacology with toxicology and environmental annotations (atorvastatin's
roles include *environmental contaminant, xenobiotic*; caffeine leads with
*mutagen*). Leading the mechanism with `roles[0]` prints nonsense, so a small
denylist (never-a-class roles) + allowlist (mechanism keywords) picks the one
headline class label; the full role list is still shown verbatim in the
identifiers section.

---

## 2. The five zones

The published page is a frozen sequence of five `## H2` zones plus Summary and
Related (the contract is in [PAGE_CONTRACT.md](PAGE_CONTRACT.md)). Twelve
collectors fold into them; three (roles, patents, salt-forms) are merged into the
Identifiers zone rather than rendered as one-line orphan sections.

| Zone `{#anchor}` | Built from | What it answers |
|---|---|---|
| **Summary** `{#summary}` | declarative lead + at-a-glance + JSON-LD | the one-line "what is this drug" |
| **Identifiers** `{#identifiers}` | §1 IDs (+ ChEBI roles, patents, salt lines) | IDs, class, chemistry |
| **Targets** `{#targets}` | §2 targets, §3 bioactivity, §8 target pathways | what it acts on, and the biology around those targets |
| **Indications & clinical** `{#indications}` | §4 indications, §5 trials, §10 CIViC | what it treats, trials, precision-oncology evidence |
| **Pharmacology** `{#pharmacology}` | §9 pharmacogenomics | genotype-guided dosing (CPIC/DPWG) |
| **Related molecules** `{#related-molecules}` | §7 related molecules | same-mechanism competitors |
| **Related Atlas pages** `{#related}` | cross-entity [mesh](mesh.md) | linked genes & diseases |

### Summary
A deterministic lead sentence (no model) — class clause (approved/phase, modality,
ChEBI-role mechanism, ATC), targets clause (top-3 target gene symbols from the
curated set), indications clause, and the top CIViC association — followed by an
at-a-glance digest and the inline JSON-LD.

### Identifiers (§1)
The identity/classification table (de-SHOUTed name, ChEMBL/PubChem/ChEBI/ATC ids,
type + phase), SMILES/IUPAC/formula/MW, the ChEBI definition, the full ChEBI-role
class list, alt-names, salt/parent navigation, and a patent-mention count.

### Targets — the load-bearing decision
This is where the drug page's correctness is decided.

- **Primary targets are GtoPdb-curated mechanism targets only.** The chain
  resolves the GtoPdb ligand (preferring the authoritative
  `>>chembl_molecule>>gtopdb_ligand` ID-join, falling back to name search for
  antibodies, whose GtoPdb ligands carry no ChEMBL id), then
  `>>gtopdb_ligand>>gtopdb_interaction` for the action + affinity, then
  `>>gtopdb>>uniprot>>hgnc` to land on the human gene. Each primary target is
  annotated with its DepMap cancer-dependency signal and links to its Atlas gene
  page.
- **The raw ChEMBL `chembl_target` set is *not* promoted to primary.** It is a
  promiscuous bioactivity cloud — treating it as curated made *Salmeterol*
  "target" TP53/SMN1. It is kept separately as `bioactivity_targets` (count +
  name sample, not gene-resolved), clearly labelled secondary. A drug with no
  GtoPdb ligand therefore has an honestly *empty* primary-target table. The same
  restriction is re-enforced in JSON-LD: a `schema:target` is asserted only for
  GtoPdb-sourced targets.

§3 adds the drug's own ChEMBL **bioactivity** rows (kept at pChembl ≥ 5), and §8
shows the **pathways the target genes sit in** by fanning the gene pathway
collector over the (small) set of target genes.

### Indications & clinical
- **§4 indications** — labelled/clinical disease uses from the ChEMBL molecule
  record, cross-walked EFO/MeSH→MONDO so each links to its Atlas disease page.
  Unmapped rows (e.g. a raw `MP:` mouse-phenotype id) are dropped, not rendered
  as fake diseases.
- **§5 trials** — `>>chembl_molecule>>clinical_trials`, with phases normalised
  (biobtree's `NaN` → "Not specified").
- **§10 CIViC** — the drug-anchored variant × indication × effect triples, via the
  ID-join `>>chembl_molecule>>civic_evidence` (not therapy-name matching), ranked
  by evidence level. Heavy for targeted therapies, absent for non-precision drugs.

This is the mirror of the **disease→drug** curated edge: on a disease page, drug
links come only from *trials + CIViC* — drugs actually tested in or curated for
the disease — never a bioactivity-inferred association. Curated provenance over
assay-cloud inference, on both ends.

### Pharmacology (§9)
Drug-level pharmacogenomics — CPIC/DPWG genotype-guided dosing for *this drug*
(e.g. atorvastatin × SLCO1B1), reached by the ID-join
`>>chembl_molecule>>pubchem>>pharmgkb>>pharmgkb_guideline` (there is no direct
ChEMBL→PharmGKB edge; the PharmGKB chemical node is reached through PubChem).

### Related molecules (§7)
Same-mechanism competitors sharing ≥1 of the drug's **primary** (curated) targets
— fanning over the bioactivity cloud instead would flood it with off-target
overlap. Two sources are merged and deduped by name (with salt/hydrate suffixes
stripped so *imatinib mesylate* collapses to *imatinib*): clinical-stage ChEMBL
candidates and approved/known PubChem drugs, ranked by shared-target count.

---

## 3. Why the data is trustworthy

- **Two-tier targets.** Every drug→gene assertion is either tier-1 *curated*
  (GtoPdb mechanism: gene-resolved, action + affinity, asserted in prose and
  JSON-LD) or tier-2 *bioactivity* (ChEMBL assay cloud: count + name only, never
  gene-resolved, never a `schema:target`). The split is structural — separate
  fields, separate render, separate JSON-LD treatment.
- **ChEBI roles over ATC names**, gated so the headline class is never a
  toxicology or environmental role.
- **Honest empty over wrong.** A drug with no curated target shows an empty
  primary table, not a bioactivity guess promoted to "target".
- **Display hygiene.** Names are de-SHOUTed for display (IMATINIB → Imatinib) but
  never gene symbols; affinities round through `fnum` so float32 artifacts don't
  leak; trial phases normalise `NaN` → "Not specified"; unmapped ontology ids are
  dropped rather than shown as diseases.
- **Resolution ranking** (exact > named > xref) and id-based seeding stop nameless
  screening entries from hijacking the canonical drug.
