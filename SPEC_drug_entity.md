# Drug entity — implementation spec

**Status:** scoped; ready to implement.
**Owner:** open.
**Target layout:** `src/atlas/drug/` — parallel implementation to
`src/atlas/gene/` and `src/atlas/disease/`. No base-class abstraction; same
"concrete, share-on-promotion" philosophy as the existing two entities.

---

## 1. Goal

Build deterministic per-drug pages for Atlas, anchored at ChEMBL molecule
IDs, covering identifiers, targets, bioactivity, indications, trials,
pharmacology, related molecules, target pathways, and pharmacogenomics —
the same deterministic-collector pattern proven on gene + disease.

The parked Hugo workflow at `/data/biobtree-content/biobtree/drug/`
defines 9 LLM-driven sections and a 20-drug curated backlog (Imatinib,
Trastuzumab, Pembrolizumab, Osimertinib, Vemurafenib, Olaparib, ...) —
nothing has been built. A single LLM-shot prototype exists for Alectinib
in pre-9-section format; it's a content reference, not a production
artifact.

---

## 2. Architecture decision

**Parallel-not-generalized**, same as gene + disease. Three reasons recapped:

1. Drug anchors carry shapes the other entities don't (ATC codes, target
   list of UniProt accessions, indication list, parent/child salt-form
   linkage, max development phase).
2. Section content differs — drug §3 is bioactivity, not protein IDs;
   drug §7 is mechanistic-competitor molecules, not pathways.
3. Reuse happens through composition (fanning gene collectors over the
   drug's targets), not inheritance.

Shared with the other entities (already in place):
- `atlas/section.py` — `Section` dataclass metadata contract
- `atlas/render_common.py` — `table()` markdown primitive
- `atlas/civic.py` — CIViC predictive-evidence aggregation
- `atlas/page/uniprot_cc.py` — evidence-code stripping for narratives
- `atlas/page/provenance.py` — UPSTREAM map + per-section walker
- `atlas/page/links.py` — (post-build link decoration; will be written
  during disease scale-out)
- `atlas/pipeline.py` — `assemble_page()` will gain a third
  `entity_type == "drug"` branch + new `run_drug()` entrypoint

---

## 3. DrugAnchors

`src/atlas/drug/anchors.py` — single resolve once, all sections consume:

```python
@dataclass(frozen=True)
class DrugAnchors:
    name: str                   # caller-given (e.g. "Imatinib")
    chembl_id: str              # CHEMBL941
    chembl_entry: dict          # full entry for traceability
    canonical_name: str         # 'IMATINIB' from entry.name
    molecule_type: str          # 'Small molecule' | 'Antibody' | 'Oligonucleotide'
    max_phase: int              # 0..4
    atc_codes: tuple            # e.g. ('L01EA01',) — pharmacology classification
    alt_names: tuple            # brand names + synonyms + chemistry names
    parent_chembl: Optional[str]    # if this is a salt/anhydrous form
    child_chembls: tuple        # if this is the parent form
    # Pre-resolved target list (the drug's pharmacological mode of action):
    targets: Tuple[TargetAnchor, ...]   # see below
    # Pre-resolved indication list (drug's disease label data):
    indications: Tuple[IndicationRecord, ...]
    xref_counts: Dict[str, int]      # from entry.xrefs
    # Chemistry descriptors — NOT on the chembl_molecule entry; sourced from
    # the pubchem + chebi dataset entries (probed: both carry the full set,
    # both open-licensed, no network-in-collect). resolve() fetches them.
    pubchem_cid: Optional[str]
    chebi_id: Optional[str]
    inchi_key: Optional[str]        # pubchem.inchi_key / chebi.inchi_key
    smiles: Optional[str]           # pubchem.smiles / chebi.smiles
    iupac_name: Optional[str]       # pubchem.iupac_name
    molecular_formula: Optional[str]
    molecular_weight: Optional[str]
    chebi_definition: Optional[str] # one-line chemical-class narrative for §1/lead
    chebi_roles: tuple              # decoded ChEBI role names (open-licensed
                                    # drug-class semantics — e.g. "tyrosine
                                    # kinase inhibitor", "antineoplastic agent")
    is_fda_approved: Optional[bool] # from the >>chembl_molecule>>pubchem row

@dataclass(frozen=True)
class TargetAnchor:
    chembl_target_id: str       # CHEMBL2107
    target_type: str            # 'SINGLE PROTEIN' | 'PROTEIN COMPLEX' | ...
    target_name: str
    uniprot: Optional[str]      # P22392 (resolved via chembl_target>>uniprot)
    gene_symbol: Optional[str]  # ALK (via uniprot>>hgnc, with workaround)
    hgnc_id: Optional[str]      # HGNC:427

@dataclass(frozen=True)
class IndicationRecord:
    efo_id: Optional[str]
    mesh_id: Optional[str]
    mondo_id: Optional[str]     # mapped from efo or mesh via biobtree
    name: Optional[str]
    max_phase: int              # phase for THIS indication specifically
    slug: Optional[str]         # for cross-link to Atlas disease page
```

**`resolve(name_or_id)` steps:**
1. If `CHEMBL...` → use directly; else `search(name, source="chembl_molecule")`
   filtered to first reviewed (`highestDevelopmentPhase >= 1`) result.
2. `entry(chembl_id, "chembl_molecule")` → fill `chembl_entry`,
   `canonical_name`, `molecule_type`, `max_phase`, `atc_codes`,
   `alt_names`, `parent_chembl`, `child_chembls`, raw `indications[]`.
3. **Salt-form normalization** — if `parent_chembl` is set, optionally
   resolve from parent instead (avoid building "DOCETAXEL" + "DOCETAXEL
   ANHYDROUS" as two pages). Surface a "this is a salt form of X" note;
   the parent gets the full page.
4. `map_all(chembl_id, ">>chembl_molecule>>chembl_target")` → for each
   target, fetch UniProt + resolve to HGNC via fallback chain
   (`>>uniprot>>hgnc` since direct `chembl_target>>hgnc` returns 0 today
   — file as biobtree #18 if not already on the list).
5. Resolve indications: `entry.indications[]` has `{efo, mesh,
   highestDevelopmentPhase}` (probed: 52 entries for Imatinib, no `mondo` —
   cross-walk each to MONDO via biobtree's `efo`/`mesh`→`mondo` edges);
   compute Atlas slug for cross-link.
6. **Chemistry block** — `map_all(chembl_id, ">>chembl_molecule>>pubchem")`
   → CID (+ `is_fda_approved` flag from the row); `entry(cid, "pubchem")`
   fills smiles / inchi_key / iupac_name / molecular_formula /
   molecular_weight. `map_all(..., ">>chembl_molecule>>chebi")` → ChEBI id;
   `entry(chebi_id, "chebi")` fills `definition` + `roles[]` (resolve each
   role `CHEBI:NNN` to its `name` — one entry call per role, ~2-4 roles).
   Biologics have no small-mol CID → block elides.

**Cost:** ~15-25 biobtree calls per drug at resolve time (added pubchem +
chebi entry + per-role name lookups).

---

## 4. Sections (12 total)

Mirror the parked 9-section design + 3 additions surfaced by the
biobtree probes. Each lives in `src/atlas/drug/sections/s01..s12_*.py`,
each exports a `SECTION = Section(...)` metadata record.

| § | Name | NEW or REUSE | Approach |
|---|---|---|---|
| 1 | drug_ids | NEW | IDs (chembl_id, pubchem CID, chebi, mesh/efo, ATC, type, max_phase, salt-form chain) + **chemistry block from pubchem/chebi entries** (InChIKey, SMILES, IUPAC name, molecular formula, molecular weight) + ChEBI one-line `definition` + filtered alt_names (brand/generic/INN — drop the IUPAC chemistry strings) |
| 2 | targets | NEW | **primary: GtoPdb curated mechanism targets** (`anchors.targets` — target + action [Inhibition/Agonist] + pAffinity; covers antibodies, which ChEMBL bioactivity misses). **secondary: `anchors.bioactivity_targets`** (raw chembl_target set — broader, bioactivity-derived). Link to gene page per target via `links.py`. Resolved via GtoPdb (`gtopdb_ligand→gtopdb_interaction→gtopdb→uniprot→hgnc`): ID-join `>>chembl_molecule>>gtopdb_ligand` first, name-search fallback for biologics; #18(b) interaction guard kept as interim. ChEMBL gene-resolution is the final fallback when a drug has no GtoPdb ligand (e.g. Metformin) |
| 3 | bioactivity | NEW | `>>chembl_molecule>>chembl_activity` rows; sort by pchembl ≥ 5; top 30 by potency; group by target |
| 4 | indications | NEW | render `anchors.indications`; group by max_phase; link to disease page per indication |
| 5 | clinical_trials | NEW | `>>chembl_molecule>>clinical_trials` direct (~200-1000 per drug); top 20 by phase; full phase distribution; status counts |
| 6 | pharmacology | NEW | **primary: ChEBI roles** — decoded role names (open-licensed drug-class semantics: "tyrosine kinase inhibitor", "antineoplastic agent"; resolved from `chebi.roles[]`). **secondary: raw ATC code(s) + whocc.no link** — NOT decoded to a hierarchy (biobtree exposes only the code `["L01EA01"]`; the WHO ATC name table is licensing-restricted, so we link out rather than reproduce it). RoA where exposed |
| 7 | related_molecules | REUSE | fan over the **GtoPdb primary targets** (curated, ~1-9 per drug — NOT all 85 bioactivity targets, which would flood the competitor list with off-target overlap): for each primary target uniprot, `>>uniprot>>chembl_target>>chembl_molecule[highestDevelopmentPhase>=2]`; aggregate top 30 competitors sharing ≥1 primary target; link to other drug pages |
| 8 | target_pathways | REUSE | for each target uniprot's gene, fan gene §7 `s07_pathways.collect(gene_anchors)`; aggregate Reactome + top GO terms across the target set |
| 9 | pharmacogenomics | NEW (blocked) | `>>chembl_molecule>>pharmgkb_drug` — currently empty (likely BIOBTREE_ISSUES #13 territory); fallback: fan target genes' `>>hgnc>>pharmgkb_gene` |
| **10** | **clinical_evidence (CIViC)** | NEW | drug-anchored view of the same `civic_evidence` table. **Use the ID-join route `>>chembl_molecule>>civic_evidence`** (probed working — same 7-col schema, properly joined; NOT therapy-name string-matching which would miss brand/salt aliases). Feed straight into `atlas.civic.aggregate_predictive`. Render `(variant, indication, effect, level)` — therapy column kept only to show combinations (e.g. "Imatinib + Dasatinib"). Completes the gene/disease/drug CIViC symmetry. |
| **11** | **patent_literature** | NEW | sum patent counts across linked `patent_compound` records (already wired on gene §10 per-molecule; surface drug-anchored here) |
| **12** | **salt_forms_and_parent** | NEW | render `parent_chembl` + `child_chembls` as inline navigation; trivial render-only |

§9 and §10 will surface heavy content for approved targeted-therapy drugs
(Imatinib has ~263 CIViC evidence rows; Sotorasib likely similar). For
non-precision-medicine drugs (Metformin, Atorvastatin) §10 elides naturally.

**Render-only derived sections** (later, if useful — analogous to disease
§15-§17):
- §13 mechanism_summary: one-paragraph "small-molecule ALK inhibitor;
  primary target Q9UM73; OFF-target on RET + ROS1" — could be inferred
  from §2 + §3 deterministically without LLM.

---

## 5. Cohort fan-out helper

`src/atlas/drug/cohort.py` — parallel to disease's `cohort.fan()`:

```python
def fan(gene_collect_fn, target_anchors: tuple[GeneAnchors]) -> list[dict]:
    """Run a gene section's collect_fn over each target's gene anchors.
    Used by §7 (related_molecules — actually shares fan pattern with §8)
    and §8 (target_pathways)."""
```

But for drug, the cohort is small (~3-100 targets per drug; mostly <30).
Fan-out cost is modest.

For §7 specifically the fan-out target is `chembl_target` IDs, not
`hgnc_id` — slightly different shape. Either generalize `cohort.fan()`
to take any anchor list + collect function, or keep two parallel helpers.
Drug side likely wants its own thin `over_targets()` helper to keep
type intent clear.

---

## 6. Render

`src/atlas/drug/render.py` — 12 `r_*` functions + a `RENDER` dict, mirror
`disease/render.py`. Each section gets a `## <title>` header + tables
following the shared `table()` shape.

`atlas/page/drug_jsonld.py` — schema.org `Drug` JSON-LD sidecar:
```json
{
  "@context": "https://schema.org",
  "@type": "Drug",
  "@id": "https://sugi.bio/atlas/drug/imatinib/",
  "name": "IMATINIB",
  "identifier": "CHEMBL941",
  "url": "https://sugi.bio/atlas/drug/imatinib/",
  "alternateName": ["Glamox", "NSC-743414", ...],
  "sameAs": [
    "https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL941/",
    "https://pubchem.ncbi.nlm.nih.gov/compound/<CID>",
    "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:45783",
    "https://www.whocc.no/atc_ddd_index/?code=L01EA01"
  ],
  "drugClass": "Antineoplastic agent",            # from ATC hierarchy
  "mechanismOfAction": "BCR-ABL tyrosine kinase inhibitor",   # from §2 primary target name
  "target": [{"@type": "Gene", "name": "ABL1", "url": "..."}, ...],
  "treats": [{"@type": "MedicalCondition", "name": "...", "url": "..."}, ...]
}
```

`atlas/page/drug_declarative.py` — one-sentence lead:
```
**Imatinib** (CHEMBL941) is a Small molecule BCR-ABL tyrosine kinase
inhibitor (ATC L01EA01); approved for 12 indications including chronic
myeloid leukemia and gastrointestinal stromal tumors; targets 84
proteins primarily ABL1, KIT, and PDGFRA.
```

`atlas/page/drug_provenance.py` — schema.org Dataset trail per section.

---

## 7. Pipeline integration

`atlas/pipeline.py`:
- `run_drug(name, dist_dir, do_summary=True, ...)` mirroring `run_disease`
- `assemble_page` gains `entity_type == "drug"` branch — picks
  `drug_declarative` + `drug_jsonld`
- `--entity drug` on the CLI; `python -m atlas.pipeline drug "Imatinib" --dist ...`

`atlas/validation/body_gate.py`:
- Already entity-agnostic via `snap_dir_for(dist_root, "drug")`. No changes
  needed beyond passing `entity="drug"` from the drug task script.

Snapshots land at `<dist>/snapshots/drug/<slug>.json`.

---

## 8. Discovery + corpus

`src/atlas/drug/corpus.py` — mirror `atlas/disease/corpus.py`:

**Source**: ChEMBL itself ships an SQLite release at
`https://chembl.gitbook.io/chembl-interface-documentation/downloads`.
Parse it for `(chembl_id, name, max_phase, indication_count,
chembl_activity_count, atc_codes)` — fields needed for filtering.

Alternative: walk biobtree's chembl_molecule index by phase filter
(`>>chembl_target>>chembl_molecule[highestDevelopmentPhase==4]`) per
target, dedupe. Less authoritative but avoids the SQLite dump.

**Tiers (proposed):**
| Tier | Filter | Estimated count |
|---|---|---|
| D1 | max_phase == 4 (approved) | ~2,500 |
| D2 | max_phase ≥ 3 | +1,500 |
| D3 | max_phase ≥ 2 | +3,000 |
| D4 | max_phase ≥ 1 (any clinical) | +5,000 |
| D5 | high-activity preclinical (≥100 chembl_activity rows w/ pchembl ≥ 5) | +10,000 |

Approved drugs (D1) are the right MVP — clinical relevance is universally
high. ~2,500 pages × ~30s each = ~21h sequential / 2-3h parallel.

**CLI:**
```bash
python -m atlas.drug.corpus build                  # parse ChEMBL release
python -m atlas.drug.corpus filter --tier D1       # emit {name, slug} list
python -m atlas.drug.corpus run --tier D1 --limit 20  # dev-batch smoke
```

---

## 9. Enju workflow

`src/atlas/drug/enju.yaml` — same 4-task DAG as disease:

```yaml
name: "Atlas drug page"
params:
  - name: drugs
    type: list<record>
    required: true
    key: slug
    fields:
      name: string
      slug: string
for_each:
  drug: "{{drugs}}"
tasks:
  - id: collect_render
    script: src/atlas/drug/tasks/collect_render.py
    writes: [build/{{drug.slug}}/bundle.json, build/{{drug.slug}}/body.md, build/{{drug.slug}}/anchors_meta.json]
  - id: body_gate
    depends_on: [collect_render]
    ...
  - id: summary
    depends_on: [collect_render]
    ...
  - id: publish
    depends_on: [body_gate, summary]
    ...
```

Mirror disease's `tasks/` scripts exactly: `collect_render.py`,
`body_gate.py`, `summary.py`, `publish.py`. Iteration variable is
`drug.{name, slug}`; `name` resolves via `DrugAnchors`, `slug` is the
caller-pinned URL key.

---

## 10. Validation set (smoke + reference)

Reuse the parked `entities.yaml` 20-drug list. Recommended smoke order
(mix of targeted-therapy + biologic + small-molecule + metabolic):

| # | Drug | ChEMBL | Why |
|---|---|---|---|
| 1 | Imatinib | CHEMBL941 | classic kinase inhibitor; rich data on every section |
| 2 | Trastuzumab | CHEMBL1201585 | antibody (test non-small-mol shape: no SMILES, no patent_compound) |
| 3 | Sotorasib | (resolve by name) | spec-flagship; tests CIViC §10 + KRAS G12C cross-link. NOTE: don't hardcode the id — `CHEMBL4439653` is empty in this biobtree; resolve-by-name finds the populated `CHEMBL4535757`. General rule: the validation set resolves by name, never pinned ChEMBL ids. |
| 4 | Osimertinib | CHEMBL3353410 | similar — EGFR T790M cross-links |
| 5 | Pembrolizumab | CHEMBL3137344 | checkpoint inhibitor antibody; tests many indications |
| 6 | Olaparib | CHEMBL521686 | PARP inhibitor; BRCAness narrative |
| 7 | Vemurafenib | CHEMBL1229517 | BRAF V600E |
| 8 | Crizotinib | CHEMBL601719 | ALK + ROS1 + MET (multi-target) |
| 9 | Alectinib | CHEMBL1738797 | matches existing prototype for content comparison |
| 10 | Venetoclax | CHEMBL3137309 | BCL2 inhibitor |
| 11 | Ibrutinib | CHEMBL1873475 | BTK |
| 12 | Palbociclib | CHEMBL2105716 | CDK4/6 |
| 13 | Bevacizumab | CHEMBL1201583 | VEGF antibody |
| 14 | Nivolumab | CHEMBL2103882 | checkpoint antibody |
| 15 | Lenalidomide | CHEMBL848 | molecular glue |
| 16 | Erlotinib | CHEMBL553 | classic EGFR-TKI |
| 17 | Gefitinib | CHEMBL939 | same; salt-form testing |
| 18 | Sorafenib | CHEMBL1336 | multi-kinase |
| 19 | Metformin | CHEMBL1431 | non-cancer; tests sparse-CIViC behavior |
| 20 | Atorvastatin | CHEMBL1487 | non-cancer; tests pharmacology+PGx (statin myopathy) |

After T1 corpus build, the canonical-drug list expands to ~2,500
approved drugs.

---

## 11. Cross-link layer (drug-side responsibility)

When `atlas/page/links.py` exists (built during disease scale-out's
post-build link pass), the drug page picks it up for three classes of
internal links:

| On a drug page | Resolves to |
|---|---|
| §2 targets table — gene_symbol column | `/atlas/gene/<symbol>/` if symbol in HGNC corpus |
| §4 indications table — mondo_id / disease name | `/atlas/disease/<slug>/` if mondo_id in disease corpus |
| §7 related_molecules table — competitor name | `/atlas/drug/<slug>/` if in drug corpus |
| §10 CIViC clinical evidence — gene + disease | both link types apply |

Symmetrically, gene §10 and disease §13 already render drug names (from
ChEMBL/CIViC) — when the drug corpus exists those should become clickable
to `/atlas/drug/<slug>/`. The drug corpus completes the three-entity
reciprocal-link mesh.

---

## 12. Open biobtree gaps (file if confirmed during implementation)

| # | Issue | Atlas impact |
|---|---|---|
| **#18** | GtoPdb drug→target: (b) `gtopdb_interaction` id-substring contamination — ✅ FIXED upstream (pending gtopdb re-update); (a) antibodies have no chembl→gtopdb key — documented biologics gap, not a bug | §2/§7 — `anchors.py` resolves ligand by ID-join (`>>chembl_molecule>>gtopdb_ligand`) first, name-search fallback for biologics; interaction guard kept as interim for (b), no-op after the re-update. |
| #13 | `chembl_molecule>>pharmgkb_drug` empty | §9 PGx — work around via gene fanout; same root cause as existing BIOBTREE_ISSUES #13 |

**Not gaps (confirmed against biobtree's edges doc, `/data/biobtree/docs`):**
- `chembl_target>>hgnc` direct returns 0 **by design** — chembl_target's xref is
  to UniProt; go `chembl_target>>uniprot>>hgnc` (what `anchors.py` does).
- `chembl_molecule>>mondo` absent — indications cross-walk via efo/mesh
  (confirmed working; handled in `_resolve_indications`).
- `opentargets` is a derived/identifier namespace (in `/api/meta` but not one of
  the 76 edge-bearing datasets) — not traversable, and not meant to be.
- indications shape confirmed `{highestDevelopmentPhase, efo, mesh}` (no `mondo`).

**Resolved during spec iteration (no biobtree request needed):**
- ~~`chembl_molecule.smiles` / `inchiKey`~~ — NOT on chembl_molecule, but
  fully available on the **pubchem** + **chebi** dataset entries (smiles,
  inchi_key, iupac_name, formula, MW). Route via `>>chembl_molecule>>pubchem`
  / `>>chebi`. §1 is richer, not poorer.
- ~~ATC hierarchy decode~~ — biobtree gives only the raw code; WHO's name
  table is license-restricted. **Use ChEBI `roles` instead** (open-licensed,
  decode to "tyrosine kinase inhibitor" etc.); ATC stays a raw code + link.

Possible remaining request (probe before filing):
- `chembl_molecule >> route_of_administration` or similar — for §6 RoA.
  File only after probing confirms; don't pre-file aspirational requests.

## Engineering notes (from spec-iteration probes)

- **Corpus-level salt-form dedup.** The D1 (max_phase==4) ChEMBL list will
  contain parent *and* child salt forms as separate rows. Fold `childs` into
  the parent at corpus-build so we don't ship "DOCETAXEL" + "DOCETAXEL
  ANHYDROUS" — reuse disease §13's existing parent/child fold logic.
- **§7 target fan cap.** Imatinib has 85 chembl_targets; fanning
  `chembl_molecule` over all 85 (many off-target/promiscuous) is the heaviest
  + noisiest section. Cap the fan to **primary targets** (by potency/activity
  count), not the full target list, or the competitor list fills with
  off-target overlap.
- **Biologics render conditionally on `molecule_type`.** Antibodies
  (Trastuzumab, Pembrolizumab, Bevacizumab, Nivolumab — half the smoke set)
  have no small-mol CID → §1 chemistry block, §3 bioactivity, §11 patents
  elide. Gate those sections on `molecule_type == "Small molecule"` so an
  antibody page doesn't render empty small-molecule scaffolding.

---

## 13. Sequencing (not a timeline)

Recommended order for the implementing agent:

1. **Probe + DrugAnchors** — implement `anchors.py`, smoke-test resolve
   on the 5 most diverse drugs in the 20-list (Imatinib, Trastuzumab,
   Sotorasib, Olaparib, Metformin). Validate the anchor shape carries
   what's needed without re-fetching.
2. **5 NEW collectors first** — §1, §3, §5, §6, §11. These are pure
   chembl-anchored, no fan-out, simpler. Self-contained smoke per drug.
3. **2 REUSE collectors** — §7, §8. Wire `cohort.fan` analog over the
   target list; verify Imatinib's pathway aggregation looks coherent
   (BCR-ABL → ABL1 → ABL1-mediated pathways).
4. **§10 CIViC integration** — reuses `atlas.civic.aggregate_predictive`;
   filter civic rows where therapy matches drug name or alt_name.
5. **§12 salt-form linkage** — trivial render-only; parent/childs from
   anchors.
6. **§9 PGx** — last; partially blocked. Wire the fallback path
   (target-gene-fanout) first; the direct chembl_molecule>>pharmgkb_drug
   route gets enabled when biobtree #13 lands.
7. **§2 + §4 + §7 cross-link decoration** — coordinate with disease
   scale-out's `atlas/page/links.py` — once that exists, plug in.
8. **Render + pipeline + workflow + task scripts** — mirror disease at
   the file-shape level.
9. **Smoke** — top-5 from the validation table end-to-end. Body_gate +
   sidecars + dist write.
10. **Corpus + tiered scale-out** — `corpus.py` builds the D1 corpus from
    ChEMBL; build all approved drugs (~2,500) in one tier.

**Out of scope for this pass:**
- LLM summary pass (defer until full corpus exists, like disease)
- Drug → drug similarity beyond shared-target (ChEMBL has chemical
  similarity; out of scope)
- Drug-drug interactions (would need a separate biobtree dataset
  ingestion)
- DrugBank ID surfacing (biobtree doesn't currently ingest DrugBank;
  ChEMBL parent IDs cover most use cases)

---

## 14. What's already done

- Parked Hugo workflow has the 9-section design + 20-drug entities.yaml +
  prompt drafts → use as content reference, not as the new implementation
- The Alectinib prototype page exists in pre-section format → useful
  as a "what does it look like in production prose" sanity check
- The CIViC aggregation logic (`atlas/civic.py`) is reusable as-is; only
  the filtering predicate (therapy-name match instead of gene/disease
  anchor) differs
- biobtree's chembl_molecule + chembl_target + chembl_activity + 
  clinical_trials + pubchem + chebi + mesh + efo edges all work today
  (probed live during this spec)

---

## 15. Reference probes (verified during spec writing)

For Imatinib (CHEMBL941), the canonical "richly-covered targeted therapy"
test case:

```
entry attrs:        type=Small molecule, max_phase=4, atc=[L01EA01],
                    indications=12 entries, childs=[CHEMBL1642, CHEMBL2386595],
                    altNames=156 entries (brand names + chemistry strings)
xrefs:              chembl_target|85, chembl_activity|424, clinical_trials|223,
                    mesh|45, efo|52, civic_evidence|263, civic_assertion|2,
                    patent_compound|3, gtopdb_ligand|1, chebi|1, pubchem|1
chembl_molecule>>chembl_target            n=85
chembl_molecule>>chembl_target>>uniprot   n=84
chembl_molecule>>chembl_target>>hgnc      n=0      ← workaround: via uniprot>>hgnc
chembl_molecule>>chembl_activity          n=300+
chembl_molecule>>clinical_trials          n=223
chembl_molecule>>patent_compound          n=3
chembl_molecule>>pubchem                  n=1
chembl_molecule>>chebi                    n=1
chembl_molecule>>mesh                     n=45
chembl_molecule>>efo                      n=52
chembl_molecule>>pharmgkb_drug            n=0      ← blocked (likely #13)
chembl_molecule>>mondo                    n=0      ← cross-walk via mesh/efo
```

This is the proven data layer the implementing agent should build against.

---

**END OF SPEC.** Next agent: clone the disease implementation's shape file
by file; substitute drug-specific data; ship the smoke 5 end-to-end before
committing to the full approved-drug build.
