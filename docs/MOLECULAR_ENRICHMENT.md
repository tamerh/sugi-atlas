# Molecular enrichment & the Gene/Protein entity contract

Status: **design + layer-A prototype** (2026-06-02). Owner-facing record of the
decision to (a) make the gene page model a *gene and its protein product(s)* as
distinct typed entities on one URL, and (b) enrich the protein/molecular layer
with drug-discovery detail biobtree already exposes.

---

## 1. The problem

The gene page answers the **identity / annotation** layer well (12 sections) but
the **molecular-mechanism** layer is thin: §4 structure is just PDB method +
resolution + AlphaFold pLDDT, and the residue-level features §3 already collects
(EGFR: 306 UniProt features — active site, 4 binding sites, 25 PTMs, 25
disulfides, 37 mutagenesis) are dumped as one flat list.

Two coupled issues:

1. **Semantic mixing.** The page is typed `@type: Gene`, but most of its content
   (structure, residues, interactions, drug binding) is a property of the
   **protein**, not the locus. Every molecular addition below is protein-level,
   so bolting them onto a `Gene` node makes the mixing acute.
2. **Latent molecular data.** biobtree exposes mechanism data the collectors
   touch only at the headline-count level (BRENDA kinetics/inhibitors), sample
   only (BindingDB), or had blocked (GtoPdb, now live).

## 2. Decision: one URL, N typed entities — do **not** split the page

The SEO-vs-semantics tension is a false choice. The semantic web is built for
many typed entities in one document, each with its own `@id`. Keep the single
consolidated URL (good for SEO and for the clinical unity users actually search
for); fix the semantics by making each **protein product a first-class,
addressable entity** rather than splitting it to its own page.

- **URL stays `/atlas/gene/SYMBOL`.** A slug is an *address keyed by the
  canonical identifier*, not a type assertion — semantics live in `@type`.
  Renaming ripples through `@id` / `sameAs` / the `@reverse` mesh / the
  `manifest` (incl. the new `canon` map) / the dist layout, for zero
  machine-semantic gain. `/protein/` is wrong for ncRNA, `/target/` excludes
  non-targets. **Reframe the human-facing label** (nav/breadcrumb →
  "Genes & proteins" or, in a drug-discovery register, "Targets") instead.
- **Protein gets its own `@id`** as a page fragment `…/gene/EGFR/#protein-P00533`
  — promotable later to a real sub-path (`/gene/EGFR/protein/P00533`) with a 301
  if one product ever outgrows the hub, without breaking the graph.

### Entity contract (JSON-LD)

```
Gene   @id …/gene/EGFR/        encodesBioChemEntity → [Protein @id, …]
Protein @id …/gene/EGFR/#protein-P00533
        isEncodedByBioChemEntity → Gene @id          (reciprocal edge)
        + Bioschemas Protein profile: hasSequenceAnnotation (residue map),
          hasRepresentation (PDB/AlphaFold), hasMolecularFunction, sameAs UniProt
        → the molecular additions (§4 below) hang off the Protein, NOT the Gene
```

A crawler/LLM then sees two correctly-typed nodes joined by a real
`encodes`/`encodedBy` edge — not one blurred `Gene`.

**Feasibility note (verified):** the collectors already stamp every molecular
feature with its source protein — `s03` stores each UniProt feature as
`{"uniprot": u, …}` (looping `for u in a.reviewed_uniprots`), and `s04` tags
AlphaFold rows the same way. The union into a flat list happens **only at
render**. So both the residue map and the per-product zoning are **render-only
over data already keyed per-protein** — including the dual-product case.
(Exception: §4 `pdb` is a flat list not yet tagged by accession — attach to the
canonical product for now; per-product PDB tagging is a small collect tweak.)

### Page body zoning (3 zones)

A pure Gene/Protein split orphans the clinical content, which is the unity users
search for — so three zones:

| Zone (`##`, anchored) | Sections |
|---|---|
| **Gene — the locus** | §1 ids · §2 transcripts · §11 expression · §9 TF-regulation · §5 orthologs · functional-genomics · GeneRIFs |
| **Protein product(s)** `#protein-<acc>` | §3 protein ids · §4 structure · **residue map** · §7 pathways · §8 interactions · §10 drugs/binding |
| **Clinical & disease** | §6 variants (tagged by layer: coding→protein, splice/regulatory→gene) · §12 disease associations |

Judgment calls (not load-bearing): orthologs→Gene (locus conservation),
pathways→Protein (Reactome is protein function), variants→Clinical bridge.

**Degradation that must hold:**
- **ncRNA** (`bundle["_noncoding"]`): the **Protein zone elides entirely**; §6/§12
  are already scrubbed (positional-inheritance) → the page collapses to the Gene
  zone + the existing "non-coding — no protein product" line. No empty scaffold.
- **Dual-product** (CDKN2A → P42771/p16 + Q8N726/p14ARF): the Protein zone
  becomes **per-product subsections**, each anchored `#protein-<acc>` to match
  its JSON-LD `@id`. Free — features are already accession-stamped.

Zone = `##` (H2, anchored) demotes sections to `###` (H3): a bigger TOC/visual
change than a reorder, but the honest hierarchy, and it produces the
`#protein-<acc>` anchors that line up with the JSON-LD `@id`.

## 3. The two-tier data rule

Part of the mission is to surface **new angles** — not to hide anything that's
real. So molecular data is *tiered, not gated*:

- **Tier 1 — consolidated.** Leads, presented plainly: UniProt curated features,
  GtoPdb curated affinities, PDB/AlphaFold.
- **Tier 2 — exploratory.** Follows, with a provenance note **that does real
  work** — e.g. "best reported Ki 8.1 nM across 479 BRENDA rows spanning
  organisms/mutants", not a bare "Ki = 8.1" presented as *the* answer, and not
  hidden either.

**The one hard line — wrong ≠ heterogeneous.** Contaminated/incorrect data
(the `mondo→clinical_trials` over-link #28; float32 artifacts #7) gets **fixed
or omitted** — you cannot caption your way out of wrong. Real-but-heterogeneous
data (BRENDA across organisms, BindingDB assay spread) gets **surfaced with a
provenance note** — that's a feature and a differentiator. The Tier-2 note must
carry the actual provenance (organism / assay count / range), or it's noise.

## 4. Enrichment backlog (all probed live on EGFR/DHFR/imatinib)

| # | Addition | Data status | Effort | Value | Tier |
|---|---|---|---|---|---|
| 1 | **Functional-residue map** — restructure §3 ufeatures into a druggable-residue view (active/catalytic, ligand-binding, PTM, disulfide, glycosylation, mutagenesis-validated); overlay ClinVar/AlphaMissense variants + PDB ligands on the same residue axis | ✅ already collected (per-accession) | Low (render) | ⭐⭐⭐ | 1 |
| 3 | **GtoPdb** target class + curated ligand affinities (`>>uniprot>>gtopdb`, `>>…>>gtopdb_interaction`) — pharmacology-grade, hand-curated | ✅ live (just unblocked) | Med | ⭐⭐⭐ | 1 |
| 2 | **BRENDA** kinetics (`>>…>>brenda_kinetics`, per-substrate Km/kcat) + inhibitors (`>>…>>brenda_inhibitor`, Ki/IC50) | ✅ live, uncovered | Med | ⭐⭐ (enzymes) | 2 |
| 4 | **DepMap/GtoPdb tractability verdict** — one-line "17.5% dependent, not common-essential" + druggable class | ✅ wired | Low | ⭐⭐ | 1 |
| 5 | **BindingDB** → potency-ranked top-N (tightest Ki/IC50 + ligand names) instead of arbitrary sample | ✅ wired (sample) | Low | ⭐⭐ | 2 |
| — | **Rhea** reaction equations (`>>uniprot>>rhea`) — substrate→product | ⚠️ biobtree gap: `equation`/ChEBI participants empty in the map projection | blocked | ⭐⭐ | — |

### Per-source provenance caveats (Tier-2 notes must carry these)
- **BRENDA** mixes organisms (human/mouse/bacterial) **and** mutant enzymes —
  prefer human/wild-type; show range + organism, never a lone extremum.
- **BindingDB** has heavy assay heterogeneity — a single "tightest Ki" is often
  one outlier assay; show assay count / range.
- **GtoPdb** is newly unblocked: verify the `>>uniprot>>gtopdb` route against the
  biobtree edges doc and spot-check 2–3 targets for contamination before
  building UI (it's gene-keyed, so structurally unlike the contaminated
  `mondo→` edges — but verify). Update the stale `COLLECTOR_NOTES` "PENDING".
- **Rhea**: file the empty-equation gap in `docs/BIOBTREE_ISSUES.md`.

## 5. Sequencing

The molecular enrichment is the *forcing function* for the protein split, so do
them as one render-change wave and verify once:

1. **Layer A** — Gene + Protein typed nodes (`@id`, `encodes`/`encodedBy`,
   Bioschemas seed). Self-contained; sets the contract the molecular sections
   plug into. *(prototype: this commit)*
2. **Layer B** — body zoning (the 3-zone table) + **#1 residue map** in the
   Protein zone + the nav label reframe.
3. New-data pickups **#3 → #5 → #2**, each in the two-tier style; file Rhea.
   *(shipped 2026-06-02)*
4. One verification pass over the whole wave vs the dense-test baseline.

### Shipped status (2026-06-02)
- **#1** residue map — ✅ layer B.
- **#3 GtoPdb** — ✅ §10: target class + top curated ligand interactions
  (Tier 1). EGFR → catalytic_receptor / ErbB; 106 interactions led by
  panitumumab/cetuximab/osimertinib-class inhibitors with action types.
- **#5 BindingDB** — ✅ §10: promoted sample → potency-ranked top-25 (Tier 2),
  human-filtered, Ki/IC50/Kd/EC50 normalized to nM, `::`-joined ligand names
  cleaned, "heterogeneous assays — not directly comparable" note.
- **#2 BRENDA kinetics** — ✅ §3 (enzymes only): per-substrate Km, ranked by
  measurement depth, Tier-2 note "aggregated across organisms/conditions". A
  max-only value from a 0-floored aggregate renders `—` (would misrepresent a
  µM substrate as mM); the inhibitor side was skipped (redundant with
  GtoPdb/BindingDB + noisy IUPAC names).
- **Rhea** — ⚠️ blocked, filed as `BIOBTREE_ISSUES #29` (equation field empty in
  the map projection; route works).
- **#4 tractability** — partial: GtoPdb druggable class now shown; DepMap
  dependency already surfaced (§3). A single combined verdict line is deferred
  (spans §3+§10).
- Deferred: variant overlay on the residue axis; web-side nav-label reframe.
