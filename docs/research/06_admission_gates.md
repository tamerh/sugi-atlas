<!-- Corpus admission-gate design — from the 1,010-page dist quality audit (agents E/F/G)
     + over-exclusion analysis. Goal: drop junk at the corpus EDGE without silently
     dropping real entities. Guiding asymmetry: UNDER-gate > OVER-gate (a little junk is
     recoverable; a missing real entity is invisible). Every gate emits a KEPT + a
     DROPPED(+reason) list we eyeball before applying; keep the ungated `dense-test`
     branch to diff against. -->

# Corpus admission gates — design + over-exclusion analysis

The 1,010-page audit found the plumbing is solid but the **corpus edge is ungated**:
junk ontology nodes, reagent "drugs", salt/subtype duplicates. These gates fix that
*at seed-build time*. The hard requirement (per the "don't miss real stuff" review):
**be surgical, prove the drop-list is pure junk, bias toward keeping.**

Meta-rules for every gate:
1. Emit `seeds/<x>.kept.txt` + `seeds/<x>.dropped.txt` (id, name, reason). Review the
   dropped list before a full run.
2. Prefer a **principled signal** (curator annotation, real pharmacological value) over
   name heuristics.
3. When in doubt, KEEP (let a little junk through; tighten iteratively).
4. The `dense-test` git branch (un-gated 1,010) stays as the before/after diff baseline.

---

## Gate 1 — Disease: drop ontology QUALIFIERS, not groupings  `[disease seed]`

**Trap (measured):** the obvious gate — Mondo's `subset: disease_grouping` (1,368 terms) —
is **far too blunt**. 20 of the top-300 seed are flagged grouping, and they include
**marquee diseases we must keep**: cardiomyopathy, acute myeloid leukemia, renal cell
carcinoma, hypertrophic/dilated cardiomyopathy, RASopathy, Charcot-Marie-Tooth, Leigh
syndrome. A disease having subtypes ≠ junk. **Do NOT gate on `disease_grouping`.**

**Surgical signal:** the junk node `inherited` (MONDO:0021152) is not a disease grouping
at all — its ancestry is `hereditary vs non-hereditary etiology` → **`disease
characteristic` (MONDO:0021125)** → PATO. That subtree is **40 nodes** of pure
qualifiers (inherited / acquired / sporadic / congenital / X-linked …). None are
diseases.

**Gate:** exclude the descendants of `disease characteristic` (≈40 ids). Over-exclusion
risk ≈ 0 (no real disease lives under a PATO quality). Implement by enhancing
`atlas.disease.corpus.parse_obo` to capture is_a parents → compute the characteristic
subtree → exclude. Optionally a tiny hand-denylist for any other obvious non-disease.

**Do NOT** drop on: `disease_grouping`, missing-OMIM/Orphanet (common/GWAS diseases lack
these), or low signal (already filtered to signal>0).

---

## Gate 2 — Drug: require a real PHARMACOLOGICAL value  `[drug seed]`

**Trap (measured):** ChEMBL `max_phase==4` is broad — water, phosphoric acid,
ammonia-solution, evans-blue, selenomethionine all pass it (water even has 24
"indications" + a route + ChEBI role "greenhouse gas"). `molecule_type` is useless
(water is "Small molecule"). The raw JSONL has **no ATC / mechanism / parent** fields.

**Signal — "did ChEMBL/biobtree populate a real pharmacological value?"** Keep a molecule
iff ANY of (generous UNION — a real therapeutic almost always hits one; reagents/dyes hit
none):
- has an **ATC code** (therapeutic classification), OR
- has a **mechanism of action** (drug_mechanism / chembl mechanism.action), OR
- has a **curated target** (GtoPdb mechanism target — the §2 primary-target path), OR
- has a **CIViC or PharmGKB** annotation, OR
- has a **pharmacological ChEBI role** (exclude when the *only* roles are
  metabolite / solvent / contaminant / food-component / greenhouse-gas etc.).

These are resolved by `resolve_drug()` (the anchor already pulls ATC, GtoPdb targets,
parent). So the gate = a one-time **seed-refinement pass**: resolve each of the 4,225
approved → keep those passing the union. ~4k calls, one-time (~10 min); emit the
dropped list to confirm it's pure reagents/dyes/supplements. Start generous; tighten
only if real drugs survive in `dropped`.

---

## Gate 3 — Drug: collapse salt forms to parent  `[drug seed]`

**Issue:** `erlotinib` and `erlotinib-hydrochloride` ship as two pages; the salt is
strictly inferior (1 CIViC vs 209) but competes for the query. Parent↔child is known to
biobtree (`chembl_moleculeparent`) and the anchor resolves `parent_chembl`.

**Gate:** in the same seed-refinement pass, if a molecule's `parent_chembl` is also in
the build, **drop the salt, keep the parent**. Over-exclusion risk ≈ 0 (the canonical
parent page remains; nothing is lost). Web layer: emit `rel=canonical` from the salt slug
→ parent (later, theme lane). Cheap pre-check: name-based salt-suffix strip catches the
common ones before the resolve pass.

---

## Gate 4 — Disease subtypes: canonical-link, do NOT drop  `[render, not admission]`

**Issue:** over-granular numbered subtypes — `lynch-syndrome-1`, `fanconi-anemia-
complementation-group-a` — compete with the clinically-canonical parent, and carry
genetically-backwards counts (the #3 edge bug, separate).

**Decision (highest over-exclusion risk → most conservative):** numbered OMIM-locus
subtypes are **real diseases**, not junk — dropping them risks losing genuine
gene/phenotype distinctions. So **do NOT admission-gate them.** Instead (render-time /
web): keep the parent as the citable page, emit `rel=canonical` subtype→parent, and show
a subtypes table on the parent. This is a dedup/SEO concern, not an inclusion one.

---

## Not admission gates (tracked separately — edge/join/render correctness)

The audit's other HIGH issues are NOT solved by inclusion filters and need their own
fixes (these corrupt facts/edges, incl. the JSON-LD mesh + answer-first lead):
- **#2 disease→trial word-match join** fabricates trials (vici-syndrome "1,156 trials").
- **#3 ungated mesh edges**: gene→disease from raw OMIM/Mondo∩has-page (pseudogene→Fanconi);
  gene→drug from promiscuous bioactivity (TP53→colchicine). Gate on GenCC/ClinGen and
  CIViC/mechanism respectively.
- **#4 non-coding genes inherit the overlapping gene's biology** (TTN-AS1 → TTN's ClinVar).
  Suppress positional ClinVar/trial/disease blocks for non-coding biotypes + add an
  affirmative "non-coding, not a drug target" line.

These run **alongside** the admission gates (Phase A′), because they poison the mesh +
structured data we pitch to AI — see 05 §6.

---

## Implementation order (tomorrow)
1. Gate 1 (disease characteristic subtree) — `parse_obo` enhancement + re-emit disease seed.
2. Gates 2+3 (drug seed-refinement: pharmacological-value union + salt→parent) — one pass.
3. Audit both `dropped.txt` lists; confirm pure junk; tighten/loosen.
4. Re-run the curated/dense slice; diff mesh + leads vs the `dense-test` baseline.
5. Then Phase A′ edge/join correctness (#2/#3/#4), then Phase B polish (#6–#13).
