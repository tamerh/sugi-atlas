"""§1 — disease_ids: federated identifier set (Mondo + EFO + MeSH + OMIM +
Orphanet) + canonical name + the per-dataset xref count table.

NEW collector (not a gene fanout) — operates directly off the DiseaseAnchors
record, which already pre-resolved the ID set during resolve(). The work
here is mostly shaping; no new biobtree calls beyond the anchor."""
import re
from atlas.section import Section
from atlas.biobtree import map_all, entry

# Cap stored ontology-family lists (children/siblings). A disease under a broad
# Mondo parent ("hereditary disease") has ~1,900 co-subtypes; storing/rendering
# all of them bloats the page (a 174 KB line). 300 keeps the list useful while
# trimming the extreme tail — true counts are kept separately.
_FAMILY_CAP = 300

# biobtree/Orphanet returns the HPO frequency band with the percentages in
# descending order — "Very frequent (99-80%)". Normalize to the conventional
# low→high reading ("80-99%"). Generic (any "N-M%" → ascending), so it also
# fixes Frequent (79-30%), Occasional (29-5%), etc.
_FREQ_RANGE = re.compile(r"(\d+)\s*-\s*(\d+)\s*%")


def _ascending_freq(s):
    if not s:
        return s
    return _FREQ_RANGE.sub(
        lambda m: "{}-{}%".format(*sorted((int(m.group(1)), int(m.group(2))))), s)


# A MeSH note is only a usable clinical description if it reads like prose, not a
# MeSH indexing/etiology fragment. Even MAIN descriptors carry junk here:
# "AOMS1 not included" (entry-term indexing note), "mutation in MPZ" (etiology
# stub). Real scope notes are full sentences (hypertension's is 200+ chars).
_NOTE_JUNK = re.compile(r"\bnot included\b|^[\w;.\-\s]{0,30}mutation in\b", re.I)


def _looks_like_description(note):
    note = (note or "").strip()
    return len(note) >= 40 and not _NOTE_JUNK.search(note)


CHAINS   = (">>mondo>>efo", ">>mondo>>mesh", ">>mondo>>mim", ">>mondo>>orphanet",
            ">>mondo>>doid", ">>mondo>>sctid", ">>mondo>>umls", ">>mondo>>ncit",
            ">>mondo>>medgen", ">>mondo>>icd10cm", ">>mondo>>icd11",
            ">>mondo>>gard", ">>mondo>>meddra", ">>mondo>>nord",
            ">>mondo>>uberon", ">>mondo>>mondochild",
            ">>mim>>hpo", ">>mondo>>hpo")
DATASETS = ("mondo", "efo", "mesh", "mim", "orphanet",
            "doid", "sctid", "umls", "ncit", "medgen",
            "icd10cm", "icd11", "gard", "meddra", "nord", "uberon", "mondochild", "hpo")

def collect(a):
    # HPO clinical features. Orphanet's curated list (with frequency bands) is the
    # primary source, but many diseases — especially Mendelian/rare ones with an
    # OMIM entry but no Orphanet phenotype list — carry HPO annotations only via
    # OMIM (and sometimes Mondo). Pull those too and merge (deduped by HPO id), so
    # the Clinical-features section isn't falsely empty: e.g. MONDO:0009056 ("cutis
    # verticis gyrata and intellectual disability") has no Orphanet phenotypes but
    # 2 HPO features via OMIM 219300 (intellectual disability, cutis gyrata).
    oa = a.orphanet_attrs or {}
    # New dicts (don't mutate the shared anchor) with the frequency band in
    # ascending percentage order.
    phenotypes = [{**p, "frequency": _ascending_freq(p.get("frequency"))}
                  for p in (oa.get("phenotypes") or [])]
    seen_hpo = {p.get("hpo_id") for p in phenotypes if p.get("hpo_id")}

    def _add_hpo(term_id, chain):
        # OMIM/Mondo HPO carry id + name + definition but no frequency.
        try:
            for r in map_all(term_id, chain):
                hid = r.get("id")
                if hid and hid not in seen_hpo:
                    seen_hpo.add(hid)
                    phenotypes.append({"hpo_id": hid, "hpo_term": r.get("name"),
                                       "frequency": "", "frequency_value": 0,
                                       "definition": r.get("definition")})
        except Exception:
            pass

    for oid in (a.omim_ids or ()):
        _add_hpo(oid, ">>mim>>hpo")
    if a.mondo_id:
        _add_hpo(a.mondo_id, ">>mondo>>hpo")
    # Orphanet (with frequency) sorts first; OMIM/Mondo extras (no frequency) after.
    phenotypes.sort(key=lambda p: float(p.get("frequency_value") or 0), reverse=True)

    # Mondo ontology family — the local neighborhood that drives the Disease
    # family section: ancestors (breadcrumb, walked up the primary is-a chain),
    # children (this term's subtypes), siblings (the parent's other children).
    # Rows carry id+name (no per-term resolution); calls are cheap and heavily
    # cached (the upper ontology is shared), and diseases aren't fanned for §1.
    ancestors = []      # nearest-first: [parent, grandparent, …, root]
    children = []
    siblings = []
    if a.mondo_id:
        seen, cur = set(), a.mondo_id           # walk up, bounded + cycle-guarded
        for _ in range(10):
            try:
                pr = map_all(cur, ">>mondo>>mondoparent")
            except Exception:
                break
            if not pr:
                break
            pid, pname = pr[0].get("id"), pr[0].get("name")
            if not pid or pid in seen:
                break
            seen.add(pid)
            ancestors.append({"id": pid, "name": pname})
            cur = pid
        try:
            children = [{"id": r.get("id"), "name": r.get("name")}
                        for r in map_all(a.mondo_id, ">>mondo>>mondochild") if r.get("id")]
        except Exception:
            pass
        if ancestors:                           # siblings = parent's children − self
            try:
                siblings = [{"id": r.get("id"), "name": r.get("name")}
                            for r in map_all(ancestors[0]["id"], ">>mondo>>mondochild")
                            if r.get("id") and r.get("id") != a.mondo_id]
            except Exception:
                pass
    parent = ancestors[0] if ancestors else None
    child_count, sibling_count = len(children), len(siblings)
    children, siblings = children[:_FAMILY_CAP], siblings[:_FAMILY_CAP]

    # MeSH scope note — a curated clinical-description paragraph. HPO (a rare/
    # Mendelian phenotype ontology) leaves the Clinical-features zone empty for
    # common multifactorial diseases (hypertension, asthma, T2D), but those carry
    # a MeSH descriptor whose scope_note IS that description. We already resolve
    # the MeSH id (for the xref); the note lives only on the entry, so fetch it.
    # A Mondo can map to several MeSH descriptors — prefer the one whose name
    # matches this disease (so we don't surface a tangential descriptor's note),
    # falling back to the first descriptor that carries a note.
    # ONLY main MeSH Descriptors (D-numbers) carry a real clinical scope note.
    # Rare diseases map to Supplementary Concept Records (C-numbers) whose same
    # field holds an indexing/etiology note instead ("Optb2 not included",
    # "mutation in prohormone convertase-1") — NOT a description, so skip those.
    def _norm(s):
        return re.sub(r"[^a-z0-9]", "", (s or "").lower())
    want = _norm(a.canonical_name or a.name)
    mesh_scope_note, mesh_note_fallback = "", ""
    for mid in (a.mesh_ids or ()):
        try:
            mat = (entry(mid, "mesh").get("Attributes") or {}).get("Mesh") or {}
        except Exception:
            continue
        if str(mat.get("is_supplementary")).lower() == "true":
            continue                          # Supplementary Concept Record — not a description
        note = (mat.get("scope_note") or "").strip()
        if not _looks_like_description(note):  # reject indexing/etiology fragments
            continue
        if _norm(mat.get("descriptor_name")) == want:
            mesh_scope_note = note            # exact-name descriptor wins outright
            break
        mesh_note_fallback = mesh_note_fallback or note
    mesh_scope_note = mesh_scope_note or mesh_note_fallback

    bundle = {
        "section": "01_disease_ids",
        "name": a.name,
        "canonical_name": a.canonical_name,
        "synonyms": list(a.synonyms or ()),
        "mondo_id": a.mondo_id,
        "efo_id": a.efo_id,
        "mesh_ids": list(a.mesh_ids),
        "omim_ids": list(a.omim_ids),
        "orphanet_ids": list(a.orphanet_ids),
        # Cross-ontology xrefs from Mondo OBO ingest. {prefix: [ids,...]} —
        # only keys with data present.
        "obo_xrefs": {k: list(v) for k, v in (a.obo_xrefs or {}).items()},
        # UBERON anatomy ids (drives schema.org `associatedAnatomy`).
        "anatomy_uberon_ids": list(a.anatomy_uberon_ids or ()),
        # Orphanet primary-entry attrs — resolved once at anchor time.
        # prevalences = multi-geography epidemiology rows (drives JSON-LD
        # `epidemiology`). phenotypes = HPO list with both label and
        # numeric frequency_value (drives JSON-LD `signOrSymptom`).
        # Empty for non-rare-disease conditions (most cancers, common dz).
        "orphanet_name": oa.get("name") or "",
        "orphanet_disorder_type": oa.get("disorder_type") or "",
        "prevalences": list(oa.get("prevalences") or []),
        "phenotypes": phenotypes,
        "phenotype_count": len(phenotypes),   # merged Orphanet + OMIM/Mondo HPO
        # Orphanet curated clinical description + inheritance/onset (biobtree
        # 2026-06-13 refresh). Definition is the license-clean, disease-specific
        # clinical paragraph; complements the MeSH scope note (common diseases).
        "orphanet_definition": (oa.get("definition") or "").strip(),
        "orphanet_inheritance": list(oa.get("inheritance") or []),
        "orphanet_onset": list(oa.get("onset") or []),
        "mesh_scope_note": mesh_scope_note,
        "is_cancer": a.is_cancer,
        "child_count": child_count,
        "sibling_count": sibling_count,
        "parent": parent,
        "ancestors": ancestors,
        "children": children,
        "siblings": siblings,
        "xref_counts": dict(a.xref_counts),
    }
    return bundle

SECTION = Section(
    id="1", name="disease_ids",
    description=("Federated disease identifiers (Mondo, EFO, MeSH, OMIM, "
                 "Orphanet) + canonical Mondo name + per-dataset xref counts "
                 "+ Orphanet epidemiology (prevalences) and clinical features "
                 "(HPO phenotype list with frequencies)."),
    needs=("mondo_id", "canonical_name", "efo_id", "mesh_ids", "omim_ids",
           "orphanet_ids", "orphanet_attrs", "obo_xrefs", "anatomy_uberon_ids",
           "xref_counts", "is_cancer"),
    produces=("mondo_id", "canonical_name", "synonyms", "efo_id", "mesh_ids", "omim_ids",
              "orphanet_ids", "obo_xrefs", "anatomy_uberon_ids",
              "orphanet_name", "orphanet_disorder_type",
              "prevalences", "phenotypes", "phenotype_count", "mesh_scope_note",
              "orphanet_definition", "orphanet_inheritance", "orphanet_onset",
              "child_count", "parent", "ancestors", "children", "siblings",
              "xref_counts", "is_cancer"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
