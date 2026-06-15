"""Cross-entity Atlas link mesh — resolve an entity reference (by ID or name)
to the internal `/atlas/<type>/<slug>/` URL of its Atlas page, IFF that page
has been built.

Design (see docs/internal/research/05 §6 + NEXT.md cross-entity work):
  - A manifest at `<dist>/atlas/manifest.json` maps every resolvable key
    (canonical IDs, child/parent ChEMBLs, normalized names + synonyms) → slug,
    per entity type. It is written by the pipeline at build time (run_gene /
    run_disease / run_drug), which knows the *real* slug (incl. disease
    slug-overrides) and every ID from the bundle — so it's never recomputed.
  - Resolution is process-global state (like atlas.biobtree.CALLS): the
    pipeline calls `load(dist_root)` before rendering; renderers + JSON-LD
    builders then call `gene_url()/disease_url()/drug_url()` freely.
  - Until a target entity is in the manifest, resolvers return None and the
    caller renders plain text (current behavior). With "render-time, 2× regen"
    the mesh is complete after a second full build; the manifest persists, so
    each run only adds.

Gene slug == HGNC symbol. Disease/drug slug == the slug written at build time.
"""
import contextlib
import json
import os
import re

try:
    import fcntl
except ImportError:                      # non-POSIX — batch is single-writer anyway
    fcntl = None


@contextlib.contextmanager
def _manifest_lock(path):
    """Serialize manifest read-modify-write across processes (audit #15: ad-hoc
    parallel per-entity runs lost updates — a concurrent reader/writer clobbered
    each other; the production batch is single-writer and unaffected). Advisory
    flock on a sidecar; no-op where fcntl is unavailable."""
    if fcntl is None:
        yield
        return
    with open(path + ".lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lk, fcntl.LOCK_UN)


# Process-global manifest. {entity_type: {resolvable_key: slug}}.
_MANIFEST = {"gene": {}, "disease": {}, "drug": {}}
# Destination canonical name per slug (audit #13): a name_key can resolve to a
# page whose canonical name differs (synonym "schizoaffective disorder" →
# /schizophrenia/). This lets a link render the page it ACTUALLY points to.
_CANON = {"gene": {}, "disease": {}, "drug": {}}
_LOADED_FROM = None

# Reverse-edge index {target_url: [[src_label, src_url, src_type, group], …]} —
# incoming cross-entity edges, so a relationship asserted from one side (TP53's
# CIViC names Venetoclax) is navigable from the other (Venetoclax ← TP53). Built
# in the batch merge phase (corpus builds only); empty for single-page builds.
_REVERSE = {}
_REVERSE_FROM = None

# Disease→drug indication index {disease_url: [{name, url, max_phase}, …]} —
# drugs whose labelled/late-stage (phase ≥3) ChEMBL indication is this disease.
# A disease-DIRECT therapeutic edge (not gene-cohort-mediated), so non-molecular
# diseases (autoimmune/clinical-only) still surface real registered drugs. Built
# in the merge phase (corpus builds only); empty for single-page builds.
_INDICATIONS = {}
_INDICATIONS_FROM = None


def _manifest_path(dist_root):
    return os.path.join(dist_root, "atlas", "manifest.json")


def _norm(s):
    """Normalize a name for keying — lowercase, collapse non-alnum to spaces."""
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9]+", " ", (s or "").lower())).strip()


def load(dist_root):
    """Load the manifest into process-global state. Idempotent per dist_root;
    safe to call at the start of every page build. Missing file → empty mesh
    (resolvers return None → plain text)."""
    global _MANIFEST, _CANON, _LOADED_FROM
    path = _manifest_path(dist_root)
    try:
        with open(path) as f:
            data = json.load(f)
        _MANIFEST = {"gene": data.get("gene") or {},
                     "disease": data.get("disease") or {},
                     "drug": data.get("drug") or {}}
        canon = data.get("canon") or {}
        _CANON = {k: canon.get(k) or {} for k in ("gene", "disease", "drug")}
    except (FileNotFoundError, json.JSONDecodeError):
        _MANIFEST = {"gene": {}, "disease": {}, "drug": {}}
        _CANON = {"gene": {}, "disease": {}, "drug": {}}
    _LOADED_FROM = dist_root
    return _MANIFEST


def load_reverse(dist_root):
    """Load the reverse-edge index (incoming cross-entity edges). Cached per
    dist_root so a worker reads the file once, not per page. Missing file →
    empty (single-page builds have no reverse mesh — it degrades silently)."""
    global _REVERSE, _REVERSE_FROM
    if _REVERSE_FROM == dist_root:
        return _REVERSE
    path = os.path.join(dist_root, "atlas", "reverse_edges.json")
    try:
        with open(path) as f:
            _REVERSE = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        _REVERSE = {}
    _REVERSE_FROM = dist_root
    return _REVERSE


def load_indications(dist_root):
    """Load the disease→drug indication index (merge sidecar). Cached per
    dist_root. Missing file → empty (single-page builds have no index)."""
    global _INDICATIONS, _INDICATIONS_FROM
    if _INDICATIONS_FROM == dist_root:
        return _INDICATIONS
    path = os.path.join(dist_root, "atlas", "indicated_drugs.json")
    try:
        with open(path) as f:
            _INDICATIONS = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        _INDICATIONS = {}
    _INDICATIONS_FROM = dist_root
    return _INDICATIONS


def indicated_drugs(dist_root, slug):
    """Drugs whose phase ≥3 ChEMBL indication is this disease — [{name, url,
    max_phase}], approved-first. Read at render time and injected into the
    disease bundle (`_indicated_drugs`) for the Therapeutics section."""
    return load_indications(dist_root).get(f"/atlas/disease/{slug}/") or []


def reset():
    """Clear the mesh (tests / fresh runs)."""
    global _MANIFEST, _CANON, _REVERSE, _LOADED_FROM, _REVERSE_FROM
    global _INDICATIONS, _INDICATIONS_FROM
    _MANIFEST = {"gene": {}, "disease": {}, "drug": {}}
    _CANON = {"gene": {}, "disease": {}, "drug": {}}
    _REVERSE = {}
    _LOADED_FROM = None
    _REVERSE_FROM = None
    _INDICATIONS = {}
    _INDICATIONS_FROM = None


def upsert(dist_root, entity, slug, id_keys=(), name_keys=(), canonical=None):
    """Read-modify-write the manifest with one entity's resolvable keys → slug.
    `id_keys` are matched verbatim (HGNC symbol, Mondo ID, ChEMBL ID, …);
    `name_keys` are stored normalized (canonical name + synonyms/aliases).
    Also updates the in-memory mesh so a build sees its own entry."""
    if not slug:
        return
    path = _manifest_path(dist_root)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    from atlas.atomicio import write_json
    # Hold the lock across the whole read-modify-write so concurrent upserts don't
    # lose each other's keys (atomic write alone only prevents a torn file).
    with _manifest_lock(path):
        try:
            with open(path) as f:
                data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            data = {}
        for k in ("gene", "disease", "drug"):
            data.setdefault(k, {})
        data.setdefault("canon", {}).setdefault(entity, {})
        bucket = data[entity]
        for k in id_keys:
            if k:
                bucket[str(k)] = slug
        for k in name_keys:
            nk = _norm(k)
            if nk:
                bucket[nk] = slug
        if canonical:
            data["canon"][entity][slug] = canonical
        write_json(path, data, indent=0, sort_keys=True)   # atomic — see atomicio
    # keep the live mesh current
    _MANIFEST.setdefault(entity, {}).update(bucket)
    if canonical:
        _CANON.setdefault(entity, {})[slug] = canonical


def _lookup(entity, *keys):
    bucket = _MANIFEST.get(entity) or {}
    for k in keys:
        if k is None:
            continue
        if k in bucket:                       # verbatim ID
            return bucket[k]
        nk = _norm(k)
        if nk and nk in bucket:               # normalized name
            return bucket[nk]
    return None


def _url(entity, slug):
    return f"/atlas/{entity}/{slug}/" if slug else None


_URL_RE = re.compile(r"^/atlas/(gene|disease|drug)/([^/]+)/$")


def canonical_label(url):
    """The destination page's canonical name for an /atlas/<entity>/<slug>/ URL,
    or None if unknown (audit #13: render the page a link actually points to)."""
    m = _URL_RE.match(url or "")
    if not m:
        return None
    return (_CANON.get(m.group(1)) or {}).get(m.group(2))


def gene_url(symbol=None, hgnc_id=None):
    return _url("gene", _lookup("gene", symbol, hgnc_id))


def disease_url(mondo_id=None, name=None):
    return _url("disease", _lookup("disease", mondo_id, name))


def drug_url(chembl_id=None, name=None):
    return _url("drug", _lookup("drug", chembl_id, name))


def uniprot_url(acc):
    """UniProtKB entry URL for an accession, else None. Shared so UniProt
    accessions link consistently across page types (the drug bioactivity table
    linked them; gene/disease cohort tables showed them as plain text)."""
    acc = (acc or "").strip()
    return f"https://www.uniprot.org/uniprotkb/{acc}/entry" if acc else None


def variant_link(v):
    """dbSNP URL for an rs-prefixed variant id, else None. Shared by every page
    type that renders variants (gene/drug PGx tables, disease GWAS/variant
    tables) so rsIDs link consistently. Non-rs identifiers (star-alleles,
    HGVS, ClinVar names) stay plain text."""
    v = (v or "").strip()
    return f"https://www.ncbi.nlm.nih.gov/snp/{v}" if v.startswith("rs") else None


def maybe_link(text, url):
    """Markdown link to an internal Atlas page when `url` is set, else plain
    text. `text` is shown verbatim (already the display name/id)."""
    if not text:
        return text or ""
    return f"[{text}]({url})" if url else str(text)


# Trailing salt/hydrate words stripped from a drug *display* name (the link
# target already resolves via id/normalized-name; this is cosmetic).
_SALT_SUFFIXES = (" anhydrous", " mesylate", " hydrochloride", " dihydrochloride",
                  " sulfate", " sulphate", " maleate", " citrate", " calcium",
                  " sodium", " potassium", " besylate", " hydrobromide", " succinate")


def _drug_display(name):
    n = (name or "").strip()
    low = n.lower()
    changed = True
    while changed:
        changed = False
        for suf in _SALT_SUFFIXES:
            if low.endswith(suf):
                n, low = n[: -len(suf)].rstrip(), low[: -len(suf)].rstrip()
                changed = True
    # De-SHOUT (audit #12): ChEMBL names are all-caps ("CISPLATIN"); the mesh
    # label must match the de-SHOUTed drug-page title, not shout in the
    # cross-entity block.
    from atlas.render_common import display_name
    return display_name(n or (name or ""))


def link_csv(cell, resolver):
    """Link each token in a comma-joined cell via resolver(token) -> url|None.
    Unmatched tokens stay plain text. For cells like "ABL1, DDR1, KIT" or a
    drug-name list — turns the resolvable ones into internal links in place."""
    if not cell:
        return cell or ""
    parts = [p.strip() for p in str(cell).split(",")]
    return ", ".join(maybe_link(p, resolver(p)) for p in parts if p)


def _gene_evidence_score(g):
    """Disease-specificity rank for a cohort gene (audit #10). Curated
    gene–disease validity (GenCC) and clinical evidence (CIViC) weigh more than
    common-variant GWAS / ClinVar overlap. Stable within a score band, so the
    upstream cohort order is preserved for ties."""
    ev = g.get("evidence") or {}
    return ((2 if ev.get("gencc") else 0)
            + (2 if ev.get("civic_evidence") else 0)
            + (1 if ev.get("clinvar") else 0)
            + (1 if ev.get("gwas") else 0))


_NCRNA_MESH_CAP = 25   # bound a non-coding RNA gene's §14 disease/drug mesh edges


def related_targets(entity_type, bundle):
    """Resolve a page's cross-entity references to BUILT Atlas targets.
    Returns {"Genes": [(label, path)], "Diseases": [...], "Drugs": [...]},
    deduped by path (e.g. "/atlas/gene/TP53/"). Single source of truth for both
    related_block (markdown) and the JSON-LD cross-entity edges."""
    from atlas.civic import therapy_label
    groups = {"Genes": [], "Diseases": [], "Trial diseases": [], "Drugs": [],
              "ncRNA diseases": [], "ncRNA drugs": [], "Pharmacogenes": []}
    seen = set()

    def add(grp, label, url):
        if url and label and url not in seen:
            seen.add(url)
            groups[grp].append((str(label), url))

    if entity_type == "gene":
        # CURATED edges only (correctness #3). Deliberately NOT b12.mondo (raw
        # OMIM/Mondo∩has-page → pseudogene→Fanconi via a shared MIM) nor
        # b10.molecules (phase≥1 bioactivity = off-target/promiscuous →
        # TP53→colchicine mislabeled schema:target). gene→disease comes from
        # GenCC/ClinGen + CIViC; gene→drug from CIViC therapies.
        b10, b12 = bundle.get("10") or {}, bundle.get("12") or {}
        for g in (b12.get("gencc") or []):
            add("Diseases", g.get("disease"), disease_url(name=g.get("disease")))
        for c in (b12.get("clingen_validity") or []):
            add("Diseases", c.get("disease"), disease_url(name=c.get("disease")))
        for r in (b10.get("civic_evidence") or []):
            th = therapy_label(r.get("therapy"))
            add("Drugs", _drug_display(th), drug_url(name=th))
            add("Diseases", r.get("disease"), disease_url(name=r.get("disease")))
        # Non-coding RNA layer (§14) — curated ncRNA→disease (LncRNADisease/HMDD)
        # and ncRNA→drug (ncRNADrug), so a lncRNA/miRNA's Related block isn't empty
        # (its protein-based curated edges are scrubbed). Manifest-gated, so
        # free-text disease names that don't match a built page (and uncovered
        # drugs) simply drop; capped to keep the mesh bounded. Interaction partners
        # are RNA names that rarely match a gene slug — deliberately not meshed.
        # Distinct groups (not "Diseases"/"Drugs") so the forward labels stay
        # honest AND the reverse index has no mapping for them (REVERSE_LABEL omits
        # ("gene","ncRNA …")) — forward-only, so a drug page isn't polluted with
        # "MALAT1 is a CIViC biomarker" and a disease page isn't told an ncRNA is a
        # GenCC-curated associated gene.
        # Label each edge with the DESTINATION's canonical name (the ncRNA
        # dataset uses free-text synonyms — "Coronary Artery Disease" resolves to
        # /coronary-artery-disorder/; the link must show the page's real title).
        b14 = bundle.get("14") or {}
        for d in (b14.get("diseases") or [])[:_NCRNA_MESH_CAP]:
            u = disease_url(name=d.get("disease_name"))
            if u:
                add("ncRNA diseases", canonical_label(u) or d.get("disease_name"), u)
        for dr in (b14.get("drugs") or [])[:_NCRNA_MESH_CAP]:
            u = drug_url(name=dr.get("drug_name"))
            if u:
                add("ncRNA drugs", canonical_label(u) or _drug_display(dr.get("drug_name")), u)
    elif entity_type == "disease":
        b4, b5, b10, b13 = (bundle.get(k) or {} for k in ("4", "5", "10", "13"))
        # Rank cohort genes by disease-specific evidence (audit #10): the cohort
        # is HGNC-id-ordered, so the top-N shown was arbitrary gwas-only genes
        # (asthma → BCR/RUNX1/RYR1). Curated gene–disease validity (GenCC) and
        # clinical evidence (CIViC) outrank common-variant (GWAS) / ClinVar.
        # CIViC somatic drivers (b4) are curated, so they lead.
        for g in (b4.get("somatic_driver_genes") or []):
            add("Genes", g.get("symbol"), gene_url(symbol=g.get("symbol")))
        for g in sorted(b5.get("genes") or [], key=_gene_evidence_score, reverse=True):
            add("Genes", g.get("symbol"), gene_url(symbol=g.get("symbol"), hgnc_id=g.get("hgnc_id")))
        # Curated/real disease→drug only (#3 extension): title-validated trial
        # drugs + CIViC therapies. Deliberately NOT b10.drugs (bioactivity hits
        # on cohort targets → off-target junk: medulloblastoma "treated by"
        # Clotrimazole/Candesartan).
        for d in (b13.get("trial_drugs") or []):
            add("Drugs", _drug_display(d.get("name") or d.get("molecule_id")),
                drug_url(chembl_id=d.get("molecule_id"), name=d.get("name")))
        for r in (b13.get("civic_evidence") or []):        # name-tier
            th = therapy_label(r.get("therapy"))
            add("Drugs", _drug_display(th), drug_url(name=th))
    elif entity_type == "drug":
        b2, b4, b7, b10 = (bundle.get(k) or {} for k in ("2", "4", "7", "10"))
        # Curated targets only (#3 extension): GtoPdb-sourced. The ChEMBL-
        # bioactivity "primary" targets are assay hits, not real targets
        # (Salmeterol "targets" KDM4E/SMN1/TP53) — they'd poison the gene mesh
        # and the "Targeted by drugs" reverse.
        for t in (b2.get("primary_targets") or []):
            if (t.get("source") or "").lower() != "gtopdb":
                continue
            add("Genes", t.get("gene_symbol"), gene_url(symbol=t.get("gene_symbol"), hgnc_id=t.get("hgnc_id")))
        # Curated ChEMBL mechanism-of-action targets — the only target for RNA
        # therapeutics (inclisiran→PCSK9, patisiran→TTR). Curated (drug_mechanism),
        # so safe for the gene mesh, unlike the excluded bioactivity assay hits;
        # this is what makes the gene→RNA-drug reverse edge exist.
        for g in (b2.get("mechanism_genes") or []):
            add("Genes", g.get("gene_symbol"), gene_url(symbol=g.get("gene_symbol"), hgnc_id=g.get("hgnc_id")))
        # Phase ≥3 only (launch decision): all-phase indications surfaced
        # investigational junk as "indicated" (asthma ← Atorvastatin, a phase-3
        # trial). The disease side renders these as a dedicated "Drugs indicated"
        # section from the indication index, NOT the generic reverse block (the
        # ("drug","Diseases") reverse label is intentionally absent below).
        # Tier the drug→disease mesh on the per-indication `approved` flag (set in
        # s04 via atlas.indication): approved indications vs investigational trials,
        # kept navigable but distinctly labelled so "in trials for cancer" never
        # reads as "treats cancer" (aspirin), while a real approval logged at phase 3
        # (imatinib→CML, an anticancer drug) is NOT demoted. Phase ≤1 omitted.
        for i in (b4.get("indications") or []):
            ph = i.get("max_phase") or 0
            if ph < 2:
                continue
            url = disease_url(mondo_id=i.get("mondo_id"), name=i.get("name"))
            if i.get("approved"):
                add("Diseases", i.get("name"), url)
            else:
                add("Trial diseases", i.get("name"), url)
        for r in (b7.get("related_molecules") or []):
            add("Drugs", _drug_display(r.get("name")), drug_url(name=r.get("name")))
        for r in (b10.get("civic_evidence") or []):        # name-tier
            add("Diseases", r.get("disease"), disease_url(name=r.get("disease")))
        # Pharmacogenes — genes with a curated PharmGKB variant/clinical
        # association for THIS drug (PGx: affects metabolism/response, NOT a drug
        # target). From the §9 pharmacogenomics bundle. Reverses on gene pages as
        # "Pharmacogenomically associated drugs".
        b9 = bundle.get("9") or {}
        pgx = {r.get("gene") for r in (b9.get("variant_annotations") or [])}
        pgx |= {r.get("gene") for r in (b9.get("clinical_annotations") or [])}
        for sym in sorted(s for s in pgx if s):
            add("Pharmacogenes", sym, gene_url(symbol=sym))

    # Disease links resolve by name and a synonym can land on a differently-named
    # page (audit #13: "schizoaffective disorder" → /schizophrenia/). Relabel to
    # the destination's canonical name so the text matches the page it opens; drop
    # un-named pages whose only label is a raw ontology accession (never link a
    # bare MONDO id as a name).
    from atlas.render_common import is_ontology_id
    for grp in ("Diseases", "Trial diseases"):
        groups[grp] = [(lbl2, url) for lbl, url in groups[grp]
                       for lbl2 in (canonical_label(url) or lbl,)
                       if not is_ontology_id(lbl2)]

    return groups


def drug_target_landscape(bundle):
    """Reverse-index reads over a drug's curated (GtoPdb) target genes — the basis
    for the drug page's "target biology reach" and "competitive landscape":
      reach      — set of corpus DISEASE urls where a target is a cohort gene
      per_target — [(symbol, n_cohort_diseases, n_drugs_targeting_incl_self)]
    Corpus build only (returns empty when no reverse index is loaded — single-page
    rebuilds / unit tests degrade silently, like the other reverse-mesh reads)."""
    targets = related_targets("drug", bundle).get("Genes") or []
    reach, per_target = set(), []
    for sym, gurl in targets:
        diseases, drugs = set(), set()
        for src_label, src_url, src_type, group in (_REVERSE.get(gurl) or []):
            if src_type == "disease" and group == "Genes":
                diseases.add(src_url)
            elif src_type == "drug" and group == "Genes":
                drugs.add(src_url)            # includes this drug (it targets the gene)
        reach |= diseases
        per_target.append((sym, len(diseases), len(drugs)))
    return {"reach": reach, "n_targets": len(targets), "per_target": per_target}


# Reverse-edge labels: a forward edge from a source of type `src_type` in its
# `group` becomes, on the TARGET page, a group of the SOURCE entities under this
# label. ONLY directions whose SOURCE edge is curated are inverted — a biomarker
# is not a target (the #3 lesson). drug→gene is now GtoPdb-curated (the ChEMBL-
# bioactivity targets are filtered out in related_targets), so "Targeted by
# drugs" is safe to invert. Disease-as-source reverses stay omitted (redundant
# with the target's own cohort/indication edges).
REVERSE_LABEL = {
    ("gene", "Drugs"):    "Biomarker genes",     # genes whose variants associate this drug (CIViC) → on drug pages
    ("drug", "Genes"):    "Targeted by drugs",   # drugs that curatedly target this gene (GtoPdb) → on gene pages
    ("gene", "Diseases"): "Associated genes",    # genes asserting association with this disease (GenCC/ClinGen/CIViC) → on disease pages
    ("disease", "Genes"): "Disease cohort memberships",  # diseases whose associated-gene cohort includes this gene → on gene pages. COHORT MEMBERSHIP (GWAS/GenCC/ClinVar/CIViC), NOT a causal claim — shown in full (coexists with the gene's forward Associated diseases; the caption flags the overlap).
    ("drug", "Pharmacogenes"): "Pharmacogenomically associated drugs",  # drugs with a PharmGKB variant/clinical association for this gene → on gene pages (PGx, not targeting)
    # ("drug","Diseases") deliberately omitted: drug→disease indications are
    # surfaced as the dedicated "Drugs indicated for this disease" Therapeutics
    # section (render r_drugs_indicated, fed by the indicated_drugs.json index),
    # which carries the development phase the flat reverse block can't.
}
_REVERSE_ORDER = ["Biomarker genes", "Targeted by drugs", "Associated genes",
                  "Disease cohort memberships", "Pharmacogenomically associated drugs"]

# Italic clarifier rendered after a group label — kills the misread that a
# predictive/biomarker or cohort-membership edge is a target/causal claim.
_GROUP_CAPTION = {
    "Biomarker drugs (CIViC)":
        "drugs whose response is associated with variants in this gene — "
        "CIViC predictive evidence, not targeting",
    "Disease cohort memberships":
        "association, not causation — diseases whose associated-gene cohort lists "
        "this gene; a subset are also under Associated diseases",
    "Pharmacogenes":
        "genes with a PharmGKB pharmacogenomic association for this drug — they "
        "affect its metabolism/response, they are not its targets",
    "Pharmacogenomically associated drugs":
        "drugs whose metabolism/response is associated with variants in this gene "
        "(PharmGKB) — not drugs that target it",
}


def _reverse_groups(my_url, forward_urls):
    """Incoming edges for `my_url` from the reverse index, grouped by predicate,
    deduped against this page's own outgoing links (`forward_urls`) and across
    reverse groups. Returns [(label, [(src_label, src_url)])] in display order."""
    from atlas.render_common import is_ontology_id
    by_label = {}
    placed = set()                     # a url shows under one reverse label only
    # These reverse relationships are DISTINCT from any forward edge of the same
    # url, so they're shown even when the url is already a forward link: a drug can
    # be both a CIViC biomarker (forward) and a GtoPdb target (reverse); a disease
    # can be a curated association (forward) and list the gene in its cohort.
    _coexist = {"Targeted by drugs", "Disease cohort memberships"}
    for entry in (_REVERSE.get(my_url) or []):
        src_label, src_url, src_type, group = entry
        if not src_url or src_url == my_url or src_url in placed:
            continue
        lab = REVERSE_LABEL.get((src_type, group))
        if not lab:
            continue
        if lab not in _coexist and src_url in forward_urls:
            continue                   # already surfaced as a forward edge of this page
        # Render the destination's canonical name; skip un-named pages whose only
        # label is a raw ontology accession (e.g. a cohort disease that resolves to
        # MONDO:0014866) — never surface a bare id as a link label (test guard).
        label = canonical_label(src_url) or src_label
        if is_ontology_id(label):
            continue
        by_label.setdefault(lab, []).append((str(label), src_url))
        placed.add(src_url)
    return [(lab, by_label[lab]) for lab in _REVERSE_ORDER if by_label.get(lab)]


def related_block(entity_type, bundle, slug=None):
    """The "## Related Atlas pages" markdown section (page end) surfacing the
    mesh as a scannable block. Forward edges (this page → others) first, then
    reverse edges (others → this page, deduped). Elides when nothing is built."""
    groups = related_targets(entity_type, bundle)
    # Never link the page to itself (a drug's related_molecules can include
    # itself; a CIViC therapy can name the drug). Drop the self url from every
    # group.
    self_url = _url(entity_type, slug) if slug else None
    if self_url:
        for grp in groups:
            groups[grp] = [(l, u) for l, u in groups[grp] if u != self_url]
    # On disease pages the gene set is the associated-gene cohort (evidence-
    # ranked, audit #10), not a curated "most relevant" shortlist — label it
    # honestly so a polygenic-disease cohort isn't read as causal genes.
    if entity_type == "disease":
        label = {"Genes": "Cohort genes"}
    elif entity_type == "gene":                      # forward gene edges: name the relationship precisely
        label = {"Diseases": "Associated diseases", "Drugs": "Biomarker drugs (CIViC)",
                 # §14 ncRNA edges live in their own groups (never mixed with the
                 # curated CIViC/GenCC ones) so they're labelled honestly, not as
                 # "CIViC biomarker drugs".
                 "ncRNA diseases": "Associated diseases (ncRNA)",
                 "ncRNA drugs": "Associated drugs (ncRNA: response / resistance)"}
    elif entity_type == "drug":                      # tier drug→disease: approved vs investigational
        label = {"Diseases": "Indicated for", "Trial diseases": "In clinical trials for"}
    else:
        label = {}
    forward_urls = {url for items in groups.values() for _lbl, url in items}
    lines = []

    def _row(lbl_text, items):
        # Show ALL cross-links — these are navigable internal edges (the point of
        # the Related section) and matter for the web team's link graph / search.
        # No "(+N more)" dead-end text. Hub pages get long rows (max ~520 incoming
        # edges); the frontend collapses them. Forward groups are bounded by the
        # cohort/target caps; reverse groups can be large for hub genes/diseases.
        row = ", ".join(maybe_link(l, u) for l, u in items)
        cap = _GROUP_CAPTION.get(lbl_text)           # italic non-causality clarifier
        head = f"**{lbl_text}** *({cap})*:" if cap else f"**{lbl_text}:**"
        lines.append(f"- {head} {row}")

    for grp in ("Genes", "Diseases", "Trial diseases", "Drugs",
                "ncRNA diseases", "ncRNA drugs", "Pharmacogenes"):
        if groups.get(grp):
            _row(label.get(grp, grp), groups[grp])
    # Reverse edges (incoming) — corpus builds only; _REVERSE is empty otherwise.
    if slug:
        for rlabel, items in _reverse_groups(_url(entity_type, slug), forward_urls):
            _row(rlabel, items)
    # Always emit the canonical #related section (PAGE_CONTRACT — every section
    # present so the TOC is identical); placeholder when nothing is built yet.
    if not lines:
        lines = ["*No linked Atlas pages yet — the cross-entity mesh grows as "
                 "the corpus expands.*"]
    return "## Related Atlas pages {#related}\n\n" + "\n".join(lines)
