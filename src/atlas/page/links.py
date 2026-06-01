"""Cross-entity Atlas link mesh — resolve an entity reference (by ID or name)
to the internal `/atlas/<type>/<slug>/` URL of its Atlas page, IFF that page
has been built.

Design (see docs/research/05 §6 + NEXT.md cross-entity work):
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
import json
import os
import re

# Process-global manifest. {entity_type: {resolvable_key: slug}}.
_MANIFEST = {"gene": {}, "disease": {}, "drug": {}}
_LOADED_FROM = None


def _manifest_path(dist_root):
    return os.path.join(dist_root, "atlas", "manifest.json")


def _norm(s):
    """Normalize a name for keying — lowercase, collapse non-alnum to spaces."""
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9]+", " ", (s or "").lower())).strip()


def load(dist_root):
    """Load the manifest into process-global state. Idempotent per dist_root;
    safe to call at the start of every page build. Missing file → empty mesh
    (resolvers return None → plain text)."""
    global _MANIFEST, _LOADED_FROM
    path = _manifest_path(dist_root)
    try:
        with open(path) as f:
            data = json.load(f)
        _MANIFEST = {"gene": data.get("gene") or {},
                     "disease": data.get("disease") or {},
                     "drug": data.get("drug") or {}}
    except (FileNotFoundError, json.JSONDecodeError):
        _MANIFEST = {"gene": {}, "disease": {}, "drug": {}}
    _LOADED_FROM = dist_root
    return _MANIFEST


def reset():
    """Clear the mesh (tests / fresh runs)."""
    global _MANIFEST, _LOADED_FROM
    _MANIFEST = {"gene": {}, "disease": {}, "drug": {}}
    _LOADED_FROM = None


def upsert(dist_root, entity, slug, id_keys=(), name_keys=()):
    """Read-modify-write the manifest with one entity's resolvable keys → slug.
    `id_keys` are matched verbatim (HGNC symbol, Mondo ID, ChEMBL ID, …);
    `name_keys` are stored normalized (canonical name + synonyms/aliases).
    Also updates the in-memory mesh so a build sees its own entry."""
    if not slug:
        return
    path = _manifest_path(dist_root)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    try:
        with open(path) as f:
            data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        data = {}
    for k in ("gene", "disease", "drug"):
        data.setdefault(k, {})
    bucket = data[entity]
    for k in id_keys:
        if k:
            bucket[str(k)] = slug
    for k in name_keys:
        nk = _norm(k)
        if nk:
            bucket[nk] = slug
    from atlas.atomicio import write_json
    write_json(path, data, indent=0, sort_keys=True)   # atomic — see atomicio
    # keep the live mesh current
    _MANIFEST.setdefault(entity, {}).update(bucket)


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


def gene_url(symbol=None, hgnc_id=None):
    return _url("gene", _lookup("gene", symbol, hgnc_id))


def disease_url(mondo_id=None, name=None):
    return _url("disease", _lookup("disease", mondo_id, name))


def drug_url(chembl_id=None, name=None):
    return _url("drug", _lookup("drug", chembl_id, name))


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
    return n or (name or "")


def link_csv(cell, resolver):
    """Link each token in a comma-joined cell via resolver(token) -> url|None.
    Unmatched tokens stay plain text. For cells like "ABL1, DDR1, KIT" or a
    drug-name list — turns the resolvable ones into internal links in place."""
    if not cell:
        return cell or ""
    parts = [p.strip() for p in str(cell).split(",")]
    return ", ".join(maybe_link(p, resolver(p)) for p in parts if p)


def related_targets(entity_type, bundle):
    """Resolve a page's cross-entity references to BUILT Atlas targets.
    Returns {"Genes": [(label, path)], "Diseases": [...], "Drugs": [...]},
    deduped by path (e.g. "/atlas/gene/TP53/"). Single source of truth for both
    related_block (markdown) and the JSON-LD cross-entity edges."""
    from atlas.civic import therapy_label
    groups = {"Genes": [], "Diseases": [], "Drugs": []}
    seen = set()

    def add(grp, label, url):
        if url and label and url not in seen:
            seen.add(url)
            groups[grp].append((str(label), url))

    if entity_type == "gene":
        b10, b12 = bundle.get("10") or {}, bundle.get("12") or {}
        for m in (b12.get("mondo") or []):
            add("Diseases", m.get("name"), disease_url(mondo_id=m.get("id"), name=m.get("name")))
        for g in (b12.get("gencc") or []):                # name-tier
            add("Diseases", g.get("disease"), disease_url(name=g.get("disease")))
        for c in (b12.get("clingen_validity") or []):     # name-tier
            add("Diseases", c.get("disease"), disease_url(name=c.get("disease")))
        for m in (b10.get("molecules") or []):
            add("Drugs", _drug_display(m.get("name")), drug_url(chembl_id=m.get("id"), name=m.get("name")))
        for r in (b10.get("civic_evidence") or []):        # name-tier
            th = therapy_label(r.get("therapy"))
            add("Drugs", _drug_display(th), drug_url(name=th))
            add("Diseases", r.get("disease"), disease_url(name=r.get("disease")))
    elif entity_type == "disease":
        b4, b5, b10, b13 = (bundle.get(k) or {} for k in ("4", "5", "10", "13"))
        for g in (b5.get("genes") or []):
            add("Genes", g.get("symbol"), gene_url(symbol=g.get("symbol"), hgnc_id=g.get("hgnc_id")))
        for g in (b4.get("somatic_driver_genes") or []):
            add("Genes", g.get("symbol"), gene_url(symbol=g.get("symbol")))
        for d in (b13.get("trial_drugs") or []):
            add("Drugs", _drug_display(d.get("name") or d.get("molecule_id")),
                drug_url(chembl_id=d.get("molecule_id"), name=d.get("name")))
        for d in (b10.get("drugs") or []):
            add("Drugs", _drug_display(d.get("name") or d.get("id")),
                drug_url(chembl_id=d.get("id"), name=d.get("name")))
        for r in (b13.get("civic_evidence") or []):        # name-tier
            th = therapy_label(r.get("therapy"))
            add("Drugs", _drug_display(th), drug_url(name=th))
    elif entity_type == "drug":
        b2, b4, b7, b10 = (bundle.get(k) or {} for k in ("2", "4", "7", "10"))
        for t in (b2.get("primary_targets") or []):
            add("Genes", t.get("gene_symbol"), gene_url(symbol=t.get("gene_symbol"), hgnc_id=t.get("hgnc_id")))
        for i in (b4.get("indications") or []):
            add("Diseases", i.get("name"), disease_url(mondo_id=i.get("mondo_id"), name=i.get("name")))
        for r in (b7.get("related_molecules") or []):
            add("Drugs", _drug_display(r.get("name")), drug_url(name=r.get("name")))
        for r in (b10.get("civic_evidence") or []):        # name-tier
            add("Diseases", r.get("disease"), disease_url(name=r.get("disease")))

    return groups


def related_block(entity_type, bundle):
    """The "## Related Atlas pages" markdown section (page end) surfacing the
    mesh as a scannable block. Elides when nothing is built."""
    groups = related_targets(entity_type, bundle)
    lines = []
    for grp in ("Genes", "Diseases", "Drugs"):
        items = groups[grp]
        if not items:
            continue
        shown = items[:12]
        row = ", ".join(maybe_link(lbl, url) for lbl, url in shown)
        extra = f" (+{len(items) - 12} more)" if len(items) > 12 else ""
        lines.append(f"- **{grp}:** {row}{extra}")
    return ("## Related Atlas pages\n\n" + "\n".join(lines)) if lines else ""
