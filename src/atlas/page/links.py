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
    with open(path, "w") as f:
        json.dump(data, f, indent=0, sort_keys=True)
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
