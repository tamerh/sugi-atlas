"""gRPC -> REST-shape projection.

The gRPC API returns proto-native dicts:
    {results: {results: [{dataset, identifier, hgnc: {...}, count, ...}, ...]}}
    {result:  {dataset, identifier, hgnc: {...}, count, entries: [...], dataset_counts: [...]}}

The REST API returns two distinct compact shapes:
    /api/search -> {schema:"id|dataset|name|xref_count", data:["pipe|encoded|rows"]}
                  (FIXED schema, per-dataset name extraction.)
    /api/map    -> {schema:"id|<compact_fields_for_target_dataset>",
                    mappings:[{source:..., targets:["pipe|encoded|rows"]}]}

Source of truth: biobtree's own config (conf/source{1,2}.dataset.json +
xref{1,2}.dataset.json). Each dataset declares `id` + `compact_fields`.
We load these once at import.

Search-name extraction matches compact.go's ExtractSourceName — per-dataset
which proto attr field provides the canonical name (most are `name`; some
use `title`/`short_name`/`gene_name`/etc; hgnc/uniprot are `names[0]`).
"""
import json
import os
from pathlib import Path

# --- biobtree config loaders ----------------------------------------------

def _load_dataset_config():
    """({dataset_name: [compact_fields...]}, {dataset_id: name})

    Path overridable via ATLAS_BIOBTREE_CONF. Default: /data/biobtree/conf."""
    conf_dir = Path(os.environ.get("ATLAS_BIOBTREE_CONF", "/data/biobtree/conf"))
    fields, id2name = {}, {}
    for fname in ("source1.dataset.json", "source2.dataset.json",
                  "xref1.dataset.json", "xref2.optional.dataset.json"):
        path = conf_dir / fname
        if not path.exists():
            continue
        with path.open() as f:
            cfg = json.load(f)
        for name, entry in cfg.items():
            cf = (entry or {}).get("compact_fields") or ""
            if cf:
                fields[name] = [s.strip() for s in cf.split(",") if s.strip()]
            ds_id = (entry or {}).get("id")
            if ds_id is not None:
                try:
                    id2name[int(ds_id)] = name
                except (TypeError, ValueError):
                    pass
    return fields, id2name


COMPACT_FIELDS, DATASET_ID_TO_NAME = _load_dataset_config()


# --- proto attr-oneof resolution ------------------------------------------

# REST dataset name -> proto oneof field name (when they diverge).
# `ontology` covers the OBO datasets (go/eco/efo/uberon/cl/oba/pato/mondo).
_ATTR_KEY_OVERRIDES = {
    "hpo":     "hpo_attr",        # proto: HPOAttr hpo_attr=79
    "string":  "stringattr",      # proto: StringAttr stringattr=28
    "go":      "ontology", "eco": "ontology", "efo": "ontology",
    "uberon":  "ontology", "cl":  "ontology", "oba": "ontology",
    "pato":    "ontology", "mondo": "ontology",
}


def _attr_dict(xref: dict, dataset_name: str) -> dict:
    """The proto attr-oneof dict for this xref's dataset (or {} if absent)."""
    key = _ATTR_KEY_OVERRIDES.get(dataset_name, dataset_name)
    return xref.get(key) or {}


# --- search-name extraction (matches compact.go's ExtractSourceName) ------

# attr-field name to read as the canonical name, per dataset.
_NAME_FIELD = {
    # default for everything not listed: "name"
    "hgnc":            ("names", 0),       # names[0]
    "uniprot":         ("names", 0),
    "pdb":             "title",
    "pubchem":         "title",
    "interpro":        "short_name",
    "bgee":            "gene_name",
    "mesh":            "descriptor_name",
    "ctd":             "chemical_name",
    "clinical_trials": "brief_title",
    "scxa":            "description",
    "intogen":         "symbol",
    "civic_evidence":  "molecular_profile",
    "civic_assertion": "molecular_profile",
    "gtopdb_interaction": "target_name",
    # chembl is special — name lives under attr.molecule.name (nested oneof);
    # if missing, the name stays empty (matches compact.go).
    "chembl":          ("molecule", "name"),
}


def _extract_name(attr: dict, dataset_name: str) -> str:
    spec = _NAME_FIELD.get(dataset_name, "name")
    if isinstance(spec, str):
        return _flatten_scalar(attr.get(spec))
    if isinstance(spec, tuple):
        cur = attr
        for s in spec:
            if isinstance(s, int):
                cur = cur[s] if isinstance(cur, list) and len(cur) > s else None
            else:
                cur = cur.get(s) if isinstance(cur, dict) else None
            if cur is None:
                return ""
        return _flatten_scalar(cur)
    return ""


# --- computed-field overrides ---------------------------------------------

def _len(v): return len(v) if v else 0


def _f_alphafold_fraction(attr):
    pl = attr.get("plddt") or []
    if not pl:
        return ""
    vh = sum(1 for x in pl if float(x) >= 90)
    return f"{vh / len(pl):.4f}"


COMPUTED_FIELDS = {
    "biogrid": {
        "interaction_count": lambda a: str(_len(a.get("interactions"))),
        "unique_partners":   lambda a: str(len(set(
            i.get("partner_b", "") for i in (a.get("interactions") or [])))),
    },
    "alphafold": {
        "fraction_plddt_very_high": _f_alphafold_fraction,
    },
    "string": {
        "interaction_count": lambda a: str(_len(a.get("interactions"))),
    },
}


# --- proto-value -> REST-string ------------------------------------------

def _flatten_scalar(v) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, list):
        return ",".join(str(x) for x in v)
    return str(v)


def _escape_pipe(s: str) -> str:
    """Match compact.go's escapePipe — REST escapes `|` inside values as `\\|`."""
    return s.replace("|", "\\|")


# --- public projections ---------------------------------------------------

def _xref_compact_row(xref: dict) -> str:
    """Encode a single Xref as `id|<compact_fields>` (matches GetCompactRow).

    compact_fields list resolved from biobtree config by xref.dataset_name.
    Missing datasets -> just `id` (matches compact.go fallback)."""
    dataset_name = xref.get("dataset_name") or ""
    identifier = xref.get("identifier") or ""
    fields = COMPACT_FIELDS.get(dataset_name, [])
    parts = [_escape_pipe(identifier)]
    if fields:
        attr = _attr_dict(xref, dataset_name)
        overrides = COMPUTED_FIELDS.get(dataset_name, {})
        for f in fields:
            val = overrides[f](attr) if f in overrides else _flatten_scalar(attr.get(f))
            parts.append(_escape_pipe(val))
    return "|".join(parts)


def _search_row(xref: dict) -> str:
    """Encode an Xref as `id|dataset|name|xref_count` (matches GetSearchCompactRow)."""
    dataset_name = xref.get("dataset_name") or ""
    identifier = xref.get("identifier") or xref.get("keyword") or ""
    attr = _attr_dict(xref, dataset_name)
    name = _extract_name(attr, dataset_name)
    xref_count = str(xref.get("count") or 0)
    return "|".join(_escape_pipe(p) for p in (identifier, dataset_name, name)) + "|" + xref_count


def search_to_rest(resp: dict) -> dict:
    """gRPC SearchResponse -> REST search shape (fixed `id|dataset|name|xref_count`)."""
    results = ((resp.get("results") or {}).get("results")) or []
    seen, rows = set(), []
    for xref in results:
        identifier = xref.get("identifier") or xref.get("keyword") or ""
        dname = xref.get("dataset_name") or ""
        key = f"{identifier}:{dname}"
        if key in seen:
            continue
        seen.add(key)
        rows.append(_search_row(xref))
    return {
        "schema": "id|dataset|name|xref_count",
        "data": rows,
    }


def mapping_to_rest(resp: dict) -> dict:
    """gRPC MappingResponse -> REST map shape.

    REST emits:
        {schema: "id|<compact_fields>", mappings: [{source, targets:[strings]}]}

    Schema is determined by the FIRST observed target's dataset (uniform per
    map call by construction — the chain has a fixed terminal dataset)."""
    results = ((resp.get("results") or {}).get("results")) or []
    schema = "id"
    target_dataset = None
    encoded_mappings = []
    for m in results:
        src_xref = m.get("source") or {}
        # Source is a single id string in REST (just `id`), not pipe-encoded.
        src = (src_xref.get("identifier") or "")
        targets_out = []
        for t in (m.get("targets") or []):
            if target_dataset is None:
                target_dataset = t.get("dataset_name") or ""
                fields = COMPACT_FIELDS.get(target_dataset, [])
                schema = "|".join(["id"] + fields) if fields else "id"
            targets_out.append(_xref_compact_row(t))
        encoded_mappings.append({"source": src, "targets": targets_out})
    # gRPC returns the cursor at top-level `nextpage`; REST hangs it under
    # `pagination.next_token`. Normalize to REST shape so client.map_all's
    # while-loop reads pagination.next_token + pagination.has_next.
    inner = resp.get("results") or {}
    raw_next = inner.get("nextpage") or ""
    pagination = dict(inner.get("pagination") or {})
    if raw_next:
        pagination.setdefault("next_token", raw_next)
        pagination.setdefault("has_next", True)
    else:
        pagination.setdefault("has_next", bool(pagination.get("next_token")))

    if not encoded_mappings:
        schema = ""
    return {
        "schema": schema,
        "mappings": encoded_mappings if encoded_mappings else None,
        "pagination": pagination,
        "nextpage": raw_next,
    }


def entry_to_rest(resp: dict) -> dict:
    """gRPC EntryResponse -> REST entry shape.

    REST entry shape:
        {dataset, identifier, Attributes:{<DatasetPascal>:{...}},
         xrefs: {schema:"dataset|count", data:["dname|count", ...]}}

    PascalCase Attributes key matches the consumer pattern in
    atlas.disease.anchors._mondo_attrs (`Attributes.Mondo` / `Attributes.Ontology`).
    """
    result = resp.get("result") or {}
    dataset_name = result.get("dataset_name") or ""
    attr = _attr_dict(result, dataset_name)

    counts = result.get("dataset_counts") or []
    xref_rows = []
    for c in counts:
        ds_id = c.get("dataset")
        dname = (c.get("dataset_name")
                 or DATASET_ID_TO_NAME.get(int(ds_id) if ds_id is not None else -1, ""))
        if not dname:
            continue  # skip rows we can't resolve (REST would have a name from the registry)
        xref_rows.append(f"{dname}|{c.get('count', 0)}")

    # PascalCase key. For ontology-family datasets, REST emits `Ontology` (singular
    # PascalCase of the proto oneof). Disease anchors already handle both.
    pascal = _to_pascal(_ATTR_KEY_OVERRIDES.get(dataset_name, dataset_name))

    return {
        "dataset":    result.get("dataset_name"),
        "identifier": result.get("identifier"),
        "Attributes": {pascal: attr} if attr else {},
        "xrefs":      {"schema": "dataset|count", "data": xref_rows},
        "entries":    result.get("entries") or [],
    }


def _to_pascal(snake: str) -> str:
    """hgnc -> Hgnc; pharmgkb_clinical -> PharmgkbClinical."""
    return "".join(p[:1].upper() + p[1:] for p in snake.split("_"))
