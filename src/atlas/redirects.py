"""Slug history → redirect map.

Page slugs derive from names (disease: Mondo label, gene: HGNC symbol, drug:
ChEMBL name), so a rename changes the URL — v1.11.10 moved 34 disease pages.
`data/corpus/slug_history.json` records every slug a stable identifier
(MONDO / HGNC / CHEMBL id) has ever been published under. After a build, any
historical slug whose identifier was rebuilt under a DIFFERENT slug becomes a
redirect old → new, written to `atlas/redirects.json` for the web layer (301s).

Only identifiers built in this run are resolved, so a sample build never claims
that the rest of the corpus vanished; retired identifiers (no page any more)
are left to a manual map.
"""
import json
import os

HISTORY = os.path.join(os.path.dirname(__file__), "..", "..", "data", "corpus", "slug_history.json")


def load_history(path=HISTORY):
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def compute(history, built):
    """history: {entity: {slug: ident}}; built: [(entity, slug, ident)] from this
    run. Returns ({entity: {old_slug: new_slug}}, updated history)."""
    current = {}                                  # (entity, ident) -> slug
    for et, slug, ident in built:
        if ident:
            current[(et, ident)] = slug
    redirects = {}
    for et, slugs in history.items():
        for old, ident in slugs.items():
            new = current.get((et, ident))
            if new and new != old:
                redirects.setdefault(et, {})[old] = new
    updated = {et: dict(v) for et, v in history.items()}
    for (et, ident), slug in current.items():
        updated.setdefault(et, {})[slug] = ident
    # A slug now published for a different identifier is no longer a redirect
    # source (the page there is real).
    live = {(et, s) for (et, _i), s in current.items()}
    redirects = {et: {o: n for o, n in m.items() if (et, o) not in live}
                 for et, m in redirects.items()}
    return {et: dict(sorted(m.items())) for et, m in redirects.items() if m}, updated


# A run this large is the full corpus — only then is "not built" evidence that a
# page was retired (a sample build would otherwise retire everything).
FULL_BUILD_MIN = 40000


def retired(history, built):
    """{entity: [old slugs]} whose identifier has no page in this (full) build —
    for the web layer to answer 410 Gone instead of 404."""
    live_ids = {(et, i) for et, _s, i in built}
    live_slugs = {(et, s) for et, s, _i in built}
    out = {}
    for et, slugs in history.items():
        for slug, ident in slugs.items():
            if (et, ident) not in live_ids and (et, slug) not in live_slugs:
                out.setdefault(et, []).append(slug)
    return {et: sorted(v) for et, v in out.items()}


def write(dist_dir, built, history_path=HISTORY):
    """Compute + write atlas/redirects.json, and the history updated with this
    run's slugs to <dist>/slug_history.json. The repo copy is never touched by a
    build (a tracked file changing under a sample/prod run blocks the release
    tag); after a full prod build, copy dist/slug_history.json over
    data/corpus/slug_history.json and commit it. Returns the redirect map."""
    history = load_history(history_path)
    redirects, updated = compute(history, built)
    with open(os.path.join(dist_dir, "atlas", "redirects.json"), "w") as f:
        json.dump(redirects, f, indent=1, sort_keys=True)
    if len(built) >= FULL_BUILD_MIN:
        with open(os.path.join(dist_dir, "atlas", "retired.json"), "w") as f:
            json.dump(retired(history, built), f, indent=1, sort_keys=True)
    with open(os.path.join(dist_dir, "slug_history.json"), "w") as f:
        json.dump(updated, f, indent=0, sort_keys=True)
    return redirects
