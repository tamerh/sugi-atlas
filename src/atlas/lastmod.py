"""Per-page `lastmod` from a content hash.

Every page's `generated_at` is the build time, so sitemaps carried no change
signal (all 52k pages "changed" every release). After render, each page's
content — minus build stamps — is hashed and compared with the registry from the
previous full build (`data/corpus/content_dates.json`): unchanged pages keep
their previous date, changed or new pages get today's. The date is written to
the frontmatter as `lastmod:`; the updated registry goes to
`<dist>/content_dates.json` (copy it over the repo file after a prod build, like
the slug history — builds never write the repo copy).
"""
import hashlib
import json
import os
import re
from datetime import datetime, timezone

REGISTRY = os.path.join(os.path.dirname(__file__), "..", "..", "data", "corpus", "content_dates.json")

# Build stamps that change every run without the content changing.
_VOLATILE = re.compile(
    r'^(generated_at|atlas_version|atlas_commit|biobtree_version|biobtree_commit|lastmod): .*$',
    re.M)
_DATES = re.compile(r"\b20\d\d-\d\d-\d\d(?:T[\d:.+Z-]+)?\b")


def content_hash(md):
    s = _VOLATILE.sub("", md)
    s = _DATES.sub("", s)
    return hashlib.sha1(s.encode("utf-8")).hexdigest()[:16]


def load(path=REGISTRY):
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def stamp(md, key, prev, today):
    """(new markdown with `lastmod:` set, registry entry)."""
    h = content_hash(md)
    old = prev.get(key)
    date = old[1] if old and old[0] == h else today
    md = re.sub(r'^lastmod: .*\n', "", md, count=1, flags=re.M)
    md = re.sub(r'^(generated_at: .*\n)', rf'\1lastmod: "{date}"\n', md, count=1, flags=re.M)
    return md, [h, date]


def apply(dist_dir, pages, registry_path=REGISTRY, today=None):
    """Stamp every rendered page [(entity, slug)]; write <dist>/content_dates.json.
    Returns (changed, unchanged) counts."""
    today = today or datetime.now(timezone.utc).date().isoformat()
    prev = load(registry_path)
    out = {}
    for et, slug in pages:
        p = os.path.join(dist_dir, "atlas", et, slug, "index.md")
        try:
            with open(p, encoding="utf-8") as f:
                md = f.read()
        except OSError:
            continue
        key = f"{et}/{slug}"
        new, entry = stamp(md, key, prev, today)
        out[key] = entry
        if new != md:
            with open(p, "w", encoding="utf-8") as f:
                f.write(new)
    with open(os.path.join(dist_dir, "content_dates.json"), "w") as f:
        json.dump(out, f, indent=0, sort_keys=True)
    changed = sum(1 for e in out.values() if e[1] == today)
    return changed, len(out) - changed
