"""Compaction of schema.org JSON-LD for INLINE page embedding.

The complete graph is always written to the `entity.jsonld` sidecar; the
inline `<script type="application/ld+json">` block in the page body only needs
a representative slice. Large arrays (the full cohort-gene list, drug list,
HPO phenotypes, ontology codes) otherwise dominate the raw markdown — audit
gap #6: 173/268 disease pages spent 700–894 lines of JSON before the readable
lead, burying the one sentence an AI agent should extract first.

`compact_for_inline` caps every list-valued field (top-level, plus one level
into `@reverse`) to `INLINE_CAP` and stamps a `comment` pointing at the sidecar
for the complete set. The sidecar serializer is untouched, so no data is lost —
only the inline copy is trimmed.
"""
import copy

INLINE_CAP = 15


def compact_for_inline(jsonld: dict, cap: int = INLINE_CAP) -> dict:
    """Return a copy of `jsonld` with every over-long list truncated to `cap`.

    Truncates top-level list fields and one level into `@reverse` (the
    associatedGene / target edge arrays). Adds a `comment` noting the trim so a
    consumer knows to fetch `entity.jsonld` for the full graph. A graph with no
    over-long array is returned unchanged (no spurious comment)."""
    out = copy.deepcopy(jsonld)
    truncated = False

    def _cap(d: dict) -> None:
        nonlocal truncated
        for k, v in list(d.items()):
            if isinstance(v, list) and len(v) > cap:
                d[k] = v[:cap]
                truncated = True

    _cap(out)
    rev = out.get("@reverse")
    if isinstance(rev, dict):
        _cap(rev)

    if truncated:
        out["comment"] = (
            f"Inline graph truncated to {cap} items per list for readability; "
            "the complete schema.org graph is in the entity.jsonld sidecar.")
    return out
