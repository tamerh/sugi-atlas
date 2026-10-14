"""A FAQPage JSON-LD block for corpus pages (gene / drug / disease / pathway),
built from the page's title + declarative description + the already-computed
TL;DR key facts. Emitted as a SECOND inline <script> (additive — it does not
touch the entity JSON-LD), because a FAQPage/QAPage is the structured-data shape
AI assistants extract most reliably.

Two always-present questions ("What is X?" and "What are the key facts about
X?") plus the intent-specific ones from atlas.page.faq_intents, each emitted only
when the page's data answers it.
"""
import json


def build_faq(title, description, tldr, url, extra=None):
    """FAQPage dict, or None when there's nothing worth asking."""
    title = (title or "").strip()
    description = (description or "").strip()
    # Strip each fact's own trailing period — they're joined with "; " and the
    # answer adds one (v1.11.7: "…no protein product..").
    facts = [f.strip().rstrip(".").rstrip() for f in (tldr or []) if f and f.strip()]
    facts = [f for f in facts if f]
    extra = [(q, a) for q, a in (extra or []) if q and a]
    if not title or not (description or facts or extra):
        return None
    qas = []
    if description:
        qas.append((f"What is {title}?", description))
    if facts:
        answer = f"Key facts about {title}: " + "; ".join(facts) + "."
        qas.append((f"What are the key facts about {title}?", answer))
    qas.extend(extra)
    if not qas:
        return None
    return {
        "@context": "https://schema.org",
        "@type": "FAQPage",
        "@id": (url or "") + "#faq",
        "mainEntity": [
            {"@type": "Question", "name": q,
             "acceptedAnswer": {"@type": "Answer", "text": a}}
            for q, a in qas],
    }


def as_script_tag(title, description, tldr, url, extra=None):
    """Inline <script> for the FAQPage, or '' when there's nothing to emit."""
    faq = build_faq(title, description, tldr, url, extra)
    if not faq:
        return ""
    # Minified inline (the pretty form was ~18% of page bytes across JSON-LD blocks).
    return f'<script type="application/ld+json">\n{json.dumps(faq, ensure_ascii=False, separators=(",", ":"))}\n</script>'
