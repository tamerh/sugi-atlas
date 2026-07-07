"""A FAQPage JSON-LD block for corpus pages (gene / drug / disease / pathway),
built from the page's title + declarative description + the already-computed
TL;DR key facts. Emitted as a SECOND inline <script> (additive — it does not
touch the entity JSON-LD), because a FAQPage/QAPage is the structured-data shape
AI assistants extract most reliably.

Two deliberately-robust, always-grammatical questions ("What is X?" and "What are
the key facts about X?") rather than fragile per-fact question generation.
"""
import json


def build_faq(title, description, tldr, url):
    """FAQPage dict, or None when there's nothing worth asking."""
    title = (title or "").strip()
    description = (description or "").strip()
    facts = [f.strip() for f in (tldr or []) if f and f.strip()]
    if not title or not (description or facts):
        return None
    qas = []
    if description:
        qas.append((f"What is {title}?", description))
    if facts:
        answer = f"Key facts about {title}: " + "; ".join(facts) + "."
        qas.append((f"What are the key facts about {title}?", answer))
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


def as_script_tag(title, description, tldr, url):
    """Inline <script> for the FAQPage, or '' when there's nothing to emit."""
    faq = build_faq(title, description, tldr, url)
    if not faq:
        return ""
    return f'<script type="application/ld+json">\n{json.dumps(faq, indent=2)}\n</script>'
