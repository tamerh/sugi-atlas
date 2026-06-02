"""Coverage canaries — catch a whole collector silently producing nothing."""
from collections import defaultdict

import pytest

from ._harness import H2_IDS, report

pytestmark = pytest.mark.integration


def test_no_canonical_section_empty_corpus_wide(pages):
    """Every canonical H2 section has REAL content (a table or bold facts) on at
    least one page of its type. A section that is the bare placeholder on EVERY
    page means its collector/renderer is broken corpus-wide."""
    populated = defaultdict(set)     # entity -> {h2_id with real content somewhere}
    present = defaultdict(set)       # entity -> {h2_id that appears at all}
    for p in pages:
        for h2id, block in p.h2_blocks().items():
            present[p.entity].add(h2id)
            if "|" in block or "**" in block:   # a table row or a bold fact line
                populated[p.entity].add(h2id)
    bad = []
    for entity, ids in H2_IDS.items():
        for h2id in ids:
            if h2id in present[entity] and h2id not in populated[entity]:
                bad.append(f"{entity} #{h2id}: placeholder on ALL pages")
    assert not bad, report(bad)
