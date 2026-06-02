"""Corpus-level integration checks over a BUILT dist (the dense set).

The unit suite guards the logic; this suite guards the actual rendered pages —
a release gate to run after a dense build, before committing to the full corpus.
Every check here corresponds to an invariant we rely on (the frozen page
contract, frontmatter schema, data-quality guards, mesh-link integrity,
JSON-LD validity), so a regression in any of them fails loudly against real
output.

Point it at a dist with ATLAS_INTEGRATION_DIST (default /data/sugi-atlas-dist).
The whole suite skips cleanly if no dist is present, so a plain `pytest` on a
machine without a build still runs the unit tests and skips these.

    pytest -m integration                      # corpus checks (needs a dist)
    pytest -m "not integration"                # unit only
    ATLAS_INTEGRATION_DIST=/tmp/x pytest -m integration
"""
import glob
import json
import os
import re

import pytest
import yaml

DIST = os.environ.get("ATLAS_INTEGRATION_DIST", "/data/sugi-atlas-dist")
ATLAS = os.path.join(DIST, "atlas")

# The FROZEN page contract (docs/PAGE_CONTRACT.md) — anchor ids in order, per
# entity. Hardcoded here independently of the generator, so a code change that
# drifts the H2 set/order fails this gate.
H2_IDS = {
    "gene":    ["summary", "identifiers", "gene-structure", "protein",
                "function", "disease", "drugs", "related"],
    "disease": ["summary", "identifiers", "genetics", "genes", "function",
                "drugs", "trials", "related"],
    "drug":    ["summary", "identifiers", "targets", "indications",
                "pharmacology", "related-molecules", "related"],
}
ID_RE = {
    "gene":    re.compile(r"^[A-Za-z0-9][A-Za-z0-9._\-]*$"),   # HGNC symbol
    "disease": re.compile(r"^MONDO:\d+$"),
    "drug":    re.compile(r"^CHEMBL\d+$"),
}

_H2_LINE = re.compile(r"^## (.+?)(?:\s*\{#([a-z0-9-]+)\})?\s*$")
_ANY_HEADING = re.compile(r"^(#{2,6}) (.+?)(?:\s*\{#([a-z0-9-]+)\})?\s*$")
_INTERNAL_LINK = re.compile(r"\]\((/atlas/(gene|disease|drug)/([^/)#]+)/[^)]*)\)")


def _maybe_skip():
    if not os.path.isdir(ATLAS) or not glob.glob(os.path.join(ATLAS, "*", "*", "page.md")):
        pytest.skip(f"no built dist at {ATLAS} (set ATLAS_INTEGRATION_DIST)",
                    allow_module_level=True)


class Page:
    def __init__(self, entity, slug, path):
        self.entity, self.slug, self.path = entity, slug, path
        self.raw = open(path).read()
        parts = self.raw.split("---\n", 2)
        if len(parts) >= 3:
            try:
                self.fm = yaml.safe_load(parts[1]) or {}
            except yaml.YAMLError:
                self.fm = {"__yaml_error__": True}
            self.body = parts[2]
        else:
            self.fm, self.body = {}, self.raw
        jp = os.path.join(os.path.dirname(path), "entity.jsonld")
        self.jsonld_raw = open(jp).read() if os.path.exists(jp) else None

    @property
    def h2(self):
        """[(label, id_or_None)] for every '## ' heading, in source order."""
        return [(m.group(1).strip(), m.group(2))
                for line in self.body.splitlines()
                if (m := _H2_LINE.match(line))]

    def jsonld(self):
        return json.loads(self.jsonld_raw) if self.jsonld_raw else None


def _load():
    pages = []
    for path in glob.glob(os.path.join(ATLAS, "*", "*", "page.md")):
        entity = os.path.basename(os.path.dirname(os.path.dirname(path)))
        slug = os.path.basename(os.path.dirname(path))
        pages.append(Page(entity, slug, path))
    return pages


def report(violations, limit=25):
    """A readable assertion message: count + the first `limit` offenders."""
    n = len(violations)
    head = "\n".join(f"  - {v}" for v in violations[:limit])
    more = f"\n  … +{n - limit} more" if n > limit else ""
    return f"{n} violation(s):\n{head}{more}"
