#!/usr/bin/env python3
"""Corpus statistics for the preprint — reproducible figures from a built corpus.

Reads the *published* page frontmatter (the shipped artifact, not an intermediate
cache) plus the indication sidecar, and emits the paper's corpus tables as
markdown so every reported number is re-derivable:

    python -m atlas.stats <dist>          # dist = a built corpus dir (contains atlas/)
    python -m atlas.stats <archive.tar.gz>

Tables produced:
  1. Corpus scale — pages per entity type.
  2. Section-level coverage — %% of pages of a type carrying real data for each
     evidence-bearing signal (component count > 0).
  3. Disease content depth — the "is it empty?" distribution.
  4. Indication coverage — diseases with a disease-direct indicated drug.

Disease-depth classification (documented here so the figure is auditable):
  rich      gene_count > 0                         an associated-gene cohort resolves,
                                                   so the molecular sections populate.
  clinical  gene_count == 0 and any of             characterized clinically/genetically
            {gwas,variant,trial,civic,             but with no cohort (e.g. antibody-
             phenotype}_count > 0                  mediated / autoimmune, anti-NMDA-like).
  thin      all of the above == 0                  identifiers + ontology family only.

Note: phenotype_count is emitted into the frontmatter only from builds at or
after the commit that added it to the disease evidence _SPEC; on older corpora
it is absent, so symptom-only diseases fall into "thin" there.
"""
import glob
import json
import os
import re
import sys
import tarfile
from collections import Counter, defaultdict

_ET_RE = re.compile(r'^entity_type:\s*"?(\w+)"?', re.M)
_COMP_BLOCK = re.compile(r'^evidence_components:\n((?:[ ]+\w+:.*\n)+)', re.M)

# Disease-depth signals: a non-cohort disease counts as "clinical/genetic" (not
# "thin") if it carries any of these.
_CLINICAL_SIGNALS = ("gwas_count", "variant_count", "trial_count", "civic_count",
                     "phenotype_count")


def _frontmatter(text):
    """The YAML block between the first two '---' fences ('' if none)."""
    parts = text.split("---", 2)
    return parts[1] if len(parts) >= 3 else ""


def _components(front):
    """The evidence_components int map from a page's frontmatter."""
    out = {}
    m = _COMP_BLOCK.search(front)
    if m:
        for line in m.group(1).splitlines():
            k, _, v = line.strip().partition(":")
            try:
                out[k] = int(v.strip())
            except ValueError:
                pass
    return out


def _iter_pages(src):
    """Yield (entity_type, components) for every page in a dist dir or .tar.gz."""
    if src.endswith((".tar.gz", ".tgz")):
        with tarfile.open(src) as t:
            for m in t:
                if m.name.endswith("/index.md") and m.isfile():
                    text = t.extractfile(m).read().decode("utf-8", "replace")
                    fm = _frontmatter(text)
                    et = _ET_RE.search(fm)
                    if et:
                        yield et.group(1), _components(fm)
    else:
        root = src if os.path.basename(src.rstrip("/")) == "atlas" else os.path.join(src, "atlas")
        for p in glob.glob(os.path.join(root, "*", "*", "index.md")):
            with open(p, encoding="utf-8") as f:
                fm = _frontmatter(f.read())
            et = _ET_RE.search(fm)
            if et:
                yield et.group(1), _components(fm)


def _disease_depth(c):
    if (c.get("gene_count") or 0) > 0:
        return "rich"
    if any((c.get(k) or 0) > 0 for k in _CLINICAL_SIGNALS):
        return "clinical"
    return "thin"


def _indication_count(src):
    """Number of diseases with >=1 disease-direct indicated/candidate drug."""
    try:
        if src.endswith((".tar.gz", ".tgz")):
            with tarfile.open(src) as t:
                for m in t:
                    if m.name.endswith("indicated_drugs.json"):
                        return len(json.load(t.extractfile(m)))
            return 0
        root = src if os.path.basename(src.rstrip("/")) == "atlas" else os.path.join(src, "atlas")
        with open(os.path.join(root, "indicated_drugs.json")) as f:
            return len(json.load(f))
    except (OSError, json.JSONDecodeError, KeyError):
        return -1   # sidecar absent


def _pct(n, total):
    return round(100 * n / total) if total else 0


def compute(src):
    scale = Counter()
    coverage = defaultdict(Counter)          # entity_type -> component -> #(>0)
    depth = Counter()
    for et, comp in _iter_pages(src):
        scale[et] += 1
        for k, v in comp.items():
            if v > 0:
                coverage[et][k] += 1
        if et == "disease":
            depth[_disease_depth(comp)] += 1
    return scale, coverage, depth, _indication_count(src)


# Pretty labels for the coverage rows (component key -> human label).
_LABELS = {
    "civic_count": "CIViC evidence", "drug_count": "Drugs", "gene_count": "Cohort genes",
    "variant_count": "Variants (ClinVar)", "gwas_count": "GWAS", "trial_count": "Clinical trials",
    "interaction_count": "Protein interactions",
}


def render_markdown(src, scale, coverage, depth, ind):
    total = sum(scale.values())
    L = [f"# Corpus statistics — {os.path.basename(src.rstrip('/'))}", ""]

    L += ["## Table 1 — Corpus scale", "",
          "| Entity type | Pages | Share |", "|---|---:|---:|"]
    for et in ("gene", "disease", "drug"):
        L.append(f"| {et.capitalize()}s | {scale[et]:,} | {_pct(scale[et], total)}% |")
    L.append(f"| **Total** | **{total:,}** | — |")

    L += ["", "## Table 2 — Section-level data coverage",
          "", "%% of pages of a type carrying real data for each signal.", ""]
    for et in ("gene", "disease", "drug"):
        n = scale[et]
        if not n:
            continue
        rows = sorted(coverage[et].items(), key=lambda kv: -kv[1])
        L += [f"### {et.capitalize()}s (n = {n:,})", "", "| Signal | % |", "|---|---:|"]
        for k, c in rows:
            L.append(f"| {_LABELS.get(k, k)} | {_pct(c, n)} |")
        L.append("")

    dn = sum(depth.values())
    L += ["## Table 3 — Disease content depth", "",
          f"Of {dn:,} disease pages, by how much content resolves:", "",
          "| Disease page | Count | Share |", "|---|---:|---:|",
          f"| **Rich** — associated-gene cohort, full molecular sections | {depth['rich']:,} | {_pct(depth['rich'], dn)}% |",
          f"| **Clinical / genetic only** — GWAS/variants/trials, no cohort (anti-NMDA-like) | {depth['clinical']:,} | {_pct(depth['clinical'], dn)}% |",
          f"| **Thin** — identifiers + ontology family only | {depth['thin']:,} | {_pct(depth['thin'], dn)}% |",
          ""]
    if ind >= 0:
        L += [f"**Disease-direct drug coverage:** {ind:,} diseases "
              f"({_pct(ind, dn)}%) carry at least one indicated or late-stage "
              f"candidate drug (independent of the gene cohort).", ""]
    return "\n".join(L)


def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]
    if not argv:
        print("usage: python -m atlas.stats <dist-dir|archive.tar.gz>", file=sys.stderr)
        return 2
    src = argv[0]
    scale, coverage, depth, ind = compute(src)
    print(render_markdown(src, scale, coverage, depth, ind))
    return 0


if __name__ == "__main__":
    sys.exit(main())
