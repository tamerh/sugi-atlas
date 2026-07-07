"""Pilot build for /atlas/variant/ pages — enumerate a gene's ClinVar variants
(P/LP + conflicting gate), render a page each + a per-gene index, into a dist
tree that slots alongside the main corpus.

Isolated from the main build on purpose (no _ENTITY_TYPES surgery yet): it loads
the EXISTING corpus manifest so variant→gene / variant→disease links gate
correctly against pages that already exist, and emits variant URLs directly.

    python -m atlas.variant.build --dist dist --out dist \
        --genes ACTA1,DACT1,ASXL1,NAA10,RPL10
"""
import argparse
import json
import os

from atlas.biobtree import search, rows
from atlas.pipeline import (build_meta, biobtree_version, atlas_version, _yaml_escape)
from atlas.page import links
from atlas.variant import collect as VC, render as VR, enrich as EN

_DATASETS = ("clinvar", "clingen_variant", "clingen_gene_validity", "clingen_dosage",
             "gencc", "mondo", "orphanet", "gard", "clinical_trials", "panelapp_gene",
             "alphamissense", "dbsnp", "gnomad_constraint", "uniprot", "pdb", "alphafold",
             "hgnc")
_CLASS_ORDER = ["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic",
                "Conflicting classifications of pathogenicity"]


def _hgnc_id(symbol):
    for r in rows(search(symbol, source="hgnc")):
        if r.get("id", "").startswith("HGNC:"):
            return r["id"]
    return None


def _frontmatter(meta, description, identifier, aliases, entity_type="variant"):
    fm = ["---"]
    ordered = [("title", meta.get("title")), ("identifier", identifier),
               ("entity_type", entity_type), ("gene", meta.get("gene")),
               ("classification", meta.get("classification")),
               ("description", description)]
    for k, v in ordered:
        if v:
            fm.append(f'{k}: "{_yaml_escape(v)}"')
    if aliases:                       # Hugo-reserved → 301 redirects (c.-form → page)
        fm.append("aliases:")
        fm += [f'  - "{a}"' for a in aliases]
    for k in ("generated_at", "atlas_version", "atlas_commit",
              "biobtree_version", "biobtree_commit", "generated_by"):
        if meta.get(k):
            fm.append(f'{k}: "{_yaml_escape(str(meta[k]))}"')
    fm.append(f"datasets: [{', '.join(_DATASETS)}]")
    fm.append("---\n")
    return "\n".join(fm)


def _write(out_root, rel_slug, page_md):
    d = os.path.join(out_root, "atlas", "variant", rel_slug)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "index.md"), "w") as f:
        f.write(page_md)


def build_gene(symbol, out_root):
    """Build every buildable variant page for one gene + its index. Returns the
    list of built variant records (for the index + summary)."""
    hgnc = _hgnc_id(symbol)
    if not hgnc:
        print(f"  {symbol}: no HGNC id — skip")
        return []
    # Pass 1 — collect every buildable variant (pure ClinVar record).
    recs, seen_slug = [], set()
    for vid, cls, name in VC.enumerate_gene(hgnc):
        rec = VC.collect(vid)
        if not rec or rec["canonical_slug"] in seen_slug:
            continue
        seen_slug.add(rec["canonical_slug"])
        recs.append(rec)
    # Per-gene enrichment context, fetched once and shared across the gene's
    # variants (AlphaMissense map, residue index, gene ACMG context, structure,
    # the primary-condition digest, diagnostic panels, and the P/LP scale count).
    pl_count = sum(1 for r in recs if VC._is_pathogenic(r.get("classification")))
    ctx = {
        "am": EN.gene_alphamissense(hgnc),
        "positions": VC.build_position_index(recs),
        "recs": recs,
        "gene_context": EN.gene_context(hgnc),
        "structure": EN.gene_structure(hgnc),
        "digest": EN.gene_condition_digest(hgnc),
        "panels": EN.gene_panelapp(hgnc),
        "pl_count": pl_count,
    }
    # Pass 2 — enrich + render each.
    for rec in recs:
        VC.attach_enrichment(rec, ctx)
        meta = build_meta("variant", rec["canonical_slug"], VR._label(rec), _DATASETS)
        meta["gene"] = rec["gene_symbol"]
        meta["classification"] = rec["classification"]
        aliases = [f"/atlas/variant/{s}/" for s in rec["slugs"][1:]]   # non-canonical → redirects
        page = (_frontmatter(meta, VR.declarative_plain(rec), f"VCV{rec['variation_id']}", aliases)
                + VR.render_body(rec))
        _write(out_root, rec["canonical_slug"], page)
    built = recs
    _write_index(symbol, hgnc, built, out_root)
    print(f"  {symbol}: {len(built)} variant pages + index")
    return built


def _write_index(symbol, hgnc, recs, out_root):
    """Per-gene variant index (crawlable hub) at /atlas/variant/gene/<sym>/."""
    by_cls = {}
    for r in recs:
        by_cls.setdefault(r["classification"], []).append(r)
    title = f"{symbol} — clinically significant variants"
    L = [f"## {symbol} variants (ClinVar)", "",
         f"**{len(recs)} pathogenic / likely-pathogenic / conflicting variants** in "
         f"{links.maybe_link(symbol, links.gene_url(symbol=symbol, hgnc_id=hgnc)) or symbol} "
         "with a dedicated page. Benign and unreviewed variants are omitted.", ""]
    for cls in _CLASS_ORDER + [c for c in by_cls if c not in _CLASS_ORDER]:
        group = sorted(by_cls.get(cls, []), key=lambda r: r["canonical_slug"])
        if not group:
            continue
        L += [f"### {cls} ({len(group)}) {{#{cls.split()[0].lower()}}}", ""]
        L.append(", ".join(f"[{VR._label(r)}](/atlas/variant/{r['canonical_slug']}/)"
                           for r in group))
        L.append("")
    meta = build_meta("variant", f"gene/{symbol.lower()}", title, _DATASETS)
    meta["gene"] = symbol
    desc = (f"{len(recs)} clinically significant ClinVar variants in {symbol} "
            "(pathogenic, likely-pathogenic, conflicting) — each with a dedicated "
            "reference page.")
    # index is a listing page (its own type → plain title, but still sitemapped)
    page = _frontmatter(meta, desc, f"{symbol}-variants", [], entity_type="variant-index") + "\n".join(L)
    _write(out_root, f"{symbol.lower()}-variants", page)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True, help="comma-separated symbols")
    ap.add_argument("--dist", default="dist", help="corpus dist (for the mesh manifest)")
    ap.add_argument("--out", default="dist", help="where to write variant pages")
    a = ap.parse_args(argv)
    links.load(a.dist)     # gate variant→gene/disease links on the real corpus
    print(f"variant pilot | biobtree {biobtree_version()} | atlas {atlas_version()}")
    total = 0
    for sym in [s.strip() for s in a.genes.split(",") if s.strip()]:
        total += len(build_gene(sym, a.out))
    print(f"done — {total} variant pages under {a.out}/atlas/variant/")


if __name__ == "__main__":
    main()
