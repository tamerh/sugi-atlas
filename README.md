# Atlas

A deterministic reference atlas for **genes, drugs, and diseases**, built on
[biobtree](https://biobtree.org) and published at `sugi.bio/atlas/`. Code gathers
and renders every fact and statistic — no model decides what to query or how to
fill a table; a single optional, gated LLM step writes only the executive summary.

> **How it works, the per-entity build details, the frozen page contract, and the
> cross-entity mesh all live in [`docs/`](docs/)** — start with
> [docs/how-it-works.md](docs/how-it-works.md).

## Quick start

```bash
pip install -e .[dev]
# Needs the biobtree Go server running (default http://127.0.0.1:9291).
python -m atlas.gene.collect TP53 all      # collect + render one gene to stdout
```

`ATLAS_BIOBTREE` overrides the base URL (transport details in
[docs/how-it-works.md](docs/how-it-works.md) § biobtree transport).

## Build, test, release — `atlas.sh`

Everything builds into a local, gitignored `./dist`.

```bash
./atlas.sh test all               # dense build + unit + integration (the pre-prod gate)
./atlas.sh prod                   # full corpus, detached → dist/logs/ + tar.gz
./atlas.sh release X.Y.Z          # gate → annotated tag vX.Y.Z → push (pipeline only)
./atlas.sh release-publish vX.Y.Z # publish the GitHub release for a pushed tag (needs gh)
```

`prod` refuses to build unless the pre-prod gate is green and returns immediately
with a PID + log path. `release` versions the **pipeline only** (the corpus is
never attached); every build stamps `atlas_version` / `atlas_commit` +
`biobtree_version` into each page's frontmatter. `./atlas.sh help` for options.

ORA enrichment backgrounds (`data/background/`) are precomputed from biobtree and
refreshed per data release: `python -m atlas.build_background`.

## Layout

| Path | Role |
|---|---|
| `src/atlas/{gene,drug,disease}/` | per-entity collectors + renderers |
| `src/atlas/page/` | shared page assembly — frontmatter, mesh links, JSON-LD |
| `src/atlas/ora.py`, `batch.py` | enrichment statistics · parallel corpus builder (collect → merge → render) |
| `src/atlas/biobtree/` | biobtree REST client (direct `/ws/` + gate `/api/`) |
| `data/background/`, `data/corpus/` | ORA backgrounds · corpus seeds (gitignored) |
| `tests/{unit,integration}/` | logic guards · corpus-level contract + data-quality checks |
| `docs/` | how-it-works, per-entity build docs, page contract, mesh |
| `dist/` | **everything generated** (gitignored): pages, intermediates, cache, logs, tarballs |

## License

Code: **MIT** (see [LICENSE](LICENSE)). The generated Atlas **content** is published
separately under **CC BY-SA 4.0**; it derives from biobtree's integrated sources,
each carrying its own upstream license + attribution (full list on sugi.bio/atlas),
so downstream reuse must comply with those terms.
