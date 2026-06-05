# Atlas

A comprehensive bio atlas for **genes, drugs, and diseases** — built on top of
[biobtree](https://biobtree.org) with a **deterministic data collector** + a
narrow **LLM synthesis** step. Pages are rendered as static markdown and served
at `sugi.bio/atlas/`.

```
biobtree (local REST API)  ──►  collect   ──►  render   ──►  page tables
  127.0.0.1:9291/ws (direct)    (no model)     (no model)
                                    │
                                    └──►  body  ──►  LLM summary  ──►  exec summary
                                                     (the only model step)
```

## Layout

| Path | Role |
|---|---|
| `src/atlas/gene/`, `drug/`, `disease/` | per-entity collectors + renderers |
| `src/atlas/biobtree/` | biobtree REST client (direct `/ws/` + gate `/api/` transports) |
| `src/atlas/page/` | shared page assembly — frontmatter, mesh links, JSON-LD |
| `src/atlas/batch.py` | parallel corpus builder (collect → merge → render) |
| `src/atlas/validation/` | `coverage.py` (25-gene fact-coverage), `stress.py`, `summary_gate.py` |
| `src/atlas/bench/` | summary-model benchmark + dataset-coverage + mining tools |
| `src/atlas/*/enju.yaml` | Enju workflow definitions (gene / drug / disease) |
| `data/corpus/dense/` | the dense reference set (~1k entities) — the test/release gate |
| `data/corpus/seeds/` | full-corpus seed lists (gitignored; regenerable) |
| `data/validation/gene/` | 25-gene reference set for fact-coverage |
| `tests/unit/`, `tests/integration/` | unit suite + corpus-level integration checks |
| `docs/` | how-it-works + per-entity build docs (gene/drug/disease) + the mesh + page contract |
| `dist/` | **everything generated** (gitignored): `atlas/` pages, `build/` corpus intermediates, `cache/` bundle cache, `logs/`, release tarballs |

## Install & run (dev)

```bash
pip install -e .[dev]
# Defaults to the biobtree Go server direct (http://127.0.0.1:9291, /ws/ + mode=lite)
# — ~8x faster than the FastAPI gate for bulk page-gen, identical response shapes.
python -m atlas.gene.collect TP53 all
python -m atlas.validation.coverage 7
```

`ATLAS_BIOBTREE` overrides the base URL. To route through the FastAPI/MCP gate
instead (e.g. the only endpoint exposed in some environments):

```bash
ATLAS_BIOBTREE=http://127.0.0.1:8000 ATLAS_BIOBTREE_TRANSPORT=gate \
  python -m atlas.gene.collect TP53 all
```

## Build & test — `atlas.sh`

The orchestrator wraps the everyday loops. Everything builds into a local,
gitignored `./dist`.

```bash
./atlas.sh test          # unit suite — fast, no build
./atlas.sh test int      # integration suite over ./dist (no rebuild)
./atlas.sh test all      # dense build + unit + integration (the pre-prod gate)
./atlas.sh prod          # full corpus, detached (nohup → dist/logs/): pre-prod check
                         #   → corpus from seeds → integration sweep → tar.gz
./atlas.sh release X.Y.Z         # gate (test all) → annotated tag vX.Y.Z → push
./atlas.sh release-publish vX.Y.Z # publish the GitHub release for a pushed tag (needs gh)
```

`prod` refuses to build the corpus unless the pre-production check is green, and
returns immediately with a PID + log path (`dist/logs/prod-<stamp>.log`).
`release` versions the **pipeline only** (the corpus is never attached); each
build stamps `atlas_version`/`atlas_commit` + the `biobtree_version` into every
page's frontmatter. Run `./atlas.sh help` for options (`--dist`, `--workers`, `--limit`).

## Two-stage testing

- **Unit** (`tests/unit/`) guards the collector/renderer logic.
- **Integration** (`tests/integration/`) guards the actual rendered pages against
  the frozen page contract (`docs/PAGE_CONTRACT.md`): H2/H3 anchors, frontmatter
  schema, data-quality, mesh-link integrity, JSON-LD. It runs over a built
  `./dist` and skips cleanly if none is present.

## History

The collector was prototyped in `/data/biobtree-collector` on branch
`deterministic-collector` (a worktree of the `biobtree-content` Hugo
repository): mining the published pages showed the tool-call trajectory for
these structured sections is mostly deterministic, so data-gathering and
table-rendering use no model at all — only the executive summary does.

## License

The code in this repository is licensed **MIT** (see [LICENSE](LICENSE)).

The generated Atlas **content** is published separately under **CC BY-SA 4.0**.
It is derived from biobtree's integrated sources, each of which carries its own
license and attribution requirements; the full per-source list is provided on
the Atlas site (sugi.bio/atlas). Downstream reuse of the content must comply
with those upstream terms.
