# Atlas

A comprehensive bio atlas for genes, drugs, and diseases — built on top of
[biobtree](https://biobtree.org) with a **deterministic data collector** + a
narrow **LLM synthesis** step. Outputs are published as a static atlas at
`sugi.bio/atlas/`.

```
biobtree (local REST API)  ──►  collect   ──►  render   ──►  page tables
  127.0.0.1:8000/api            (no model)     (no model)
                                    │
                                    └──►  body  ──►  LLM summary  ──►  exec summary
                                                     (the only model step)
```

## Layout

| Path | Role |
|---|---|
| `src/atlas/gene/` | gene collectors + renderer (12 sections) |
| `src/atlas/drug/`, `src/atlas/disease/` | future entity types |
| `src/atlas/validation/` | `coverage.py` (25-gene fact-coverage), `stress.py` (diversity stress), `summary_gate.py` (hardened LLM judge) |
| `src/atlas/bench/` | summary-model benchmark + dataset-coverage + mining tools |
| `src/atlas/cli.py` | thin CLI entry point used by Enju workflows |
| `validation_data/gene/` | 25-gene reference set for fact-coverage |
| `workflows/` | Enju workflow YAMLs (genes / drugs / diseases) |
| `examples/` | sample generated pages (committed for OSS reviewers) |
| `docs/BIOBTREE_ISSUES.md` | upstream biobtree API issues / improvement requests |
| `docs/COLLECTOR_NOTES.md` | findings + design notes from the collector prototype |

## Install & run (dev)

```bash
pip install -e .
ATLAS_BIOBTREE=http://127.0.0.1:8000 python -m atlas.gene.collect TP53 all
python -m atlas.validation.coverage 7
python -m atlas.validation.stress
```

## Publishing

Atlas is the generator. Final pages are written to `$ATLAS_DIST/atlas/{gene,drug,disease}/...`
and committed to a private `sugi-atlas-dist` repo, which is mounted by the
`biobtree-content` Hugo site at build time and served at `sugi.bio/atlas/`.

## History

The collector was prototyped in `/data/biobtree-collector` on branch
`deterministic-collector` (a worktree of the `biobtree-content` Hugo
repository). See `docs/COLLECTOR_NOTES.md` for the rationale, design findings,
and per-section validation summary inherited from that work.
