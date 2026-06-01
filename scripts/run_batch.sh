#!/usr/bin/env bash
# Run the parallel Atlas batch over the curated corpus lists.
#
#   scripts/run_batch.sh [DIST] [WORKERS]
#
# Defaults: DIST=/data/sugi-atlas-dist, WORKERS=(cores-2).
# Corpus lists live in corpus/{genes,diseases,drugs}.txt (one entity per line,
# '#' comments allowed). Edit those to change the corpus. The batch is
# resilient — a bad/unresolvable entity is skipped+logged, not fatal.
set -euo pipefail
cd "$(dirname "$0")/.."

DIST="${1:-/data/sugi-atlas-dist}"
WORKERS="${2:-$(( $(nproc) - 2 ))}"
[ "$WORKERS" -lt 1 ] && WORKERS=1

echo "corpus: genes=$(grep -cvE '^\s*#|^\s*$' corpus/genes.txt) \
diseases=$(grep -cvE '^\s*#|^\s*$' corpus/diseases.txt) \
drugs=$(grep -cvE '^\s*#|^\s*$' corpus/drugs.txt)"

python -m atlas.batch \
  --dist "$DIST" \
  --workers "$WORKERS" \
  --genes    "@corpus/genes.txt" \
  --diseases "@corpus/diseases.txt" \
  --drugs    "@corpus/drugs.txt"
