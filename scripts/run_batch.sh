#!/usr/bin/env bash
# Run the parallel Atlas batch over corpus lists, optionally archiving the dist.
#
#   scripts/run_batch.sh [DIST] [WORKERS] [ARCHIVE_DIR]
#
# Defaults: DIST=/data/sugi-atlas-dist, WORKERS=(cores-2), ARCHIVE_DIR=unset.
# Corpus lists: env CORPUS_DIR (default corpus/) → {genes,diseases,drugs}.txt
#   (one entity per line, '#' comments allowed). For the full run, point
#   CORPUS_DIR at a dir of seed lists (e.g. cp corpus/seeds/genes_hgnc.txt ...).
# If ARCHIVE_DIR is given, the dist's atlas/ tree is tar.gz'd there on success
# (full corpus → archive, NOT git; see docs). The batch is resilient — a
# bad/unresolvable entity is skipped+logged, not fatal.
set -euo pipefail
cd "$(dirname "$0")/.."

DIST="${1:-/data/sugi-atlas-dist}"
WORKERS="${2:-$(( $(nproc) - 2 ))}"
ARCHIVE_DIR="${3:-}"
CORPUS_DIR="${CORPUS_DIR:-corpus}"
[ "$WORKERS" -lt 1 ] && WORKERS=1

echo "corpus($CORPUS_DIR): genes=$(grep -cvE '^\s*#|^\s*$' "$CORPUS_DIR/genes.txt") \
diseases=$(grep -cvE '^\s*#|^\s*$' "$CORPUS_DIR/diseases.txt") \
drugs=$(grep -cvE '^\s*#|^\s*$' "$CORPUS_DIR/drugs.txt") | workers=$WORKERS"

python -m atlas.batch \
  --dist "$DIST" \
  --workers "$WORKERS" \
  --genes    "@$CORPUS_DIR/genes.txt" \
  --diseases "@$CORPUS_DIR/diseases.txt" \
  --drugs    "@$CORPUS_DIR/drugs.txt"

# Archive the published corpus (deploy/cold-storage artifact — the full corpus
# is NOT committed to git; the pipeline + versioned inputs make it reproducible).
if [ -n "$ARCHIVE_DIR" ]; then
  mkdir -p "$ARCHIVE_DIR"
  STAMP="$(date +%Y%m%d-%H%M%S)"
  OUT="$ARCHIVE_DIR/atlas-corpus-$STAMP.tar.gz"
  echo "[archive] tarring $DIST/atlas → $OUT …"
  tar -C "$DIST" -czf "$OUT" atlas
  echo "[archive] $(du -h "$OUT" | cut -f1)  $OUT"
  ( cd "$DIST" && find atlas -mindepth 2 -maxdepth 2 -type d | sed 's#/[^/]*$##' \
      | sort | uniq -c | sed 's/^/[archive] pages: /' )
fi
