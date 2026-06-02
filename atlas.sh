#!/usr/bin/env bash
#
# atlas.sh — Sugi Atlas build & test orchestrator.
#
# Two verbs cover the whole loop. `test` is the inner loop (validate the
# generator); `prod` is the outer loop (ship the real corpus). Both build into
# a local, gitignored ./dist — nothing is pushed to a separate dist repo.
#
#   ./atlas.sh test [int|all]
#   ./atlas.sh prod
#   ./atlas.sh help
#
#   test          Unit suite only — fast, no build. The default.
#   test int      Integration suite over ./dist (must already be built).
#   test all      Pre-production check: build the dense set into ./dist, then
#                 run unit + integration against it. Run this before `prod`.
#
#   prod          Production build, detached via nohup (logs → logs/):
#                   1. pre-production check  (== test all: dense build + tests)
#                   2. full corpus from corpus/seeds/ into ./dist
#                   3. integration sweep over the full output
#                   4. versioned tar.gz archive (dist/atlas-corpus-<stamp>-<sha>)
#                 Refuses to build the corpus if any test fails. Returns
#                 immediately with a PID + log path; tail the log to follow.
#
# Options (also read from the environment):
#   -d, --dist DIR    output dist     (default $ATLAS_DIST or ./dist)
#   -w, --workers N   parallel workers (default: nproc-2; env ATLAS_WORKERS)
#   -n, --limit N     build only first N of each type — quick smoke
#   -h, --help        show help
#
# Examples:
#   ./atlas.sh test            # fast unit suite
#   ./atlas.sh test all        # dense build + unit + integration
#   ./atlas.sh test all -n 5   # same, 5-of-each smoke
#   ./atlas.sh prod            # full corpus, detached, logged, archived
#
set -euo pipefail
cd "$(dirname "$0")"

# ---- defaults ---------------------------------------------------------------
DIST="${ATLAS_DIST:-$(pwd)/dist}"
WORKERS="${ATLAS_WORKERS:-}"
LIMIT="${ATLAS_LIMIT:-}"
BIOBTREE="${ATLAS_BIOBTREE:-http://127.0.0.1:9291}"
DENSE_DIR="corpus/dense"
SEED_DIR="corpus/seeds"

# ---- pretty -----------------------------------------------------------------
if [ -t 1 ]; then B=$'\e[1m'; G=$'\e[32m'; Y=$'\e[33m'; R=$'\e[31m'; D=$'\e[2m'; N=$'\e[0m'
else B=""; G=""; Y=""; R=""; D=""; N=""; fi
say()  { printf "%s\n" "${B}▸ $*${N}"; }
ok()   { printf "%s\n" "${G}✓ $*${N}"; }
warn() { printf "%s\n" "${Y}! $*${N}"; }
die()  { printf "%s\n" "${R}✗ $*${N}" >&2; exit 1; }
usage() { sed -n '3,/^set -euo/p' "$0" | sed '$d;s/^# \{0,1\}//'; }

# ---- helpers ----------------------------------------------------------------
workers() { [ -n "$WORKERS" ] && { echo "$WORKERS"; return; }; local n=$(( $(nproc) - 2 )); [ "$n" -lt 1 ] && n=1; echo "$n"; }
count()   { grep -cvE '^\s*#|^\s*$' "$1" 2>/dev/null || echo 0; }

preflight() {
  if ! curl -fsS -o /dev/null --max-time 10 "$BIOBTREE/ws/?i=TP53&mode=lite" 2>/dev/null; then
    die "biobtree unreachable at $BIOBTREE — start it or set ATLAS_BIOBTREE"
  fi
  ok "biobtree reachable at $BIOBTREE"
}

# build <label> <genes-list> <diseases-list> <drugs-list> → into $DIST
build() {
  local label="$1" g="$2" s="$3" r="$4" w
  for f in "$g" "$s" "$r"; do [ -f "$f" ] || die "missing list $f"; done
  w=$(workers)
  say "$label → $DIST  ${D}(genes=$(count "$g") diseases=$(count "$s") drugs=$(count "$r") | workers=$w${LIMIT:+ limit=$LIMIT})${N}"
  local t0=$SECONDS
  python -m atlas.batch --dist "$DIST" --workers "$w" \
    --genes "@$g" --diseases "@$s" --drugs "@$r" ${LIMIT:+--limit "$LIMIT"}
  ok "$label built in $((SECONDS - t0))s"
}

build_dense() { build "dense set" "$DENSE_DIR/genes.txt" "$DENSE_DIR/diseases.txt" "$DENSE_DIR/drugs.txt"; }

build_full() {
  [ -f "$SEED_DIR/genes_hgnc.txt" ] || { say "seeds missing → regenerating"; python scripts/build_corpus.py; }
  build "full corpus" "$SEED_DIR/genes_hgnc.txt" "$SEED_DIR/diseases_mondo_ranked.txt" "$SEED_DIR/drugs_chembl_approved.txt"
}

# Scope each pass to its own directory so neither reports the other as
# "deselected" — what runs is exactly what's collected.
unit()        { say "unit tests"; python -m pytest tests/unit; }
integration() {
  [ -d "$DIST/atlas" ] || die "no built dist at $DIST/atlas — run './atlas.sh test all' first"
  say "integration tests over $DIST"
  ATLAS_INTEGRATION_DIST="$DIST" python -m pytest tests/integration
}

archive() {
  local sha stamp out
  sha=$(git rev-parse --short HEAD 2>/dev/null || echo nogit)
  stamp=$(date +%Y%m%d-%H%M%S); out="$DIST/atlas-corpus-$stamp-$sha.tar.gz"
  say "archiving $DIST/atlas → $out"
  tar -C "$DIST" -czf "$out" atlas
  ok "$(du -h "$out" | cut -f1)  $out"
}

# ---- commands ---------------------------------------------------------------
cmd_test() {
  case "${1:-unit}" in
    unit)             unit ;;
    int|integration)  integration ;;
    all)              preflight; build_dense; echo; unit; echo; integration ;;
    *)                die "test: unknown mode '$1' (use: int | all, or no arg for unit)" ;;
  esac
  ok "tests passed"
}

# The production pipeline — runs inside the detached worker.
prod_run() {
  say "PRODUCTION BUILD — $(date)"
  preflight
  say "[1/4] pre-production check (dense build + tests)"
  build_dense; echo; unit; echo; integration
  ok "pre-production check green"
  echo; say "[2/4] full corpus"
  build_full
  echo; say "[3/4] integration sweep over full output"
  integration
  echo; say "[4/4] archive"
  archive
  ok "PRODUCTION BUILD COMPLETE — $(date)"
}

cmd_prod() {
  if [ "${ATLAS_PROD_WORKER:-}" != "1" ]; then          # launcher: detach + log
    mkdir -p logs
    local stamp log pid
    stamp=$(date +%Y%m%d-%H%M%S); log="logs/prod-$stamp.log"
    ln -sf "prod-$stamp.log" logs/prod.latest
    ATLAS_PROD_WORKER=1 ATLAS_DIST="$DIST" ATLAS_WORKERS="$WORKERS" ATLAS_LIMIT="$LIMIT" \
      nohup "$0" prod >"$log" 2>&1 &
    pid=$!
    ok "prod started — pid $pid"
    say "log: $log"
    printf "%s\n" "  ${D}follow:  tail -f $log   (or logs/prod.latest)${N}"
    return 0
  fi
  prod_run                                              # worker: run the pipeline
}

# ---- arg parse --------------------------------------------------------------
[ $# -ge 1 ] || { usage; exit 0; }
CMD="$1"; shift
POS=()
while [ $# -gt 0 ]; do
  case "$1" in
    -d|--dist)    DIST="$2"; shift 2 ;;
    -w|--workers) WORKERS="$2"; shift 2 ;;
    -n|--limit)   LIMIT="$2"; shift 2 ;;
    -h|--help)    usage; exit 0 ;;
    -*)           die "unknown option: $1" ;;
    *)            POS+=("$1"); shift ;;
  esac
done
set -- "${POS[@]:-}"

case "$CMD" in
  test)           cmd_test "${1:-}" ;;
  prod)           cmd_prod ;;
  help|-h|--help) usage ;;
  *)              die "unknown command: $CMD  (run './atlas.sh help')" ;;
esac
