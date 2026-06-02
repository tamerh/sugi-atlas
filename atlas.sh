#!/usr/bin/env bash
#
# atlas.sh — Sugi Atlas orchestrator.
#
# One entry point for the everyday loops: run the tests, build the dense
# reference set, gate a release, refresh the corpus seeds, or kick off the full
# corpus. It only orchestrates — the real work lives in `python -m atlas.batch`,
# `pytest`, and `scripts/build_corpus.py`; this just wires them together with
# sane defaults and a uniform interface.
#
#   ./atlas.sh <command> [options]
#
# Commands:
#   test [unit|int|all]   Run tests (default: unit).
#                           unit → pytest -m "not integration"  (fast, no dist)
#                           int  → pytest -m integration         (over $DIST)
#                           all  → unit then int
#   dense                 Build the dense reference set into $DIST
#                           (corpus/dense/{genes,diseases,drugs}.txt).
#   gate                  Release gate: dense build → integration tests.
#                           Green here is the precondition for a full corpus.
#   seeds [genes|drugs|diseases]
#                         Regenerate the full-corpus seed lists from biobtree's
#                           own ingest sources (corpus/seeds/, gitignored).
#   corpus                Full corpus build from corpus/seeds/ → tar.gz archive
#                           (the full corpus is archived, NOT committed to git).
#   health                Probe biobtree reachability ($ATLAS_BIOBTREE).
#   help                  This message.
#
# Options (override the defaults; also honoured from the environment):
#   -d, --dist DIR        output dist            (default $ATLAS_DIST or /data/sugi-atlas-dist)
#   -w, --workers N       parallel workers       (default: nproc-2)
#   -n, --limit N         build only first N of each type (test slice)
#   -a, --archive DIR     archive dir for `corpus` tar.gz (default: $DIST/../atlas-archives)
#   -h, --help            show help
#
# Examples:
#   ./atlas.sh test                 # fast unit suite
#   ./atlas.sh test all             # unit + integration over the current dist
#   ./atlas.sh dense -n 20          # quick 20-of-each smoke build
#   ./atlas.sh gate                 # the pre-release loop: dense then int tests
#   ./atlas.sh seeds drugs          # refresh just the drug seed list
#   ./atlas.sh corpus -w 24         # full corpus, 24 workers, archived
#
set -euo pipefail
cd "$(dirname "$0")"

# ---- defaults ---------------------------------------------------------------
DIST="${ATLAS_DIST:-/data/sugi-atlas-dist}"
WORKERS=""
LIMIT=""
ARCHIVE=""
BIOBTREE="${ATLAS_BIOBTREE:-http://127.0.0.1:9291}"

CORPUS_DENSE="corpus/dense"
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
default_workers() { local n; n=$(( $(nproc) - 2 )); [ "$n" -lt 1 ] && n=1; echo "$n"; }

count() { grep -cvE '^\s*#|^\s*$' "$1" 2>/dev/null || echo 0; }

# Build a list-dir's {genes,diseases,drugs}.txt into $DIST via the parallel
# batch driver. $1=list dir, $2=human label.
build_lists() {
  local dir="$1" label="$2"
  local g="$dir/genes.txt" s="$dir/diseases.txt" r="$dir/drugs.txt"
  for f in "$g" "$s" "$r"; do
    [ -f "$f" ] || die "missing list $f — (regenerate with: ./atlas.sh seeds)"
  done
  local w="${WORKERS:-$(default_workers)}"
  say "$label → ${B}$DIST${N}"
  printf "%s\n" "  ${D}genes=$(count "$g")  diseases=$(count "$s")  drugs=$(count "$r")  | workers=$w${LIMIT:+  limit=$LIMIT}${N}"
  local t0 t1
  t0=$(date +%s)
  python -m atlas.batch \
    --dist "$DIST" --workers "$w" \
    --genes "@$g" --diseases "@$s" --drugs "@$r" \
    ${LIMIT:+--limit "$LIMIT"}
  t1=$(date +%s)
  ok "$label built in $((t1 - t0))s"
}

run_integration() {
  [ -d "$DIST/atlas" ] || die "no built dist at $DIST/atlas — run './atlas.sh dense' first"
  say "integration tests over ${B}$DIST${N}"
  ATLAS_INTEGRATION_DIST="$DIST" python -m pytest -m integration
}

run_unit() {
  say "unit tests"
  python -m pytest -m "not integration"
}

# ---- commands ---------------------------------------------------------------
cmd_test() {
  case "${1:-unit}" in
    unit)        run_unit ;;
    int|integration) run_integration ;;
    all)         run_unit; echo; run_integration ;;
    *)           die "test: unknown mode '${1}' (use: unit | int | all)" ;;
  esac
  ok "tests passed"
}

cmd_dense() { build_lists "$CORPUS_DENSE" "dense set"; }

cmd_gate() {
  say "release gate: dense build → integration tests"
  build_lists "$CORPUS_DENSE" "dense set"
  echo
  run_integration
  ok "gate green — dense build + integration tests pass"
}

cmd_seeds() {
  say "regenerating corpus seeds → ${B}$SEED_DIR${N}"
  python scripts/build_corpus.py ${1:+--only "$1"}
  ok "seeds written (gitignored; regenerable)"
}

cmd_corpus() {
  local g="$SEED_DIR/genes_hgnc.txt"
  local s="$SEED_DIR/diseases_mondo_ranked.txt"
  local r="$SEED_DIR/drugs_chembl_approved.txt"
  for f in "$g" "$s" "$r"; do
    [ -f "$f" ] || die "missing seed $f — run './atlas.sh seeds' first"
  done
  local w="${WORKERS:-$(default_workers)}"
  local arch="${ARCHIVE:-$(dirname "$DIST")/atlas-archives}"
  warn "FULL CORPUS — genes=$(count "$g") diseases=$(count "$s") drugs=$(count "$r") | workers=$w"
  warn "output → $DIST   archive → $arch   (NOT committed to git)"
  local t0 t1
  t0=$(date +%s)
  python -m atlas.batch \
    --dist "$DIST" --workers "$w" \
    --genes "@$g" --diseases "@$s" --drugs "@$r" \
    ${LIMIT:+--limit "$LIMIT"}
  t1=$(date +%s)
  ok "corpus built in $(( (t1 - t0) / 60 ))m$(( (t1 - t0) % 60 ))s"
  mkdir -p "$arch"
  local out="$arch/atlas-corpus-$(date +%Y%m%d-%H%M%S).tar.gz"
  say "archiving $DIST/atlas → $out"
  tar -C "$DIST" -czf "$out" atlas
  ok "$(du -h "$out" | cut -f1)  $out"
}

cmd_health() {
  say "biobtree at ${B}$BIOBTREE${N}"
  local url="$BIOBTREE/ws/?i=TP53&mode=lite"
  if code=$(curl -fsS -o /dev/null -w '%{http_code}' --max-time 10 "$url" 2>/dev/null); then
    ok "reachable (HTTP $code, probe TP53)"
  else
    die "unreachable — start biobtree or set ATLAS_BIOBTREE (gate: ATLAS_BIOBTREE=http://127.0.0.1:8000)"
  fi
}

# ---- arg parse --------------------------------------------------------------
[ $# -ge 1 ] || { usage; exit 0; }
CMD="$1"; shift

POSITIONAL=()
while [ $# -gt 0 ]; do
  case "$1" in
    -d|--dist)    DIST="$2"; shift 2 ;;
    -w|--workers) WORKERS="$2"; shift 2 ;;
    -n|--limit)   LIMIT="$2"; shift 2 ;;
    -a|--archive) ARCHIVE="$2"; shift 2 ;;
    -h|--help)    usage; exit 0 ;;
    -*)           die "unknown option: $1" ;;
    *)            POSITIONAL+=("$1"); shift ;;
  esac
done
set -- "${POSITIONAL[@]:-}"

case "$CMD" in
  test)            cmd_test "${1:-}" ;;
  dense)           cmd_dense ;;
  gate)            cmd_gate ;;
  seeds)           cmd_seeds "${1:-}" ;;
  corpus|full)     cmd_corpus ;;
  health)          cmd_health ;;
  help|-h|--help)  usage ;;
  *)               die "unknown command: $CMD  (run './atlas.sh help')" ;;
esac
