#!/usr/bin/env bash
#
# atlas.sh — Sugi Atlas build & test orchestrator.
#
# Two verbs cover the whole loop. `test` is the inner loop (validate the
# generator); `prod` is the outer loop (ship the real corpus). Both build into
# a local, gitignored ./dist — nothing is pushed to a separate dist repo.
#
#   ./atlas.sh test [int|all]
#   ./atlas.sh sample [light|dense|both]
#   ./atlas.sh prod
#   ./atlas.sh release vX.Y.Z
#   ./atlas.sh release-publish vX.Y.Z
#   ./atlas.sh help
#
#   sample        Self-contained review corpus for the web/test team — built into
#                 dist/sample/ AND archived there as a tarball, so a tester can pull
#                 it without the full corpus. Two variants, same build+archive
#                 convention so both land in dist/sample/ consistently:
#                   light  (~30 pages, data/corpus/sample/*) → atlas-corpus-sample-light-<stamp>-<ver>.tar.gz
#                   dense (~1k pages, data/corpus/dense/*)   → atlas-corpus-sample-dense-<stamp>-<ver>.tar.gz
#                   both   builds light then dense; both tars persist in dist/sample/
#                 Default is light. Mesh is sample-local (links resolve only among
#                 the built pages) — expected for a sample.
#
#   test          Unit suite only — fast, no build. The default.
#   test int      Integration suite over ./dist (must already be built).
#   test all      Pre-production check: build the dense set into ./dist, then
#                 run unit + integration against it. Run this before `prod`.
#
#   release vX.Y.Z  Tag + GitHub release of the PIPELINE (code only — the corpus
#                 is never attached; we share the generator, not our outputs).
#                 Gate is `test all`; then `git tag` + push (over the existing SSH
#                 key). gh is OPTIONAL: if present it runs `gh release create`
#                 --generate-notes; if absent the tag is still pushed and the gh
#                 command is printed to run elsewhere (no token needed on prod).
#                 Requires a committed tree. Build the corpus off the tag (`prod`)
#                 to stamp pages with vX.Y.Z.
#
#   release-publish vX.Y.Z  Create the GitHub release for an already-pushed tag.
#                 Run this where gh is installed + authenticated (e.g. a laptop);
#                 the build box only needs `release` (tag + push, no token). Tag
#                 must already exist on origin.
#
#   prod          Production build, detached via nohup (logs → dist/logs/):
#                   1. pre-production check  (== test all: dense build + tests)
#                   2. full corpus from data/corpus/seeds/ into ./dist
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
DENSE_DIR="data/corpus/dense"
SEED_DIR="data/corpus/seeds"
SAMPLE_DIR="data/corpus/sample"

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

# Capture the pipeline version + commit ONCE at build start (git-derived,
# enju-style) and export them so the batch stamps every page with the code that
# actually built it, and the archive name matches. A tagged commit → a clean
# "vX.Y.Z"; between tags → "vX.Y.Z-N-gSHA"; no git → "dev"/"nogit". Idempotent
# (`:=` only sets when unset), so all phases of one run share the same stamp.
stamp_version() {
  : "${ATLAS_BUILD_VERSION:=$(git describe --tags --always 2>/dev/null || echo dev)}"
  : "${ATLAS_BUILD_COMMIT:=$(git rev-parse --short HEAD 2>/dev/null || echo nogit)}"
  export ATLAS_BUILD_VERSION ATLAS_BUILD_COMMIT
}

preflight() {
  if ! curl -fsS -o /dev/null --max-time 10 "$BIOBTREE/ws/?i=TP53&mode=lite" 2>/dev/null; then
    die "biobtree unreachable at $BIOBTREE — start it or set ATLAS_BIOBTREE"
  fi
  ok "biobtree reachable at $BIOBTREE"
}

# build <label> <genes-list> <diseases-list> <drugs-list> [<pathways-list>] → into $DIST
build() {
  local label="$1" g="$2" s="$3" r="$4" p="${5:-}" w pcount=0
  for f in "$g" "$s" "$r"; do [ -f "$f" ] || die "missing list $f"; done
  [ -n "$p" ] && [ -f "$p" ] && pcount=$(count "$p")
  w=$(workers)
  stamp_version                       # version/commit captured at build start
  # Start every build from a clean slate so pages from a prior phase/run can't
  # linger as orphans. This matters in `prod`: the [1/4] dense check populates
  # dist/atlas, then [2/4] full corpus must NOT inherit those dense-only pages
  # (e.g. GBA, which the full HGNC seed lists as GBA1, or salt-form drugs) — they
  # would fail the manifest⇄pages bijection.
  rm -rf "$DIST/atlas" "$DIST/cache"
  say "$label → $DIST  ${D}(genes=$(count "$g") diseases=$(count "$s") drugs=$(count "$r") pathways=$pcount | workers=$w${LIMIT:+ limit=$LIMIT})${N}"
  local t0=$SECONDS
  python -m atlas.batch --dist "$DIST" --workers "$w" \
    --genes "@$g" --diseases "@$s" --drugs "@$r" \
    ${p:+--pathways "@$p"} ${LIMIT:+--limit "$LIMIT"}
  ok "$label built in $((SECONDS - t0))s"
}

# NOTE: pathways are SHELVED (see docs/internal/IMPROVEMENT_ROADMAP.md thread D) —
# the vertical code is retained but dormant; not seeded into the dense set or prod.
# To revive, add a pathways list as the 5th build() arg (build() still supports it).
build_dense() { build "dense set" "$DENSE_DIR/genes.txt" "$DENSE_DIR/diseases.txt" "$DENSE_DIR/drugs.txt"; }

build_full() {
  [ -f "$SEED_DIR/drugs_chembl_seed.txt" ] || { say "seeds missing → regenerating"; python -m atlas.build_corpus; }
  build "full corpus" "$SEED_DIR/genes_hgnc.txt" "$SEED_DIR/diseases_mondo_ranked.txt" "$SEED_DIR/drugs_chembl_seed.txt"
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
  local infix="${1:+$1-}" stamp out   # optional name infix, e.g. "sample" → atlas-corpus-sample-…
  stamp_version                       # same stamp as the build (not archive-time HEAD)
  stamp=$(date +%Y%m%d-%H%M%S); out="$DIST/atlas-corpus-${infix}$stamp-$ATLAS_BUILD_VERSION.tar.gz"
  say "archiving $DIST/atlas → $out"
  tar -C "$DIST" -czf "$out" atlas
  ok "$(du -h "$out" | cut -f1)  $out"
}

# Sample corpus for the web/test team — built into dist/sample/ AND archived
# there as a tarball, so a tester can pull a self-contained set without the full
# corpus. Two variants share one build+archive convention so both tars land in
# dist/sample/ consistently (light ~30 pages, dense ~1k pages); `both` produces
# both tars side by side. The cross-entity mesh is sample-local (links resolve
# only among the built pages) — expected for a sample.
_sample_variant() {   # <label> <archive-infix> <genes> <diseases> <drugs>  (DIST pre-scoped to dist/sample)
  build "$1" "$3" "$4" "$5"
  echo; integration
  echo; archive "$2"
}

cmd_sample() {
  local mode="${1:-light}"
  preflight
  local DIST="$DIST/sample"           # build + archive under dist/sample/ (function-scoped)
  case "$mode" in
    light) _sample_variant "sample set (light, ~30 pages)" sample-light \
             "$SAMPLE_DIR/genes.txt" "$SAMPLE_DIR/diseases.txt" "$SAMPLE_DIR/drugs.txt" ;;
    dense) _sample_variant "sample set (dense)" sample-dense \
             "$DENSE_DIR/genes.txt" "$DENSE_DIR/diseases.txt" "$DENSE_DIR/drugs.txt" ;;
    both)  _sample_variant "sample set (light, ~30 pages)" sample-light \
             "$SAMPLE_DIR/genes.txt" "$SAMPLE_DIR/diseases.txt" "$SAMPLE_DIR/drugs.txt"
           echo
           _sample_variant "sample set (dense)" sample-dense \
             "$DENSE_DIR/genes.txt" "$DENSE_DIR/diseases.txt" "$DENSE_DIR/drugs.txt" ;;
    *)     die "usage: ./atlas.sh sample [light|dense|both]" ;;
  esac
  ok "sample corpus ready ($mode) — tar(s) in $DIST/"
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
    local logdir="$DIST/logs" stamp log pid
    mkdir -p "$logdir"
    stamp=$(date +%Y%m%d-%H%M%S); log="$logdir/prod-$stamp.log"
    ln -sf "prod-$stamp.log" "$logdir/prod.latest"
    ATLAS_PROD_WORKER=1 ATLAS_DIST="$DIST" ATLAS_WORKERS="$WORKERS" ATLAS_LIMIT="$LIMIT" \
      nohup "$0" prod >"$log" 2>&1 &
    pid=$!
    ok "prod started — pid $pid"
    say "log: $log"
    printf "%s\n" "  ${D}follow:  tail -f $log   (or $logdir/prod.latest)${N}"
    return 0
  fi
  prod_run                                              # worker: run the pipeline
}

# release vX.Y.Z — tag + GitHub release of the PIPELINE (code only; the corpus is
# never attached — we share the generator, not our outputs). Gate is `test all`
# (dense build + unit + integration) so a release provably builds end-to-end. The
# version flows from the tag: a subsequent `prod` off this commit stamps pages
# with the clean tag (see atlas_version() — no version-bump commit moves HEAD off
# the tag, so the corpus reads a pristine vX.Y.Z).
cmd_release() {
  local version="${1:-}"
  [ -n "$version" ] || die "usage: ./atlas.sh release <version>   (e.g. v1.0.0)"
  [[ "$version" == v* ]] || version="v$version"
  # gh is OPTIONAL: the tag + push (what actually versions the corpus) runs here
  # with the existing SSH key, so no GitHub token needs to live on this box. The
  # GitHub-release publish can run on any machine that has gh (e.g. a laptop).
  # The tag captures HEAD, so the tree must be committed. The untracked tmp
  # working doc (corpus-stats.tmp.md) is allowed to stay dirty.
  local dirty; dirty=$(git status --porcelain | grep -vE 'corpus-stats\.tmp\.md' || true)
  [ -z "$dirty" ] || die "uncommitted changes (commit before releasing — the tag captures HEAD):
$dirty"
  git rev-parse "$version" >/dev/null 2>&1 && die "tag $version already exists"
  say "RELEASE $version — gate: test all (dense build + unit + integration)"
  preflight; build_dense; echo; unit; echo; integration
  ok "gate green"
  echo; say "tagging $version at $(git rev-parse --short HEAD)"
  git tag -a "$version" -m "Sugi Atlas pipeline $version"
  git push origin "$version"
  ok "tag $version pushed to origin"
  echo
  if command -v gh >/dev/null 2>&1; then
    say "creating GitHub release (pipeline source only)"
    gh release create "$version" --title "$version" --generate-notes
    ok "RELEASE $version published — pipeline only, no corpus attached"
  else
    warn "gh not installed here — the tag is pushed. Publish the GitHub release on a gh machine (e.g. your laptop):"
    printf "%s\n" "    ${B}./atlas.sh release-publish $version${N}"
  fi
  echo; say "build the corpus off this tag (./atlas.sh prod) to stamp pages with $version"
}

# release-publish vX.Y.Z — create the GitHub release for an ALREADY-pushed tag.
# The counterpart to `release` for the split workflow: `release` tags + pushes on
# the build box (where biobtree is, no gh/token needed); `release-publish` runs
# wherever gh is installed + authenticated (e.g. a laptop). It refuses to invent a
# tag — gh would otherwise create one at the default-branch HEAD if it's missing.
cmd_release_publish() {
  local version="${1:-}"
  [ -n "$version" ] || die "usage: ./atlas.sh release-publish <version>   (e.g. v1.0.0)"
  [[ "$version" == v* ]] || version="v$version"
  command -v gh >/dev/null 2>&1 || die "gh (GitHub CLI) not found — run this where gh is installed (e.g. your laptop)"
  gh auth status >/dev/null 2>&1 || die "gh not authenticated — run 'gh auth login' first"
  git ls-remote --tags origin "refs/tags/$version" | grep -q . \
    || die "tag $version not on origin — run './atlas.sh release $version' on the build box first"
  gh release view "$version" >/dev/null 2>&1 && die "GitHub release $version already exists"
  say "publishing GitHub release $version (pipeline source only — no corpus attached)"
  gh release create "$version" --title "$version" --generate-notes
  ok "RELEASE $version published"
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
  test)            cmd_test "${1:-}" ;;
  sample)          cmd_sample "${1:-}" ;;
  prod)            cmd_prod ;;
  release)         cmd_release "${1:-}" ;;
  release-publish) cmd_release_publish "${1:-}" ;;
  help|-h|--help)  usage ;;
  *)              die "unknown command: $CMD  (run './atlas.sh help')" ;;
esac
