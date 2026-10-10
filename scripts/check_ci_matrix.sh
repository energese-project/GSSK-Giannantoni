#!/bin/sh
# DEM-POR-01, its structural half: CI runs the whole catalogue on every toolchain, with -Werror.
#
# Verifies: NFR-DET-002, NFR-POR-001 (DEM-POR-01)
#
# The demonstration itself is a green CI run on the build matrix. This check makes sure that run
# means what it says. It fails when:
#   1. the build matrix loses one of macOS clang, Linux clang, Linux GCC;
#   2. CFLAGS stops carrying -std=c99 -Wall -Wextra -Werror;
#   3. a suite in the Makefile's CI_TESTS has no `make` step in deploy.yml. CI_TESTS is the list
#      contributors run locally (`make ci-local`); a suite in it that CI never runs is a test that
#      a green build does not cover. test-advanced is exempt: the Makefile records that the
#      pre-push checks run it and deploy.yml deliberately does not.
# POSIX sh + awk + grep.

set -eu
LC_ALL=C; export LC_ALL
cd "$(dirname "$0")/.."

WF=${CI_WORKFLOW:-.github/workflows/deploy.yml}
MK=${CI_MAKEFILE:-Makefile}
EXEMPT="test-advanced"
fail=0
err() { printf 'check-ci-matrix: %s\n' "$*" >&2; fail=1; }

for combo in "macos-latest clang" "ubuntu-latest clang" "ubuntu-latest gcc"; do
    os=${combo% *}; cc=${combo#* }
    grep -qE "\{ *os: *$os *, *cc: *$cc *\}" "$WF" || err "the build matrix has no $cc on $os (NFR-POR-001)"
done

flags=$(awk '/^CFLAGS[ \t]*=/ { sub(/^CFLAGS[ \t]*=[ \t]*/, ""); print; exit }' "$MK")
for f in -std=c99 -Wall -Wextra -Werror; do
    case " $flags " in *" $f "*) ;; *) err "CFLAGS lacks $f (NFR-POR-001)" ;; esac
done

suites=$(awk '/^CI_TESTS[ \t]*=/ { on = 1; sub(/^CI_TESTS[ \t]*=/, "") }
              on { line = $0; cont = sub(/\\[ \t]*$/, "", line); printf "%s ", line; if (!cont) exit }' "$MK")
[ -n "$suites" ] || err "no CI_TESTS list in $MK"
n=0
for t in $suites; do
    n=$((n + 1))
    case " $EXEMPT " in *" $t "*) continue ;; esac
    grep -qE "run:.*make[^#]*[ \t]$t([ \t]|\$)" "$WF" || err "$t is in CI_TESTS but no step in $WF runs it (NFR-DET-002)"
done

[ "$fail" -eq 0 ] || { echo "check-ci-matrix: FAILED" >&2; exit 1; }
echo "check-ci-matrix: OK -- 3 toolchains, -Werror, $n CI_TESTS suites each run by $WF"
