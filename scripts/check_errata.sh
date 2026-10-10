#!/bin/sh
# VAL-08 (BR-010): every erratum X1..Xn that PLAN.md §5 lists has a test that asserts the printed
# form fails its oracle, named for the erratum: an ok("Xn: ..."), close_c(...) or ord_case(...) whose label names Xn.
#
# Verifies: BR-010 (VAL-08)
#
# Fails when PLAN.md lists an erratum no test assertion names, or lists none at all.
# POSIX sh + grep + sed.

set -eu
LC_ALL=C; export LC_ALL
cd "$(dirname "$0")/.."

TESTS=${ERRATA_TESTS:-"tests/test_mop.c tests/test_mop_emergence.c tests/test_giannantoni.c"}
fail=0; n=0
ids=$(grep -oE '^\| X[0-9]+ \|' PLAN.md | grep -oE 'X[0-9]+' | sort -u -t X -k 2 -n)
[ -n "$ids" ] || { echo "check-errata: PLAN.md lists no errata" >&2; exit 1; }
for x in $ids; do
    n=$((n + 1))
    # shellcheck disable=SC2086
    if ! grep -hE '(ok|close_c|ord_case)\("' $TESTS | grep -qE "\"[^\"]*\b$x\b"; then
        echo "check-errata: $x has no test assertion naming it (VAL-08)" >&2; fail=1
    fi
done
[ "$fail" -eq 0 ] || { echo "check-errata: FAILED" >&2; exit 1; }
echo "check-errata: OK -- $n errata, each with a named assertion"
