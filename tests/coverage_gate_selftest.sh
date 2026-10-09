#!/bin/sh
# Self-test of scripts/coverage_gate.sh: the gate must bite.
# Verifies: NFR-COV-001 (T-COV-01)
#
# The old gate printed OK on input it could not parse and parsed with `grep -P` (PLAN §1 B3). Each
# case below is a report the gate must reject or accept; the full run on the real units is
# `make coverage-gia`, in the Linux quality gate.

set -eu
LC_ALL=C; export LC_ALL
cd "$(dirname "$0")/.."
G=scripts/coverage_gate.sh

tmp=$(mktemp -d "${TMPDIR:-/tmp}/covgate.XXXXXX")
trap 'rm -rf "$tmp"' EXIT INT TERM
failures=0
expect() {  # expect <what> <want: pass|fail> <report-text>
    printf '%s\n' "$3" > "$tmp/r"
    if sh "$G" 90 "$tmp/r" > "$tmp/out" 2>&1; then got=pass; else got=fail; fi
    if [ "$got" = "$2" ]; then printf '  %-58s PASS\n' "$1"
    else printf '  %-58s FAIL (gate said %s)\n' "$1" "$got"; cat "$tmp/out"; failures=$((failures + 1)); fi
}

echo "=== coverage gate self-test ==="
expect "90.00% passes a 90% gate"           pass "File 'src/engine.c'
Lines executed:90.00% of 1000"
expect "89.00% fails a 90% gate"            fail "File 'src/engine.c'
Lines executed:89.00% of 1000"
expect "weighted by line count, not averaged" fail "File 'src/engine.c'
Lines executed:95.00% of 100
File 'src/mop.c'
Lines executed:85.00% of 900"
expect "two units together above the gate"  pass "File 'src/engine.c'
Lines executed:92.00% of 500
File 'src/mop.c'
Lines executed:90.00% of 500"
expect "garbage fails"                      fail "this is not gcov output"
expect "empty report fails"                 fail ""
expect "a File with no Lines line fails"    fail "File 'src/engine.c'"
expect "a non-numeric percentage fails"     fail "File 'src/engine.c'
Lines executed:n/a% of 1000"

echo
if [ "$failures" -eq 0 ]; then echo "ALL PASS"; else echo "FAILURES PRESENT"; fi
echo "failures: $failures"
[ "$failures" -eq 0 ]
