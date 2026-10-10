#!/bin/sh
# Self-test of the no-skip guard (PLAN.md G1, ADR 0018 rule 1, task guard-no-skip).
#
# `make test` used to print SKIPPED and exit 0 for a model with no golden file (PLAN §1 B1), so a
# model "passed" by having no expected output. This runs `make test` against doctored copies of
# tests/expected and the allowlist and checks each way the guard must respond.

set -eu
LC_ALL=C; export LC_ALL
cd "$(dirname "$0")/.."

tmp=$(mktemp -d "${TMPDIR:-/tmp}/guard_no_skip.XXXXXX")
trap 'rm -rf "$tmp"' EXIT INT TERM
failures=0
ok() {
    if [ "$2" -eq 0 ]; then printf '  %-62s PASS\n' "$1"
    else printf '  %-62s FAIL\n' "$1"; failures=$((failures + 1)); fi
}
mt() {  # mt <expected-dir> <allowlist> : runs make test, output in $tmp/log, status echoed
    set +e
    ${MAKE:-make} -s test EXPECTED_DIR="$1" SKIP_ALLOWLIST="$2" > "$tmp/log" 2>&1
    echo $?
    set -e
}

echo "=== guard-no-skip self-test ==="
VICTIM=decay_model
cp -R tests/expected "$tmp/exp"
rm "$tmp/exp/$VICTIM.csv"
# Every doctored allowlist starts from the real one, so only the victim's entry differs.

# 1. A missing golden file fails make test.
cp tests/skip_allowlist.txt "$tmp/real"
rc=$(mt "$tmp/exp" "$tmp/real")
st=1; [ "$rc" -ne 0 ] && grep -q "Testing $VICTIM... FAILED (no expected output" "$tmp/log" && st=0
ok "missing golden file, not allowlisted: make test fails" "$st"

# 2. Allowlisted with a reason: skipped, and make test passes.
{ cat tests/skip_allowlist.txt; printf '%s  # self-test: allowlisted on purpose\n' "$VICTIM"; } > "$tmp/allow"
rc=$(mt "$tmp/exp" "$tmp/allow")
st=1; [ "$rc" -eq 0 ] && grep -q "Testing $VICTIM... SKIPPED (allowlisted" "$tmp/log" && st=0
ok "missing golden file, allowlisted with a reason: skipped" "$st"

# 3. An allowlist entry without a reason is not an allowlist entry.
{ cat tests/skip_allowlist.txt; printf '%s\n' "$VICTIM"; } > "$tmp/allow_bare"
rc=$(mt "$tmp/exp" "$tmp/allow_bare")
st=1; [ "$rc" -ne 0 ] && grep -q "Testing $VICTIM... FAILED (no expected output, and the" "$tmp/log" && st=0
ok "allowlist entry with no reason: make test fails" "$st"

# 4. A stale entry -- allowlisted, yet the golden file exists -- fails, so the list cannot grow
#    silently past what it needs.
rc=$(mt tests/expected "$tmp/allow")
st=1; [ "$rc" -ne 0 ] && grep -q "Testing $VICTIM... FAILED (stale entry" "$tmp/log" && st=0
ok "allowlisted model that has a golden file: make test fails" "$st"

# 5. The real configuration passes.
rc=$(mt tests/expected tests/skip_allowlist.txt)
ok "repository's own expected files and allowlist: make test passes" "$rc"

echo
if [ "$failures" -eq 0 ]; then echo "ALL PASS"; else echo "FAILURES PRESENT"; cat "$tmp/log"; fi
echo "failures: $failures"
[ "$failures" -eq 0 ]
