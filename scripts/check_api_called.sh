#!/bin/sh
# Every declared function of the Giannantoni API is called by a test, and none is a stub
# (NFR-API-001, PLAN.md G3, ADR 0018 rule 3; vv-plan.md T-API-02).
#
# Verifies: NFR-API-001, IF-API-002, IF-API-003, IF-API-004, IF-API-005 (T-API-02)
#
# Fails when:
#   1. a function declared in one of the headers is named by no test source;
#   2. a unit source contains a not-implemented marker (AGENTS.md's `(void)`-stub allowance does
#      not cover these headers).
# A header or source that does not exist yet is skipped by name, so W2-W8 join the check by landing.
#
# Overridable for the self-test (tests/api_called_selftest.sh): API_HEADERS, API_SOURCES, API_TESTS.
# POSIX sh + awk + grep.

set -eu
LC_ALL=C; export LC_ALL
cd "$(dirname "$0")/.."

API_HEADERS=${API_HEADERS:-"include/gia_status.h include/idc.h include/relational.h include/mop.h include/mop_seed.h"}
# engine.h predates the baseline and is not NFR-API-001's; only the functions IF-API-005 adds to it
# are checked, and each only once it is declared.
API_ENGINE=${API_ENGINE:-include/engine.h}
API_ENGINE_FUNCS=${API_ENGINE_FUNCS:-"gia_ordinality_record gia_at_maximum_ordinality gia_closure gia_solution_drift gia_drift_projection gia_emergy_source_term gia_emergy_check_limits gia_emergy_carried gia_emergy_global_balance gia_emergy_balance_solve gia_oform_scalar gia_oform_binary gia_oform_duet gia_oform_duet_binary gia_circle_product gia_circle_reduce"}
API_SOURCES=${API_SOURCES:-"src/idc.c src/relational.c src/mop.c src/mop_seed.c src/harmony.c"}
API_TESTS=${API_TESTS:-"tests/test_giannantoni.c tests/test_mop.c tests/test_mop_emergence.c tests/test_mop_threads.c"}

fail=0
err() { printf 'check-api-called: %s\n' "$*" >&2; fail=1; }

tests=""
for t in $API_TESTS; do [ -f "$t" ] && tests="$tests $t"; done

# Function names declared in a header: comments and preprocessor lines removed, statements split
# on ';' outside braces, typedefs skipped, and the identifier right before the first '(' taken.
declared() {
    awk '
        { line = $0 }
        incom { if (match(line, /\*\//)) { line = substr(line, RSTART + 2); incom = 0 } else next }
        {
            while (match(line, /\/\*/)) {
                pre = substr(line, 1, RSTART - 1); rest = substr(line, RSTART + 2)
                if (match(rest, /\*\//)) line = pre " " substr(rest, RSTART + 2)
                else { line = pre; incom = 1; break }
            }
            sub(/\/\/.*/, "", line)
            if (line ~ /^[ \t]*#/) next
            buf = buf " " line
        }
        END {
            n = length(buf); stmt = ""; depth = 0
            for (i = 1; i <= n; i++) {
                c = substr(buf, i, 1)
                if (c == "{") { depth++; continue }
                if (c == "}") { depth--; continue }
                if (depth > 0) continue
                if (c == ";") { emit(stmt); stmt = ""; continue }
                stmt = stmt c
            }
        }
        function emit(s) {
            sub(/^[ \t]+/, "", s)
            if (s ~ /^typedef[ \t]/ || s ~ /^extern[ \t]+"C"/) return
            if (!match(s, /[A-Za-z_][A-Za-z0-9_]*[ \t]*\(/)) return
            name = substr(s, RSTART, RLENGTH); sub(/[ \t]*\($/, "", name)
            print name
        }
    ' "$1"
}

n_fun=0
for h in $API_HEADERS; do
    [ -f "$h" ] || continue
    for f in $(declared "$h" | sort -u); do
        n_fun=$((n_fun + 1))
        # shellcheck disable=SC2086
        if [ -z "$tests" ] || ! grep -qw "$f" $tests; then
            err "$h declares $f, and no test calls it (NFR-API-001)"
        fi
    done
done

if [ -f "$API_ENGINE" ]; then
    declared "$API_ENGINE" | sort -u > "${TMPDIR:-/tmp}/api_engine.$$"
    for f in $API_ENGINE_FUNCS; do
        grep -qx "$f" "${TMPDIR:-/tmp}/api_engine.$$" || continue
        n_fun=$((n_fun + 1))
        # shellcheck disable=SC2086
        if [ -z "$tests" ] || ! grep -qw "$f" $tests; then
            err "$API_ENGINE declares $f (IF-API-005), and no test calls it (NFR-API-001)"
        fi
    done
    rm -f "${TMPDIR:-/tmp}/api_engine.$$"
fi

for s in $API_SOURCES; do
    [ -f "$s" ] || continue
    if grep -nEi 'UNIMPLEMENTED|NOT_?IMPLEMENTED|not yet implemented|ENOSYS' "$s" > /dev/null; then
        grep -nEi 'UNIMPLEMENTED|NOT_?IMPLEMENTED|not yet implemented|ENOSYS' "$s" | while IFS= read -r l; do
            printf 'check-api-called: %s:%s -- a stub in planned API (ADR 0018 rule 3)\n' "$s" "$l" >&2
        done
        fail=1
    fi
done

if [ "$fail" -ne 0 ]; then echo "check-api-called: FAILED" >&2; exit 1; fi
echo "check-api-called: OK -- $n_fun declared function(s), each called by a test"
