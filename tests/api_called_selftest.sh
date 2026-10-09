#!/bin/sh
# Self-test of scripts/check_api_called.sh (T-API-02): the check must bite.
#
# The fixture headers are the C blocks of IF-API-002, IF-API-003 and IF-API-004, copied from
# docs/requirements/icd.md, so the declaration parser is tested on the headers W2-W8 will write.

set -eu
LC_ALL=C; export LC_ALL
cd "$(dirname "$0")/.."

tmp=$(mktemp -d "${TMPDIR:-/tmp}/apisel.XXXXXX")
trap 'rm -rf "$tmp"' EXIT INT TERM
failures=0
ok() {
    if [ "$2" -eq 0 ]; then printf '  %-62s PASS\n' "$1"
    else printf '  %-62s FAIL\n' "$1"; failures=$((failures + 1)); fi
}
check() {  # check <headers> <sources> <tests> : status of the check, its stderr in $tmp/err
    set +e
    API_HEADERS="$1" API_SOURCES="$2" API_TESTS="$3" API_ENGINE=/nonexistent \
        sh scripts/check_api_called.sh > "$tmp/out" 2> "$tmp/err"
    rc=$?
    set -e
    return $rc
}

# The three API blocks of icd.md, as headers.
for id in IF-API-002 IF-API-003 IF-API-004; do
    awk -v id="### $id" '
        index($0, id) == 1 { want = 1; next }
        want && /^```c/   { inblk = 1; next }
        inblk && /^```/   { exit }
        inblk             { print }
    ' docs/requirements/icd.md > "$tmp/$id.h"
done

echo "=== check-api-called self-test ==="

# Every function each block declares, by hand from icd.md.
cat > "$tmp/want" <<'NAMES'
gia_idc_of gia_lde2_solve gia_lde2_eval gia_lde2_terms gia_lde2_free gia_binary_solve gia_binary_eval
gia_riccati_solve gia_riccati_eval gia_nl1410_roots gia_idc_taylor gia_idc_refuse
rel_mul rel_mul3 rel_exp rel_root rel_root_pow rel_mul_pow
gia_mop_couple gia_mop_solve gia_matrioska_free gia_mop_second gia_eqs gia_mop_network_beta
gia_verdict_str gia_harmony_residual gia_harmony_verdict gia_harmony_observed
NAMES
tr ' ' '\n' < "$tmp/want" | grep . | sort > "$tmp/want.sorted"

# A test file naming every one of them passes.
{ echo "/* calls */"; tr '\n' ' ' < "$tmp/want.sorted"; echo; } > "$tmp/all_called.c"
st=0; check "$tmp/IF-API-002.h $tmp/IF-API-003.h $tmp/IF-API-004.h" "" "$tmp/all_called.c" || st=1
ok "icd.md headers, every function called: passes" "$st"
grep -q "28 declared" "$tmp/out" || { echo "    parser found: $(cat "$tmp/out")" >&2; st=1; }
ok "parser finds exactly the 28 functions icd.md declares" "$st"

# Drop each name in turn: the check must name it.
st=0
for f in $(cat "$tmp/want.sorted"); do
    grep -vx "$f" "$tmp/want.sorted" | tr '\n' ' ' > "$tmp/one_missing.c"
    if check "$tmp/IF-API-002.h $tmp/IF-API-003.h $tmp/IF-API-004.h" "" "$tmp/one_missing.c"; then
        echo "    $f uncalled, check passed" >&2; st=1
    elif ! grep -q "declares $f," "$tmp/err"; then
        echo "    $f uncalled, check did not name it" >&2; st=1
    fi
done
ok "each of the 28, left uncalled, fails the check by name" "$st"

# Typedefs, including function-pointer typedefs, are not functions.
grep -q "gia_cfn\|gia_construction" "$tmp/err" && st=1 || st=0
ok "typedefs (gia_cfn, gia_construction) are not counted" "$st"

# A stub marker in a unit source fails.
printf 'gia_status gia_idc_of(void) { return GIA_E_UNIMPLEMENTED; }\n' > "$tmp/stub.c"
st=1; check "$tmp/IF-API-002.h" "$tmp/stub.c" "$tmp/all_called.c" || st=0
ok "a not-implemented return in a unit fails" "$st"

# With no test file at all, a declared function fails.
st=1; check "$tmp/IF-API-003.h" "" "$tmp/no_such_test.c" || st=0
ok "no test file: every declaration fails" "$st"

# The repository as it stands passes.
st=0; sh scripts/check_api_called.sh > /dev/null 2>&1 || st=1
ok "repository: passes" "$st"

echo
if [ "$failures" -eq 0 ]; then echo "ALL PASS"; else echo "FAILURES PRESENT"; fi
echo "failures: $failures"
[ "$failures" -eq 0 ]
