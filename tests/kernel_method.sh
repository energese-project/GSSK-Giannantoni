#!/bin/sh
# T-KER-01 (FR-KER-001, ADR 0022): the kernel's matrix-exponential method is `expm`, and
# `incipient` is a deprecated alias for the same code path.
#
# Verifies: FR-KER-001 (T-KER-01)
#
# Oracle (ADR 0022 decision 5): the same model run with "incipient" and with "expm" writes a
# byte-identical CSV; the "incipient" run prints the deprecation notice exactly once on stderr,
# naming expm, and the "expm" run prints none; each serialises back the spelling it was loaded
# with. POSIX sh.

set -eu
LC_ALL=C; export LC_ALL
cd "$(dirname "$0")/.."

GSSK=${GSSK:-bin/gssk}
DUMP=${DUMP:-bin/dump_serialized}
for b in "$GSSK" "$DUMP"; do [ -x "$b" ] || { echo "kernel_method: $b not built" >&2; exit 1; }; done

tmp=$(mktemp -d "${TMPDIR:-/tmp}/kernel_method.XXXXXX")
trap 'rm -rf "$tmp"' EXIT INT TERM
failures=0
ok() {
    if [ "$2" -eq 0 ]; then printf '  %-62s PASS\n' "$1"
    else printf '  %-62s FAIL\n' "$1"; failures=$((failures + 1)); fi
}

echo "=== kernel method label (T-KER-01) ==="
for model in examples/decay_model.json examples/diffusion_model.json; do
    name=$(basename "$model" .json)
    for m in incipient expm; do
        sed 's/"method":[[:space:]]*"rk4"/"method": "'"$m"'"/' "$model" > "$tmp/$name.$m.json"
        grep -q "\"method\": \"$m\"" "$tmp/$name.$m.json" || { echo "    could not set method in $model" >&2; exit 1; }
        set +e
        "$GSSK" "$tmp/$name.$m.json" "$tmp/$name.$m.csv" > /dev/null 2> "$tmp/$name.$m.err"
        echo $? > "$tmp/$name.$m.rc"
        set -e
    done
    st=0
    [ "$(cat "$tmp/$name.incipient.rc")" -eq 0 ] && [ "$(cat "$tmp/$name.expm.rc")" -eq 0 ] || st=1
    ok "$name: both spellings load and run" "$st"
    st=0; cmp -s "$tmp/$name.incipient.csv" "$tmp/$name.expm.csv" || st=1
    [ -s "$tmp/$name.expm.csv" ] || st=1
    ok "$name: incipient and expm CSVs are byte-identical" "$st"
    n=$(grep -c 'deprecated' "$tmp/$name.incipient.err" || true)
    st=0; [ "$n" -eq 1 ] || { echo "    $n notices" >&2; st=1; }
    grep 'deprecated' "$tmp/$name.incipient.err" | grep -q 'expm' || st=1
    grep 'deprecated' "$tmp/$name.incipient.err" | grep -q 'matrix exponential' || st=1
    ok "$name: incipient prints one notice naming expm" "$st"
    st=0; ! grep -q 'deprecated' "$tmp/$name.expm.err" || st=1
    ok "$name: expm prints no notice" "$st"
done

# Serialisation writes back the spelling the model was loaded with (ADR 0022 decision 4).
for m in incipient expm; do
    mkdir -p "$tmp/ser.$m"
    "$DUMP" "$tmp/ser.$m" "$tmp/decay_model.$m.json" > /dev/null 2>&1 || true
    st=0
    grep -q "\"method\":[[:space:]]*\"$m\"" "$tmp/ser.$m/decay_model.$m.model.json" 2>/dev/null || st=1
    ok "a model loaded as $m serialises as $m" "$st"
done

echo
if [ "$failures" -eq 0 ]; then echo "ALL PASS"; else echo "FAILURES PRESENT"; fi
echo "failures: $failures"
[ "$failures" -eq 0 ]
