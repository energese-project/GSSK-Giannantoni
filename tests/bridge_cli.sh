#!/bin/sh
# System tests for bin/gia_bridge: one GSSK-schema model through both engines.
#
# Verifies: FR-BRG-001 (T-BRG-01)
# Verifies: FR-BRG-002 (T-BRG-02)
#
# The oracle is closed form, never a previous run. examples/decay_model.json is one linear pathway,
# dQ/dt = -k Q, with k = 0.05, dt = 0.5, t_end = 20, Q0 = 100, so h = k dt = 0.025 and n = 40 steps:
#   the Giannantoni engine (exp(A t)) must give   Q0 e^{-k t}
#   the kernel's Euler gives exactly              Q0 (1 - h)^n
#   the kernel's RK4 gives exactly                Q0 R^n,  R = 1 - h + h^2/2 - h^3/6 + h^4/24
# so the bridge's reported difference at t_end is fixed by arithmetic, and must match it. A model
# the projection carries only in part is refused (exit 2): its difference would be the missing
# pathways, not integration error (TODO.md, the bridge item). POSIX sh + awk.

set -eu
LC_ALL=C; export LC_ALL
cd "$(dirname "$0")/.."

BRIDGE=${BRIDGE:-bin/gia_bridge}
[ -x "$BRIDGE" ] || { echo "bridge_cli: $BRIDGE not built" >&2; exit 1; }
tmp=$(mktemp -d "${TMPDIR:-/tmp}/bridge_cli.XXXXXX")
trap 'rm -rf "$tmp"' EXIT INT TERM
failures=0
ok() {
    if [ "$2" -eq 0 ]; then printf '  %-62s PASS\n' "$1"
    else printf '  %-62s FAIL\n' "$1"; failures=$((failures + 1)); fi
}
brun() {  # brun <name> <args...> : $tmp/<name>.{out,err,csv}, status in $tmp/<name>.rc
    n=$1; shift
    set +e
    "$BRIDGE" "$@" --csv "$tmp/$n.csv" > "$tmp/$n.out" 2> "$tmp/$n.err"
    echo $? > "$tmp/$n.rc"
    set -e
}
# value <file> <key>: the value of a `key: value` report line.
value() { awk -F': ' -v k="$2" '$1 == k { print $2; exit }' "$1"; }

echo "=== gia_bridge system tests ==="

for m in euler rk4 expm; do
    brun "$m" examples/decay_model.json --method "$m"
    st=0; [ "$(cat "$tmp/$m.rc")" -eq 0 ] || { cat "$tmp/$m.err" >&2; st=1; }
    ok "decay_model, kernel $m: exits 0" "$st"
    got=$(value "$tmp/$m.out" "bridge.max_rel_difference.biomass")
    st=0
    awk -v m="$m" -v got="$got" 'BEGIN {
        h = 0.025; n = 40; ex = exp(-1)
        if (m == "euler") k = (1 - h) ^ n
        else if (m == "rk4") k = (1 - h + h*h/2 - h*h*h/6 + h*h*h*h/24) ^ n
        else k = ex
        want = (k > ex ? k - ex : ex - k) / ex
        if (m == "expm") { if (!(got + 0 <= 1e-9)) exit 1; exit 0 }
        # Each difference is a difference of two numbers near e^-1 = 0.37, so it is good to
        # ~1e-16/want relative at best: the Euler difference (1.3e-2) to 1e-9, RK4 (3.3e-9) to 1e-5.
        tol = (m == "euler") ? 1e-9 : 1e-5
        d = got - want; if (d < 0) d = -d
        if (!(got != "" && d <= tol * want)) exit 1 }' || { echo "    $m: got $got" >&2; st=1; }
    ok "decay_model, kernel $m: max rel difference = the closed form's" "$st"
done

# The last CSV row: t = 20, the Giannantoni column is 100 e^{-1}, the kernel's Euler column 100 (1-h)^40.
st=0
head -n 1 "$tmp/euler.csv" | grep -qx 'time,biomass_kernel,biomass_gia,biomass_diff,environment_kernel,environment_gia,environment_diff' \
    || { head -n 1 "$tmp/euler.csv" >&2; st=1; }
tail -n 1 "$tmp/euler.csv" | awk -F, '{
    ex = 100 * exp(-1); eu = 100 * (1 - 0.025) ^ 40
    if ($1 + 0 != 20) exit 1
    if ((($3 - ex) ^ 2) ^ 0.5 > 1e-10 * ex) exit 1
    if ((($2 - eu) ^ 2) ^ 0.5 > 1e-10 * eu) exit 1
    if (((($4 - ($2 - $3)) ^ 2) ^ 0.5) > 1e-12) exit 1 }' || { tail -n 1 "$tmp/euler.csv" >&2; st=1; }
ok "CSV: header by node; t = 20 row is 100(1-h)^40, 100 e^-1, their diff" "$st"
[ "$(wc -l < "$tmp/euler.csv" | tr -d ' ')" -eq 42 ] && st=0 || st=1
ok "CSV: one row per kernel step, t = 0 .. 20 (41 rows + header)" "$st"

# Honest labels: both trajectories are classical, so the difference is integration error.
st=0
grep -qx 'label.bridge_difference: classical' "$tmp/euler.out" || st=1
grep -q 'not incipient drift' "$tmp/euler.out" || st=1
[ "$(value "$tmp/euler.out" "bridge.kernel_method")" = "euler" ] || st=1
[ "$(value "$tmp/euler.out" "bridge.coverage")" = "100.0%" ] || st=1
ok "report: coverage 100.0%, kernel method, difference labelled classical" "$st"

# Three stores joined by two reversible (diffusive) pathways: a linear network, which the
# Giannantoni engine solves exactly and the kernel's Padé expm per step to far below 1e-6.
brun diff examples/diffusion_model.json --method expm
st=0; [ "$(cat "$tmp/diff.rc")" -eq 0 ] || st=1
got=$(value "$tmp/diff.out" "bridge.max_rel_difference")
awk -v g="$got" 'BEGIN { exit !(g != "" && g + 0 <= 1e-6) }' || { echo "    got $got" >&2; st=1; }
ok "diffusion_model, kernel expm: the engines agree to 1e-6" "$st"

# A source driven by a sine (ADR 0006), projected into this engine's waveform vocabulary on its own
# clock: linear in the state, so under RK4 the engines agree to 1e-6, the forced source itself
# exactly. forced_source_model adds a square and a step on its edges, which no generator carries:
# refused, each named.
brun sine tests/fixtures/bridge_sine_source.json --method rk4
st=0; [ "$(cat "$tmp/sine.rc")" -eq 0 ] || st=1
got=$(value "$tmp/sine.out" "bridge.max_rel_difference")
awk -v g="$got" 'BEGIN { exit !(g != "" && g + 0 <= 1e-6) }' || { echo "    got $got" >&2; st=1; }
got=$(value "$tmp/sine.out" "bridge.max_rel_difference.sun")
awk -v g="$got" 'BEGIN { exit !(g != "" && g + 0 <= 1e-12) }' || { echo "    sun $got" >&2; st=1; }
ok "a sine-forced source, kernel rk4: the engines agree to 1e-6" "$st"
brun forced examples/forced_source_model.json
st=0; [ "$(cat "$tmp/forced.rc")" -eq 2 ] || st=1
grep -q "edge 'cropping': forcing not carried .*square" "$tmp/forced.out" || st=1
grep -q "edge 'respiration': forcing not carried .*step" "$tmp/forced.out" || st=1
ok "forced_source_model: square and step forcing named, refused (exit 2)" "$st"

# Refusals and load errors.
brun partial examples/atwood_model.json
st=0; [ "$(cat "$tmp/partial.rc")" -eq 2 ] || st=1
grep -q 'coverage' "$tmp/partial.err" || st=1
[ ! -s "$tmp/partial.csv" ] || st=1
ok "a partly projected model (atwood): exit 2, names coverage, no CSV" "$st"
brun missing /nonexistent.json
[ "$(cat "$tmp/missing.rc")" -eq 1 ] && st=0 || st=1
ok "an unreadable model: exit 1" "$st"
brun badm examples/decay_model.json --method leapfrog
[ "$(cat "$tmp/badm.rc")" -eq 1 ] && grep -q 'method' "$tmp/badm.err" && st=0 || st=1
ok "an unknown --method: exit 1, naming it" "$st"

echo
if [ "$failures" -eq 0 ]; then echo "ALL PASS"; else echo "FAILURES PRESENT"; fi
echo "failures: $failures"
[ "$failures" -eq 0 ]
