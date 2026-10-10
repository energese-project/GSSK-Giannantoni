#!/bin/sh
# Self-test for the engine comparison page (docs/results/engine_comparison.md
# and its SVG plots) and the CLI flag it depends on.
#
# Verifies: FR-OUT-005 (T-OUT-05)
#
# Regenerating the page takes about two minutes (supply_chain_30 through the
# Giannantoni engine), so `make check-engine-comparison` runs once, in the
# quality gate. This runs on every toolchain and checks what can be checked
# fast: that the committed page is complete and well formed, that one of its
# numbers agrees with an independent tool, and that --full-precision prints
# the same trajectory as the default output, only with every digit.
set -u
cd "$(dirname "$0")/.." || exit 2
TMP=$(mktemp -d) || exit 2
trap 'rm -rf "$TMP"' EXIT
fail=0
check() { if eval "$2"; then echo "  ok   $1"; else echo "  FAIL $1"; fail=1; fi; }
PAGE=docs/results/engine_comparison.md

# 1. --full-precision: same rows, rounded the same, with 17 significant digits.
./bin/gssk examples/decay_model.json "$TMP/short.csv" > /dev/null 2>&1
./bin/gssk examples/decay_model.json "$TMP/full.csv" --full-precision > /dev/null 2>&1
check "--full-precision rounds to the default output, row for row" \
      'python3 -I - "$TMP/short.csv" "$TMP/full.csv" <<PY
import sys
a = [l.split(",") for l in open(sys.argv[1]).read().split()]
b = [l.split(",") for l in open(sys.argv[2]).read().split()]
assert a[0] == b[0] and len(a) == len(b) > 2
for ra, rb in zip(a[1:], b[1:]):
    assert "%.4f" % float(rb[0]) == ra[0]
    assert all("%.6f" % float(y) == x for x, y in zip(ra[1:], rb[1:]))
assert max(len(v.replace(".", "").lstrip("0")) for r in b[1:] for v in r[1:]) >= 15
PY'
check "the default output still matches its golden file" \
      './bin/csv_compare tests/expected/decay_model.csv "$TMP/short.csv" 2>/dev/null'

# 2. The page is complete: one plot per model, every plot well-formed SVG.
check "a plot for every kernel model (all but invalid_model)" \
      '( for m in examples/*.json; do n=$(basename $m .json); [ $n = invalid_model ] && continue; [ -f docs/results/plots/$n.svg ] || exit 1; done )'
check "a plot for every Giannantoni-format model" \
      '( for m in examples/giannantoni/*.json; do [ -f docs/results/plots/gia_$(basename $m .json).svg ] || exit 1; done )'
check "every plot parses as SVG, with a title and a description" \
      'python3 -I -c "
import glob, xml.etree.ElementTree as E
ns = \"{http://www.w3.org/2000/svg}\"
for p in glob.glob(\"docs/results/plots/*.svg\"):
    r = E.parse(p).getroot()
    assert r.tag == ns + \"svg\" and r.find(ns + \"title\").text and r.find(ns + \"desc\").text, p
"'
check "the page links every plot it has" \
      '[ $(grep -c "](plots/" $PAGE) -eq $(ls docs/results/plots/*.svg | wc -l) ]'

# 3. One number against an independent tool: decay_model's RK4 difference is
#    the bridge's own max_scaled_difference for its one storage, biomass.
want=$(./bin/gia_bridge examples/decay_model.json --method rk4 2>/dev/null |
       awk -F': ' '/^bridge.max_scaled_difference.biomass:/ { printf "%.2e", $2 }')
check "decay_model RK4 on the page ($want) is the bridge's figure" \
      'grep -q "^| \[\`decay_model\`\](#decay_model) | 100.0% | Giannantoni exp(At) | exact | rk4 | $want |" $PAGE'
check "decay_model expm on the page is exact to rounding" \
      'grep "^| \[\`decay_model\`\]" $PAGE | grep -q "| < 1e-13 |\$"'

[ $fail -eq 0 ] && echo "engine_compare: OK" || echo "engine_compare: FAILED"
exit $fail
