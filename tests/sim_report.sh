#!/bin/sh
# Self-test for scripts/sim_report.py, the simulation report CI writes to the
# job summary. A report that always says "matches golden" would be worse than
# none, so this proves each verdict can fire.
#
# Verifies: FR-OUT-004 (T-OUT-04)
#
# The oracle is a perturbation of known size: decay_model's golden file is
# copied and one value moved by a chosen amount, so the deviation the report
# must print is fixed in advance.
set -u
cd "$(dirname "$0")/.." || exit 2
REPORT="python3 -I scripts/sim_report.py"
TMP=$(mktemp -d) || exit 2
trap 'rm -rf "$TMP"' EXIT
fail=0
check() { if eval "$2"; then echo "  ok   $1"; else echo "  FAIL $1"; fail=1; fi; }

# Perturbs row 2 (t = first step), column 2 (biomass) of decay_model by $1.
perturbed() {
    rm -rf "$TMP/exp" && cp -R tests/expected "$TMP/exp"
    awk -F, -v OFS=, -v d="$1" 'NR == 3 { $2 = sprintf("%.10f", $2 + d) } { print }' \
        tests/expected/decay_model.csv > "$TMP/exp/decay_model.csv"
}

$REPORT > "$TMP/clean.md"; rc=$?
check "the committed goldens: exit 0" '[ $rc -eq 0 ]'
check "the committed goldens: none failed" 'grep -q "\*\*0\*\* failed" "$TMP/clean.md"'
check "every model has a row" \
      '[ $(grep -c -E "^\| \`.*(✅ matches golden|❌|⏭️ skipped)" "$TMP/clean.md") -eq $(ls examples/*.json | wc -l) ]'
check "the allowlisted model is reported skipped, with its reason" \
      'grep -q "invalid_model.*skipped: deliberately malformed" "$TMP/clean.md"'
# decay is monotone: its sparkline must fall from the top bar to the bottom one.
check "decay_model's biomass sparkline falls from █ to ▁" \
      'grep "^| \`decay_model\` " "$TMP/clean.md" | grep "biomass █" | grep -q "▁ |\$"'

perturbed 0.001
$REPORT --expected "$TMP/exp" > "$TMP/off.md"; rc=$?
check "a golden off by 1e-3: exit 1" '[ $rc -eq 1 ]'
check "a golden off by 1e-3: decay_model differs, deviation 1.0e-03" \
      'grep -q "^| \`decay_model\` | rk4 | 40 | 1.0e-03 | ❌ differs from golden" "$TMP/off.md"'
check "a golden off by 1e-3: the others still match" \
      'grep -q "\*\*1\*\* failed" "$TMP/off.md"'

perturbed 0.0000005
$REPORT --expected "$TMP/exp" > "$TMP/near.md"; rc=$?
check "a golden off by 5e-7 (inside 1e-6): exit 0, deviation shown" \
      '[ $rc -eq 0 ] && grep -q "^| \`decay_model\` | rk4 | 40 | 5.0e-07 | ✅ matches golden" "$TMP/near.md"'

rm -rf "$TMP/exp" && cp -R tests/expected "$TMP/exp" && rm "$TMP/exp/simple_model.csv"
$REPORT --expected "$TMP/exp" > "$TMP/gone.md"; rc=$?
check "a missing golden: exit 1, named" \
      '[ $rc -eq 1 ] && grep -q "^| \`simple_model\` .*❌ no golden file" "$TMP/gone.md"'

[ $fail -eq 0 ] && echo "sim_report: OK" || echo "sim_report: FAILED"
exit $fail
