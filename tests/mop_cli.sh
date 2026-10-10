#!/bin/sh
# System tests for bin/giannantoni_sim (docs/requirements/vv-plan.md §1, level "System").
#
# Each check runs the CLI end to end and inspects what it wrote: the run report on stdout,
# the trajectory CSV, and the exit status. The oracle is the interface specification
# (docs/requirements/icd.md), never a previous run's output.
#
# POSIX sh + awk, so it behaves the same under macOS awk and Linux mawk.

set -eu
LC_ALL=C; export LC_ALL
cd "$(dirname "$0")/.."

SIM=${SIM:-bin/giannantoni_sim}
[ -x "$SIM" ] || { echo "mop_cli: $SIM not built" >&2; exit 1; }

tmp=$(mktemp -d "${TMPDIR:-/tmp}/mop_cli.XXXXXX")
trap 'rm -rf "$tmp"' EXIT INT TERM
failures=0

ok() {   # ok <what> <status: 0 pass>
    if [ "$2" -eq 0 ]; then printf '  %-62s PASS\n' "$1"
    else printf '  %-62s FAIL\n' "$1"; failures=$((failures + 1)); fi
}

run() {  # run <name> <model> : writes $tmp/<name>.{out,err,csv,json}, status in $tmp/<name>.rc
    set +e
    "$SIM" "$2" --steps 4 --csv "$tmp/$1.csv" --out "$tmp/$1.json" \
        >"$tmp/$1.out" 2>"$tmp/$1.err"
    echo $? > "$tmp/$1.rc"
    set -e
}

echo "=== giannantoni_sim system tests ==="

# ---------------------------------------------------------------------------
# Verifies: FR-OUT-001, BR-009 (T-OUT-02)
#
# Every output column, and every reported quantity, carries exactly one label from
# the closed set of FR-OUT-001. The labels are read from the report and matched
# against the CSV header the same run wrote, so a column added without a label, or a
# label left behind for a removed column, fails.
# ---------------------------------------------------------------------------
LABELS='implemented classical assumed proxy illustrative'

for model in examples/giannantoni/*.json; do
    name=$(basename "$model" .json)
    run "$name" "$model"
    ok "$name: exits 0" "$(cat "$tmp/$name.rc")"

    head -n 1 "$tmp/$name.csv" | tr ',' '\n' | sort > "$tmp/$name.cols"
    awk -F': ' '/^label\./ { sub(/^label\./, "", $1); print $1 }' "$tmp/$name.out" \
        | sort > "$tmp/$name.labelled"

    # One label line per column, and no column without one.
    st=0; [ -s "$tmp/$name.cols" ] || st=1
    for c in $(cat "$tmp/$name.cols"); do
        n=$(grep -cx "$c" "$tmp/$name.labelled" || true)
        [ "$n" -eq 1 ] || { echo "    column $c has $n label lines" >&2; st=1; }
    done
    ok "$name: one label line per CSV column" "$st"

    # Every label is one of the five.
    st=0
    awk -F': ' '/^label\./ { print $2 }' "$tmp/$name.out" > "$tmp/$name.values"
    [ -s "$tmp/$name.values" ] || st=1
    while IFS= read -r v; do
        case " $LABELS " in *" $v "*) ;; *) echo "    unknown label '$v'" >&2; st=1 ;; esac
    done < "$tmp/$name.values"
    ok "$name: every label is from the closed set" "$st"

    # The reported quantities that are not CSV columns are labelled too.
    st=0
    for q in ordinality maximum_ordinality closure harmony generative_step solution_drift; do
        grep -q "^label\.$q: " "$tmp/$name.out" || { echo "    no label for $q" >&2; st=1; }
    done
    ok "$name: reported quantities labelled" "$st"

    # E3 in PLAN §2: the trajectories are the matrix exponential, and may not be
    # labelled as an implementation of Giannantoni's calculus.
    st=0
    awk -F': ' '/^label\..*_Q: / && $2 != "classical" { bad = 1 }
                /^label\.closure: / && $2 != "proxy" { bad = 1 }
                /^label\.harmony: / && $2 != "assumed" { bad = 1 }
                END { exit bad }' "$tmp/$name.out" || st=1
    ok "$name: classical, proxy and assumed outputs say so" "$st"
done

# ---------------------------------------------------------------------------
# Verifies: IF-OUT-001 (T-OUT-01)
#
# The CSV header is exactly icd.md IF-OUT-001's: time; per node in seed order
# <id>_Q, <id>_Em, <id>_Tr, <id>_drift_proj; then conservation, emergy_excess. The
# node order is the seed file's, listed here by hand. Every value is a finite number.
# ---------------------------------------------------------------------------
want_header() {  # want_header <node ids...>
    printf 'time'
    for id in "$@"; do printf ',%s_Q,%s_Em,%s_Tr,%s_drift_proj' "$id" "$id" "$id" "$id"; done
    printf ',conservation,emergy_excess\n'
}
for spec in "input:source_1 interaction_1 store_1 consumer_1 heat_1" \
            "trophic_chain:sun producer herbivore decomposer heat" \
            "coproduction:sun soil nutrients growth wood seeds heat" \
            "harmonic_couples:sun a b c" \
            "closed_loop:source_1 interaction_1 store_1 consumer_1 recycler_1 heat_1"; do
    name=${spec%%:*}; ids=${spec#*:}
    run "$name" "examples/giannantoni/$name.json"
    # shellcheck disable=SC2086
    want_header $ids > "$tmp/$name.want"
    head -n 1 "$tmp/$name.csv" > "$tmp/$name.got"
    st=0; cmp -s "$tmp/$name.want" "$tmp/$name.got" || { echo "    header differs:" >&2; cat "$tmp/$name.got" >&2; st=1; }
    ok "$name: CSV header is exactly IF-OUT-001's" "$st"
    st=0
    awk -F, 'NR > 1 { for (i = 1; i <= NF; i++)
                          if ($i !~ /^-?[0-9]+(\.[0-9]*)?([eE][-+]?[0-9]+)?$/) { bad = 1; print "    row " NR ", field " i ": [" $i "]" > "/dev/stderr" } }
             END { exit bad }' "$tmp/$name.csv" || st=1
    ok "$name: every value is a finite number" "$st"
done

# ---------------------------------------------------------------------------
# Verifies: NFR-DET-001 (T-DET-01)
#
# The same binary on the same input writes byte-identical output: CSV, generated
# graph and report. Both runs use the same paths, since the report names them.
# ---------------------------------------------------------------------------
for model in examples/giannantoni/*.json; do
    name=$(basename "$model" .json)
    st=0
    for r in 1 2; do
        run det "$model"
        for x in out csv json; do cp "$tmp/det.$x" "$tmp/det$r.$x"; done
    done
    for x in out csv json; do
        cmp -s "$tmp/det1.$x" "$tmp/det2.$x" || { echo "    .$x differs between runs" >&2; st=1; }
    done
    ok "$name: two runs are byte-identical (report, CSV, graph)" "$st"
done

# ---------------------------------------------------------------------------
# Verifies: FR-ORD-002, FR-ORD-003, FR-ORD-004, IF-OUT-003 (T-ORD-03)
#
# The report carries ADR 0021's record and verdict, and closure only as a proxy.
# The values are ADR 0021 §3's table, hand-derived there: input.json is
# {2, 0, 0, 0, 1}, below, closure 0.000; closed_loop.json is {3, 3, 0, 0, 0}, at
# maximum, closure 1.000. The catalogue mutation "gate on closure" is
# closed_loop's old 0.800 with the sink counted.
# ---------------------------------------------------------------------------
# The three later seeds, by hand from each drawing:
#   trophic_chain    producer -> herbivore -> decomposer, producer -> decomposer: every couple one-way,
#                    no module, no replicate: {3, 0, 0, 0, 3}, no cycle.
#   coproduction     soil (energy) and nutrients (drawn control) both feed `growth`: (nutrients, soil)
#                    is 2; wood and seeds are its replicated products: 1/2; the other four unrelated.
#   harmonic_couples a -> b -> c -> a: all three couples 2/2, closure 1.
for spec in "input:{2, 0, 0, 0, 1}:no:0.000" "closed_loop:{3, 3, 0, 0, 0}:yes:1.000" \
            "trophic_chain:{3, 0, 0, 0, 3}:no:0.000" "coproduction:{4, 0, 1, 1, 4}:no:0.000" \
            "harmonic_couples:{3, 3, 0, 0, 0}:yes:1.000"; do
    name=${spec%%:*}; rest=${spec#*:}; rec=${rest%%:*}; rest=${rest#*:}
    max=${rest%%:*}; clo=${rest#*:}
    run "$name" "examples/giannantoni/$name.json"
    st=0
    grep -qxF "ordinality: $rec" "$tmp/$name.out" || { grep '^ordinality' "$tmp/$name.out" >&2 || true; st=1; }
    grep -qxF "maximum_ordinality: $max" "$tmp/$name.out" || st=1
    grep -qxF "closure (proxy): $clo" "$tmp/$name.out" || { grep '^closure' "$tmp/$name.out" >&2 || true; st=1; }
    ok "$name: ordinality $rec, maximum $max, closure (proxy) $clo" "$st"
done

# ---------------------------------------------------------------------------
# Verifies: FR-OUT-002, IF-CLI-001, BR-009 (T-OUT-03)
#
# icd.md IF-CLI-001: exit 0 on success, 1 on a load or validation error, 2 on a
# refusal, with the reason on stderr. FR-OUT-002: the reason names the feature and
# its source. The seeds are the T-ROB-01 corpus's; the expected status and the
# source each refusal must name are the interface's, by hand.
# ---------------------------------------------------------------------------
mop_run() {  # mop_run <corpus case> : $tmp/mop.{out,err,csv}, status in $tmp/mop.rc
    rm -f "$tmp/mop.csv"
    set +e
    "$SIM" "tests/mop_fuzz/$1.json" --steps 3 --csv "$tmp/mopt.csv" --out "$tmp/mop.json" \
        --mop-out "$tmp/mop.csv" >"$tmp/mop.out" 2>"$tmp/mop.err"
    echo $? > "$tmp/mop.rc"
    set -e
}
for spec in "ok_affine:0:" "ok_full:0:" \
            "ok_refused_domain:2:FR-MOP-002" "ok_network:2:FR-MOP-008" \
            "ok_refused_second:2:FR-MOP-005" "ok_refused_eqs_x11:2:X11" \
            "bad_k_negative:1:mop.k" "bad_couple_module_free_sink:1:mop.beta[0]" \
            "bad_truncated_json:1:" "ok_none:1:--mop-out"; do
    case_=${spec%%:*}; rest=${spec#*:}; want=${rest%%:*}; names=${rest#*:}
    mop_run "$case_"
    st=0; [ "$(cat "$tmp/mop.rc")" -eq "$want" ] || { echo "    exit $(cat "$tmp/mop.rc"), want $want" >&2; st=1; }
    if [ "$want" -ne 0 ]; then
        [ -s "$tmp/mop.err" ] || { echo "    nothing on stderr" >&2; st=1; }
        [ -z "$names" ] || grep -qF -- "$names" "$tmp/mop.err" || { echo "    stderr does not name $names:" >&2; cat "$tmp/mop.err" >&2; st=1; }
        [ ! -e "$tmp/mop.csv" ] || { echo "    a MOP CSV was written" >&2; st=1; }
    fi
    ok "$case_: exits $want${names:+, stderr names $names}" "$st"
done
# A refusal says "refused"; a load error does not.
mop_run ok_refused_domain; st=0; grep -q 'refused' "$tmp/mop.err" || st=1
mop_run bad_k_negative;    grep -q 'refused' "$tmp/mop.err" && st=1
ok "a refusal says refused on stderr; a load error does not" "$st"

# IF-OUT-002 through the CLI: the header, and the run report's label.
mop_run ok_full
st=0
[ "$(head -n 1 "$tmp/mop.csv")" = "time,a__b_re,a__b_im,R_H" ] || { head -n 1 "$tmp/mop.csv" >&2; st=1; }
[ "$(wc -l < "$tmp/mop.csv" | tr -d ' ')" -eq 5 ] || st=1
grep -qx 'label.mop_alpha: implemented' "$tmp/mop.out" || st=1
grep -qx 'label.second_equation: assumed' "$tmp/mop.out" || st=1
ok "--mop-out: IF-OUT-002's header, 1 + steps rows, labelled" "$st"

# ---------------------------------------------------------------------------
# Verifies: FR-HAR-004 (T-HAR-08)
#
# The constructor's matrix has one row per component (ADRs 0014, 0021): input.json has five
# nodes but two components (store_1, consumer_1), closed_loop.json three. It was once sized from
# every node, modules and habitat included, under a "components N" label.
# ---------------------------------------------------------------------------
for spec in "input:2" "closed_loop:3" "trophic_chain:3" "coproduction:4" "harmonic_couples:3"; do
    name=${spec%%:*}; want=${spec#*:}
    run "$name" "examples/giannantoni/$name.json"
    st=0
    grep -qE "^  components N +$want\$" "$tmp/$name.out" || { grep 'components N' "$tmp/$name.out" >&2 || true; st=1; }
    ok "$name: the harmony constructor has N = $want rows, one per component" "$st"
done

# ---------------------------------------------------------------------------
# Verifies: FR-OUT-003, IF-OUT-003, FR-HAR-002 (T-OUT-02)
#
# One `harmony.<construction>: <verdict>` line per evaluated construction, verbatim from the
# detector, and as docs/results/harmony_verdicts.md records for each seed. A "not evaluated" row
# means the report has no harmony line for that construction.
# ---------------------------------------------------------------------------
st=0; rows=0
grep '^| `examples/' docs/results/harmony_verdicts.md \
    | awk -F'|' '{ gsub(/[ `]/, "", $2); gsub(/^ +| +$/, "", $3); gsub(/^ +| +$/, "", $4); print $2 "|" $3 "|" $4 }' \
    > "$tmp/verdicts"
while IFS='|' read -r seed key verdict; do
    rows=$((rows + 1))
    run hv "$seed"
    case "$verdict" in
        "not evaluated")
            ! grep -q "^harmony\.$key: " "$tmp/hv.out" || { echo "    $seed: $key has a verdict, the table says not evaluated" >&2; st=1; } ;;
        *)
            grep -qx "harmony\.$key: $verdict" "$tmp/hv.out" || { echo "    $seed: harmony.$key is not $verdict" >&2; st=1; } ;;
    esac
done < "$tmp/verdicts"
[ "$rows" -ge 8 ] || { echo "    only $rows table rows read" >&2; st=1; }
ok "harmony verdicts match docs/results/harmony_verdicts.md ($rows rows)" "$st"
st=0
for v in $(sed -n 's/^harmony\.[a-z_]*: //p' "$tmp/hv.out"); do
    case "$v" in imposed|transported|present|absent) ;; *) echo "    verdict '$v'" >&2; st=1 ;; esac
done
ok "every harmony verdict is one of FR-HAR-002's four words" "$st"

# ---------------------------------------------------------------------------
# Verifies: IF-OUT-002, FR-HAR-001 (T-OUT-01)
#
# harmonic_couples.json's mop block gives beta_ac = -beta_ab, k = 1, so alpha_ac = -alpha_ab at every
# t, and with N = 3 the one root is -1: R_H = 0 (to rounding) for t > 0, and the cell is empty at
# t = 0 where alpha_12 = 0. By hand: alpha_ab(1) = 1 + 0.25/2 = 1.125.
# ---------------------------------------------------------------------------
set +e
"$SIM" examples/giannantoni/harmonic_couples.json --steps 4 --csv "$tmp/hc.csv" --out "$tmp/hc.json" \
    --mop-out "$tmp/hc.mop.csv" > "$tmp/hc.out" 2>&1
rc=$?
set -e
st=0
[ "$rc" -eq 0 ] || st=1
[ "$(head -n 1 "$tmp/hc.mop.csv")" = "time,a__b_re,a__b_im,a__c_re,a__c_im,R_H" ] || st=1
awk -F, 'NR == 2 && $6 != "" { bad = 1 }
         NR > 2 { r = $6 + 0; if ($6 == "" || r > 1e-12 || r < -1e-12) bad = 1 }
         END { if (NR != 6) bad = 1; exit bad }' "$tmp/hc.mop.csv" || st=1
awk -F, 'END { if ($2 != 1.125 || $4 != -1.125) exit 1 }' "$tmp/hc.mop.csv" || st=1
ok "harmonic_couples: R_H = 0 for t > 0, empty at t = 0; alpha_ab(1) = 1.125" "$st"

echo
if [ "$failures" -eq 0 ]; then echo "ALL PASS"; else echo "FAILURES PRESENT"; fi
echo "failures: $failures"
[ "$failures" -eq 0 ]
