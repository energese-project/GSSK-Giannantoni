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
# Verifies: FR-OUT-001 (T-OUT-02)
#
# Every output column, and every reported quantity, carries exactly one label from
# the closed set of FR-OUT-001. The labels are read from the report and matched
# against the CSV header the same run wrote, so a column added without a label, or a
# label left behind for a removed column, fails.
# ---------------------------------------------------------------------------
LABELS='implemented classical assumed proxy illustrative'

for model in examples/giannantoni/input.json examples/giannantoni/closed_loop.json; do
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
    for q in ordinality harmony generative_step solution_drift; do
        grep -q "^label\.$q: " "$tmp/$name.out" || { echo "    no label for $q" >&2; st=1; }
    done
    ok "$name: reported quantities labelled" "$st"

    # E3 in PLAN §2: the trajectories are the matrix exponential, and may not be
    # labelled as an implementation of Giannantoni's calculus.
    st=0
    awk -F': ' '/^label\..*_Q: / && $2 != "classical" { bad = 1 }
                /^label\.ordinality: / && $2 != "proxy" { bad = 1 }
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
for model in examples/giannantoni/input.json examples/giannantoni/closed_loop.json; do
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
# Verifies: FR-OUT-002, IF-CLI-001 (T-OUT-03)
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
[ "$(head -n 1 "$tmp/mop.csv")" = "time,a__b_re,a__b_im" ] || { head -n 1 "$tmp/mop.csv" >&2; st=1; }
[ "$(wc -l < "$tmp/mop.csv" | tr -d ' ')" -eq 5 ] || st=1
grep -qx 'label.mop_alpha: implemented' "$tmp/mop.out" || st=1
grep -qx 'label.second_equation: assumed' "$tmp/mop.out" || st=1
ok "--mop-out: IF-OUT-002's header, 1 + steps rows, labelled" "$st"

echo
if [ "$failures" -eq 0 ]; then echo "ALL PASS"; else echo "FAILURES PRESENT"; fi
echo "failures: $failures"
[ "$failures" -eq 0 ]
