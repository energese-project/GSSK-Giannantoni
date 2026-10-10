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

echo
if [ "$failures" -eq 0 ]; then echo "ALL PASS"; else echo "FAILURES PRESENT"; fi
echo "failures: $failures"
[ "$failures" -eq 0 ]
