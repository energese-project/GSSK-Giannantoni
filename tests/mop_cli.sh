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
    for q in ordinality harmony generative_step duet; do
        grep -q "^label\.$q: " "$tmp/$name.out" || { echo "    no label for $q" >&2; st=1; }
    done
    ok "$name: reported quantities labelled" "$st"

    # E3/E4 in PLAN §2: the trajectories are the matrix exponential, and the per-node
    # calculus columns apply the drift identity to an invented phi. Neither may be
    # labelled as an implementation of Giannantoni's calculus.
    st=0
    awk -F': ' '/^label\..*_Q: / && $2 != "classical" { bad = 1 }
                /^label\..*_(idc|tdc|drift): / && $2 != "illustrative" { bad = 1 }
                /^label\.ordinality: / && $2 != "proxy" { bad = 1 }
                /^label\.harmony: / && $2 != "assumed" { bad = 1 }
                END { exit bad }' "$tmp/$name.out" || st=1
    ok "$name: classical, proxy and assumed outputs say so" "$st"
done

echo
if [ "$failures" -eq 0 ]; then echo "ALL PASS"; else echo "FAILURES PRESENT"; fi
echo "failures: $failures"
[ "$failures" -eq 0 ]
