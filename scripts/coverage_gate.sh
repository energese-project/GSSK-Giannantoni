#!/bin/sh
# Line-coverage gate over gcov's report (NFR-COV-001, PLAN.md G2, ADR 0018 rule 2).
#
#   coverage_gate.sh <min-percent> <gcov-report>
#
# <gcov-report> is the stdout of `gcov -n` for the units to be gated: pairs of
#     File 'src/engine.c'
#     Lines executed:90.28% of 1790
# Lines executed are summed across files (each file's percentage times its line count), so the gate
# is on the units together, as NFR-COV-001 states it. It FAILS, rather than passing, when:
#   - the report has no File/Lines pair it can parse (the old gate printed OK here: PLAN §1 B3);
#   - a File line has no Lines line after it, or a Lines line is not a number;
#   - the total is below <min-percent>.
# The parse is POSIX awk; the old one used `grep -oP`, which macOS grep does not have.

set -eu
LC_ALL=C; export LC_ALL

[ $# -eq 2 ] || { echo "usage: $0 <min-percent> <gcov-report>" >&2; exit 2; }
min=$1; report=$2
[ -r "$report" ] || { echo "coverage-gate: cannot read $report" >&2; exit 1; }

awk -v min="$min" '
    function bad(msg) { printf "coverage-gate: unparsable report: %s\n", msg > "/dev/stderr"; broken = 1 }
    /^File / {
        if (file != "") bad("no Lines line for " file)
        file = $0; sub(/^File \047?/, "", file); sub(/\047$/, "", file)
        next
    }
    /^Lines executed:/ {
        if (file == "") { bad("Lines line with no File line: " $0); next }
        s = $0; sub(/^Lines executed:/, "", s)
        pct = s; sub(/%.*/, "", pct)
        n = s; sub(/^[^ ]* of /, "", n)
        if (pct !~ /^[0-9]+(\.[0-9]+)?$/ || n !~ /^[0-9]+$/) { bad("not a number: " $0); file = ""; next }
        hit += pct * n / 100; total += n; files++
        printf "  %-40s %7.2f%% of %d\n", file, pct, n
        file = ""
        next
    }
    /^No executable lines/ { file = ""; next }
    END {
        if (file != "") bad("no Lines line for " file)
        if (files == 0) bad("no File/Lines pair found")
        if (broken) { print "coverage-gate: FAILED (report could not be read; a gate that cannot read its input does not pass)"; exit 1 }
        pct = 100 * hit / total
        printf "coverage-gate: %.2f%% of %d lines in %d unit(s) (gate: %s%%)\n", pct, total, files, min
        if (pct + 1e-9 < min) { print "coverage-gate: FAILED"; exit 1 }
        print "coverage-gate: OK"
    }
' "$report"
