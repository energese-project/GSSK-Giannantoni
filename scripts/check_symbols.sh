#!/bin/sh
# Structural checks on the Giannantoni library units (docs/requirements/vv-plan.md §7.8).
#
#   T-REE-02   no writable data symbols: no file-scope or static mutable state (NFR-REE-001)
#   T-ERR-01   no exit, abort or assert-to-abort in a library object (NFR-ERR-001)
#   INS-SEP-01 no Giannantoni unit includes gssk.h: the two engines stay separate (NFR-SEP-001)
#
# Verifies: NFR-REE-001 (T-REE-02)
# Verifies: NFR-ERR-001 (T-ERR-01)
# Verifies: NFR-SEP-001 (INS-SEP-01)
# Verifies: FR-HAR-003 (T-HAR-01)
#
#   T-HAR-01   no solver or detector unit calls the harmony constructor: no undefined reference
#              to gia_harmony_assume_* in mop, mop_seed, relational, idc or engine (FR-HAR-003).
#              bin/test_mop_emergence, linked without harmony.o, is the link-time half.
#
# The objects are read from lib/, so build first (`make check-symbols` does). The unit list is
# fixed here rather than globbed so a new unit is added on purpose; a listed unit whose source does
# not exist yet is skipped by name, so W2-W8 can land in any order.
#
# POSIX sh + awk + nm: GNU nm's System V format on Linux, `nm -m` on macOS, each of which names
# the section a symbol lives in.

set -eu
LC_ALL=C; export LC_ALL
cd "$(dirname "$0")/.."

UNITS="engine validation projection idc mop relational mop_seed harmony"
HEADERS="include/engine.h include/idc.h include/mop.h include/mop_seed.h include/relational.h include/gia_status.h"
LIB=${LIB_DIR:-lib}

tmpf=$(mktemp "${TMPDIR:-/tmp}/check_symbols.XXXXXX")
trap 'rm -f "$tmpf"' EXIT INT TERM

# writable_symbols <obj>: one "<name> (<section>)" line per symbol in a writable section.
writable_symbols() {
    if [ "$(uname -s)" = Darwin ]; then
        # nm -m: "<addr> (__DATA,__bss) non-external _sort_model"
        nm -m "$1" | awk '/\(__DATA,__(data|bss|common)\)/ { print $NF, $2 }'
    else
        # GNU nm, System V format: name|value|class|type|size|line|section
        nm --format=sysv "$1" | awk -F'|' '
            { c = $3; gsub(/ /, "", c); sec = $7; gsub(/ /, "", sec); n = $1; gsub(/ /, "", n) }
            c ~ /^[bBdD]$/ && sec !~ /^\.data\.rel\.ro/ { printf "%s (%s)\n", n, sec }'
    fi
}

fail=0
err() { printf 'check-symbols: %s\n' "$*" >&2; fail=1; }
checked=""

for u in $UNITS; do
    [ -f "src/$u.c" ] || continue
    obj="$LIB/$u.o"
    [ -f "$obj" ] || { err "$obj not built (run make all first)"; continue; }
    checked="$checked $u"

    # T-REE-02: writable data. A symbol in a writable data or bss section is mutable state.
    # Read-only data with relocations (a `static const char *const tbl[]`) is in
    # .data.rel.ro on ELF and __DATA,__const on Mach-O, which nm also letters 'd', so the
    # section, not the letter, decides.
    writable_symbols "$obj" > "$tmpf" || true
    while IFS= read -r line; do
        [ -n "$line" ] && err "$obj: writable data symbol $line -- mutable state is shared across models and threads (NFR-REE-001)"
    done < "$tmpf"

    # T-ERR-01: process exits. The names cover glibc, musl and Darwin.
    for s in $(nm -u "$obj" | awk '{ print $NF }'); do
        case "${s#_}" in
            exit|_exit|_Exit|abort|__assert_fail|_assert_rtn|_assert|quick_exit)
                err "$obj calls $s -- a library unit returns a status, it never ends the process (NFR-ERR-001)" ;;
        esac
    done

    # T-HAR-01: the constructor is called by no solver and no detector.
    case " mop mop_seed relational idc engine " in *" $u "*)
        for s in $(nm -u "$obj" | awk '{ print $NF }'); do
            case "${s#_}" in gia_harmony_assume_*)
                err "$obj calls $s -- harmony is detected, never constructed, in a solver (FR-HAR-003)" ;;
            esac
        done ;;
    esac

    # INS-SEP-01: the kernel's header.
    if grep -nE '^[[:space:]]*#[[:space:]]*include[[:space:]]*[<"]gssk\.h[>"]' "src/$u.c" >/dev/null; then
        err "src/$u.c includes gssk.h -- the two engines share no headers (ADR 0011, NFR-SEP-001)"
    fi
done
for h in $HEADERS; do
    [ -f "$h" ] || continue
    if grep -nE '^[[:space:]]*#[[:space:]]*include[[:space:]]*[<"]gssk\.h[>"]' "$h" >/dev/null; then
        err "$h includes gssk.h (ADR 0011, NFR-SEP-001)"
    fi
done

if [ "$fail" -ne 0 ]; then echo "check-symbols: FAILED" >&2; exit 1; fi
echo "check-symbols: OK --$checked"
