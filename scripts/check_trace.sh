#!/bin/sh
# Verifies: NFR-TRC-001 (T-TRC-01)
#
# Requirements traceability check for the Giannantoni kernel. See
# docs/requirements/README.md, "The trace check". Fails when:
#   1. a requirement ID is duplicated or malformed;
#   2. a requirement names no verification ID, or one not in the catalogue;
#   3. a catalogue entry verifies nothing, or names a requirement that does not exist;
#   4. a requirement is `implemented` but no tag under tests/ or scripts/ names it;
#   5. a tag names a requirement, or a verification ID, that does not exist.
#
# POSIX sh + awk + grep -E: it must behave the same under macOS awk and Linux mawk.

set -eu
LC_ALL=C; export LC_ALL          # sort and comm must agree on collation
cd "$(dirname "$0")/.."

REQ_DIR=docs/requirements
REQ_FILES="$REQ_DIR/stakeholder.md $REQ_DIR/srs.md $REQ_DIR/icd.md"
CATALOGUE="$REQ_DIR/vv-plan.md"
TAG_DIRS="tests scripts"

REQ_RE='(BR|FR|NFR|IF)(-[A-Z]+)?-[0-9][0-9][0-9]'
VER_RE='((T|INS|ANA|DEM)-[A-Z]+|VAL)-[0-9][0-9]'

tmp=$(mktemp -d "${TMPDIR:-/tmp}/check_trace.XXXXXX")
trap 'rm -rf "$tmp"' EXIT INT TERM
fail=0
err() { printf 'check-trace: %s\n' "$*" >&2; fail=1; }

for f in $REQ_FILES $CATALOGUE; do
    [ -f "$f" ] || { err "missing $f"; exit 1; }
done

# ---- Requirements: one block per "### ID — title" heading ---------------------
# Emits:  REQ <id> <status>      and  VER <id> <verification-id>
#         BAD <file>:<line> <text> for a heading that looks like a requirement but is malformed
awk -v req_re="^$REQ_RE\$" -v ver_re="$VER_RE" '
    function flush() { if (cur != "") printf "REQ %s %s\n", cur, (status == "" ? "none" : status) }
    /^### / {
        flush(); cur = ""; status = ""
        id = $2
        if (id ~ /^(BR|FR|NFR|IF)/) {
            if (id ~ req_re) cur = id
            else printf "BAD %s:%d %s\n", FILENAME, FNR, $0
        }
        next
    }
    cur != "" && index($0, "**Verification:**") {
        s = substr($0, index($0, "**Verification:**") + 17)
        cut = index(s, "**Priority:**"); if (cut) s = substr(s, 1, cut - 1)
        while (match(s, ver_re)) {
            printf "VER %s %s\n", cur, substr(s, RSTART, RLENGTH)
            s = substr(s, RSTART + RLENGTH)
        }
        if (match($0, /\*\*Status:\*\* [a-z]+/)) status = substr($0, RSTART + 12, RLENGTH - 12)
    }
    END { flush() }
' $REQ_FILES > "$tmp/req"

grep '^BAD ' "$tmp/req" | while IFS= read -r line; do
    printf 'check-trace: malformed requirement heading: %s\n' "${line#BAD }" >&2
done
grep -q '^BAD ' "$tmp/req" && fail=1

awk '$1 == "REQ" { print $2 }' "$tmp/req" | sort > "$tmp/req_ids"
awk '$1 == "REQ" && $3 == "implemented" { print $2 }' "$tmp/req" | sort > "$tmp/req_impl"

# Rule 1: duplicates.
for id in $(uniq -d "$tmp/req_ids"); do err "duplicate requirement ID $id"; done
sort -u "$tmp/req_ids" > "$tmp/req_set"

# Status must be one of the three.
awk '$1 == "REQ" && $3 !~ /^(planned|implemented|withdrawn)$/ { print $2, $3 }' "$tmp/req" |
while read -r id st; do
    printf 'check-trace: %s has status "%s" (expected planned, implemented or withdrawn)\n' "$id" "$st" >&2
done
awk '$1 == "REQ" && $3 !~ /^(planned|implemented|withdrawn)$/' "$tmp/req" | grep -q . && fail=1

# ---- Catalogue: table rows whose first cell is a verification ID ---------------
# Emits:  CAT <ver-id>   and   CATV <ver-id> <req-id>
awk -F'|' -v ver_re="^$VER_RE\$" -v req_re="$REQ_RE" '
    /^\|/ {
        id = $2; gsub(/^[ \t]+|[ \t]+$/, "", id)
        if (id !~ ver_re) next
        print "CAT " id
        s = $3
        while (match(s, req_re)) {
            printf "CATV %s %s\n", id, substr(s, RSTART, RLENGTH)
            s = substr(s, RSTART + RLENGTH)
        }
    }
' "$CATALOGUE" > "$tmp/cat"
awk '$1 == "CAT" { print $2 }' "$tmp/cat" | sort > "$tmp/cat_ids"
for id in $(uniq -d "$tmp/cat_ids"); do err "duplicate catalogue ID $id"; done
sort -u "$tmp/cat_ids" > "$tmp/cat_set"

# Rule 2: every requirement names at least one catalogued verification ID.
awk '$1 == "VER" { print $2 }' "$tmp/req" | sort -u > "$tmp/req_with_ver"
for id in $(comm -23 "$tmp/req_set" "$tmp/req_with_ver"); do
    err "$id names no verification ID"
done
awk '$1 == "VER" { print $3, $2 }' "$tmp/req" | sort -u | while read -r v id; do
    grep -qx "$v" "$tmp/cat_set" || printf 'check-trace: %s names %s, which is not in the catalogue\n' "$id" "$v" >&2
done
awk '$1 == "VER" { print $3 }' "$tmp/req" | sort -u > "$tmp/ver_named"
[ -z "$(comm -23 "$tmp/ver_named" "$tmp/cat_set")" ] || fail=1

# Rule 3: every catalogue entry verifies an existing requirement, and at least one.
awk '$1 == "CATV" { print $2 }' "$tmp/cat" | sort -u > "$tmp/cat_with_req"
for v in $(comm -23 "$tmp/cat_set" "$tmp/cat_with_req"); do
    err "catalogue entry $v verifies no requirement"
done
awk '$1 == "CATV" { print $3 }' "$tmp/cat" | sort -u > "$tmp/cat_reqs"
for id in $(comm -23 "$tmp/cat_reqs" "$tmp/req_set"); do
    err "catalogue names $id, which is not a requirement"
done

# ---- Tags in tests/ and scripts/ ---------------------------------------------
# A tag is the text after "Verifies:" on one line.
grep -rhoE 'Verifies:[^*]*' $TAG_DIRS 2>/dev/null > "$tmp/tags" || true
grep -oE "$REQ_RE" "$tmp/tags" | sort -u > "$tmp/tag_reqs" || true
grep -oE "$VER_RE" "$tmp/tags" | sort -u > "$tmp/tag_vers" || true

# Rule 5: tags name only existing requirements and catalogue IDs.
for id in $(comm -23 "$tmp/tag_reqs" "$tmp/req_set"); do
    err "a Verifies: tag names $id, which is not a requirement"
done
for v in $(comm -23 "$tmp/tag_vers" "$tmp/cat_set"); do
    err "a Verifies: tag names $v, which is not in the catalogue"
done

# Rule 4: implemented requirements are tagged somewhere.
for id in $(comm -23 "$tmp/req_impl" "$tmp/tag_reqs"); do
    err "$id is implemented but no Verifies: tag under $TAG_DIRS names it"
done

n_req=$(wc -l < "$tmp/req_set" | tr -d ' ')
n_impl=$(wc -l < "$tmp/req_impl" | tr -d ' ')
n_cat=$(wc -l < "$tmp/cat_set" | tr -d ' ')
if [ "$fail" -ne 0 ]; then
    echo "check-trace: FAILED" >&2
    exit 1
fi
echo "check-trace: OK — $n_req requirements ($n_impl implemented), $n_cat verification entries"
