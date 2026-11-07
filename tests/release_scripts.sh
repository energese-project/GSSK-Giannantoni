#!/usr/bin/env bash
# The scripts behind the automatic release on a merge to main:
#   scripts/next_version.sh   the version after a bump
#   scripts/release_level.sh  how far to bump, from the PR's labels or the commit message
#   scripts/auto_release.sh   the release commit (through scripts/release.sh)
#
# The oracle is semantic versioning and the CHANGELOG format scripts/release.sh
# writes, not the scripts' own output. auto_release.sh runs in a scratch git
# repository holding copies of the real version files, so nothing here touches
# this checkout.
set -euo pipefail
cd "$(dirname "$0")/.."
root=$PWD
fail=0
check() {
  local what=$1 got=$2 want=$3
  if [ "$got" = "$want" ]; then
    printf '  %-62s PASS\n' "$what"
  else
    printf '  %-62s FAIL (got %s, want %s)\n' "$what" "$got" "$want"
    fail=1
  fi
}
status() { "$@" > /dev/null 2>&1 && echo 0 || echo $?; }

echo "next_version.sh"
check "6.0.1 patch -> 6.0.2" "$(scripts/next_version.sh 6.0.1 patch)" 6.0.2
check "6.0.1 minor -> 6.1.0 (patch resets)" "$(scripts/next_version.sh 6.0.1 minor)" 6.1.0
check "6.3.9 major -> 7.0.0 (minor and patch reset)" "$(scripts/next_version.sh 6.3.9 major)" 7.0.0
check "9.9.9 patch -> 9.9.10 (no carry)" "$(scripts/next_version.sh 9.9.9 patch)" 9.9.10
check "refuses 6.0 (two parts)" "$(status scripts/next_version.sh 6.0 patch)" 2
check "refuses v6.0.1" "$(status scripts/next_version.sh v6.0.1 patch)" 2
check "refuses an unknown level" "$(status scripts/next_version.sh 6.0.1 huge)" 2

echo "release_level.sh"
level() { local msg=$1; shift; printf '%s\n' "$msg" | scripts/release_level.sh "$@"; }
check "no labels, plain message: patch" "$(level 'Merge pull request #1')" patch
check "release:minor" "$(level '' bug release:minor)" minor
check "release:major" "$(level '' release:major)" major
check "release:skip" "$(level '' release:skip docs)" skip
check "two bump labels take the larger" "$(level '' release:minor release:major release:patch)" major
check "a label beats the commit message" "$(level 'Release: major' release:patch)" patch
check "'Release: minor' in the message, no label" "$(level $'Title\n\nRelease: minor')" minor
check "'Release: skip' in the message, no label" "$(level $'Title\n\nRelease: skip')" skip
check "a mention of release: in prose is not a trailer" "$(level 'This fixes the release: notes')" patch
check "release:skip with release:minor fails" "$(status level '' release:skip release:minor)" 2
# Found running the workflow step: with a label deciding, the script exited
# without reading stdin, the `git log | release_level.sh` writer got SIGPIPE,
# and pipefail failed the release (exit 141). A long message makes it certain.
check "reads a long message to the end under pipefail" \
  "$(set -o pipefail; head -c 300000 /dev/zero | tr '\0' 'x' | fold -w 80 | scripts/release_level.sh release:skip > /dev/null; echo $?)" 0

echo "auto_release.sh"
scratch() {
  local dir
  dir=$(mktemp -d)
  mkdir -p "$dir/include" "$dir/docs" "$dir/scripts"
  cp include/gssk.h "$dir/include/"
  cp package.json CITATION.cff "$dir/"
  cp scripts/release.sh scripts/next_version.sh scripts/auto_release.sh scripts/check_version_sync.py "$dir/scripts/"
  cp gssk.schema.json "$dir/" 2> /dev/null || true
  # A CHANGELOG with the given [Unreleased] body above a released section.
  printf '# Changelog\n\n## [Unreleased]\n%s\n---\n\n## [1.0.0] — 2020-01-01\n\n- old\n' "$1" > "$dir/docs/CHANGELOG.md"
  (cd "$dir" && git init -q && git -c user.name=t -c user.email=t@t add -A && git -c user.name=t -c user.email=t@t commit -qm base)
  echo "$dir"
}
current=$(sed -n 's/^#define GSK_VERSION_STRING "\(.*\)"$/\1/p' include/gssk.h)
expected=$(scripts/next_version.sh "$current" minor)

dir=$(scratch $'\n### Added\n\n- **A feature.** It does things.\n\n### Fixed\n')
new=$(cd "$dir" && GIT_AUTHOR_NAME=t GIT_AUTHOR_EMAIL=t@t GIT_COMMITTER_NAME=t GIT_COMMITTER_EMAIL=t@t scripts/auto_release.sh minor 'Unused note (#9)')
check "prints the new version" "$new" "$expected"
check "gssk.h, package.json and CITATION.cff agree" "$(cd "$dir" && python3 scripts/check_version_sync.py > /dev/null 2>&1 && echo ok)" ok
log=$(git -C "$dir" log -1 --format=%s)
check "commits 'chore: release v$expected'" "$log" "chore: release v$expected"
cl=$(cat "$dir/docs/CHANGELOG.md")
check "dates the release heading" "$(grep -c "^## \[$expected\] — $(date +%Y-%m-%d)$" <<< "$cl")" 1
check "keeps the entry" "$(grep -c '^- \*\*A feature' <<< "$cl")" 1
release=$(sed -n "/^## \[$expected\]/,/^---$/p" <<< "$cl")
check "drops the empty ### Fixed heading" "$(grep -c '^### Fixed' <<< "$release")" 0
check "does not add the note when there are entries" "$(grep -c 'Unused note' <<< "$cl")" 0
check "opens a fresh [Unreleased] above it" "$(grep -n '^## \[' <<< "$cl" | head -1 | cut -d: -f2-)" "## [Unreleased]"
rm -rf "$dir"

dir=$(scratch $'\n### Added\n')
new=$(cd "$dir" && GIT_AUTHOR_NAME=t GIT_AUTHOR_EMAIL=t@t GIT_COMMITTER_NAME=t GIT_COMMITTER_EMAIL=t@t scripts/auto_release.sh patch 'Tidy the docs (#12)')
release=$(sed -n "/^## \[$new\]/,/^---$/p" "$dir/docs/CHANGELOG.md")
check "with no entries, the PR's title becomes one" "$(grep -c '^- Tidy the docs (#12)$' <<< "$release")" 1
check "under ### Changed" "$(grep -c '^### Changed$' <<< "$release")" 1
check "with no entries and no note, it refuses" "$(cd "$dir" && git checkout -q HEAD~1 -- . && status scripts/auto_release.sh patch)" 1
rm -rf "$dir"

[ "$fail" -eq 0 ] && echo "OK (0 failures)" || { echo "FAILED"; exit 1; }
