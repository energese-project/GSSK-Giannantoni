#!/usr/bin/env bash
# scripts/auto_release.sh — cut the release a merge to main makes.
#
# Usage: scripts/auto_release.sh <patch|minor|major> [changelog note]
#
# Bumps the current version (include/gssk.h) by <level> and commits it with
# scripts/release.sh, so the version files, CITATION.cff and the dated
# CHANGELOG heading are exactly what a hand-cut release has. Before that it
# tidies `## [Unreleased]`: headings with no entries are dropped, and if no
# entry is left, the note (the PR's title and number) becomes one, so a
# release never ships an empty section. Prints the new version.
#
# Run by deploy.yml's tag-release job, which pushes the commit to main and
# tags it. Tested by tests/release_scripts.sh.
set -euo pipefail

level=${1:-}
note=${2:-}
dir=$(cd "$(dirname "$0")" && pwd)

current=$(sed -n 's/^#define GSK_VERSION_STRING "\(.*\)"$/\1/p' include/gssk.h)
next=$("$dir/next_version.sh" "$current" "$level")

python3 - "$note" <<'PY'
import re, sys
note = sys.argv[1]
path = 'docs/CHANGELOG.md'
text = open(path, encoding='utf-8').read()
start = text.index('## [Unreleased]')
end = text.index('\n---\n', start)
lines = text[start:end].split('\n')
heading, body = lines[0], '\n'.join(lines[1:])
# The block's ### sections; drop any that has a heading and nothing under it.
kept = []
for section in re.split(r'(?m)^(?=### )', body):
    section = section.strip('\n')
    if not section:
        continue
    title, _, content = section.partition('\n')
    if title.startswith('### ') and not content.strip():
        continue
    kept.append(section)
if not any(re.search(r'(?m)^- ', s) for s in kept):
    if not note:
        sys.exit('auto_release: [Unreleased] has no entries and no note was given')
    kept = [f'### Changed\n\n- {note}']
block = heading + '\n\n' + '\n\n'.join(kept) + '\n'
open(path, 'w', encoding='utf-8').write(text[:start] + block + text[end:])
PY

"$dir/release.sh" "$next" > /dev/null
echo "$next"
