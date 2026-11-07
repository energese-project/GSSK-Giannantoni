#!/bin/sh
# scripts/release_level.sh — how far a merge to main bumps the version.
#
# Usage: scripts/release_level.sh [label ...]
#   with the merged PR's labels as arguments, and the merge commit's message
#   on stdin (for a push with no PR).
# Prints patch, minor, major or skip.
#
#   release:major / release:minor / release:patch   on the PR: that bump
#   release:skip                                    on the PR: no release
#   "Release: major|minor|patch|skip"               a line in the commit message,
#                                                   when the PR has no release label
#   nothing                                         patch
#
# Two bump labels take the larger. release:skip with a bump label is a
# contradiction and fails, rather than guessing which one was meant.
# Used by deploy.yml's tag-release job; tested by tests/release_scripts.sh.
set -eu

# Read the whole commit message first, even when a label decides: exiting
# with stdin unread kills the writer with SIGPIPE, and under the workflow's
# pipefail (`git log ... | release_level.sh`) that fails the release step.
message=$(cat)

rank() {
  case $1 in patch) echo 1 ;; minor) echo 2 ;; major) echo 3 ;; *) echo 0 ;; esac
}

best=''
skip=0
for label in "$@"; do
  case $label in
    release:skip) skip=1 ;;
    release:patch | release:minor | release:major)
      level=${label#release:}
      if [ "$(rank "$level")" -gt "$(rank "${best:-none}")" ]; then best=$level; fi ;;
  esac
done

if [ "$skip" -eq 1 ] && [ -n "$best" ]; then
  echo "release_level: the PR has both release:skip and release:$best" >&2
  exit 2
fi
if [ "$skip" -eq 1 ]; then echo skip; exit 0; fi
if [ -n "$best" ]; then echo "$best"; exit 0; fi

# No label: a "Release: <level>" line in the commit message, else patch.
trailer=$(printf '%s\n' "$message" | (grep -E '^Release: *(major|minor|patch|skip) *$' || true) | tail -n 1 | sed -E 's/^Release: *([a-z]+) *$/\1/')
echo "${trailer:-patch}"
