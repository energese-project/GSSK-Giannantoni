#!/bin/sh
# scripts/next_version.sh — the version after <version> for a bump <level>.
#
# Usage: scripts/next_version.sh <major.minor.patch> <patch|minor|major>
#   6.0.1 patch -> 6.0.2    6.0.1 minor -> 6.1.0    6.0.1 major -> 7.0.0
#
# Used by scripts/auto_release.sh; tested by tests/release_scripts.sh.
set -eu

version=${1:-}
level=${2:-}

case $version in
  *[!0-9.]* | .* | *. | *..* | '')
    echo "next_version: '$version' is not major.minor.patch" >&2; exit 2 ;;
esac
IFS=. read -r major minor patch extra <<EOF
$version
EOF
if [ -z "$major" ] || [ -z "$minor" ] || [ -z "$patch" ] || [ -n "${extra:-}" ]; then
  echo "next_version: '$version' is not major.minor.patch" >&2; exit 2
fi

case $level in
  patch) echo "$major.$minor.$((patch + 1))" ;;
  minor) echo "$major.$((minor + 1)).0" ;;
  major) echo "$((major + 1)).0.0" ;;
  *) echo "next_version: level must be patch, minor or major, not '$level'" >&2; exit 2 ;;
esac
