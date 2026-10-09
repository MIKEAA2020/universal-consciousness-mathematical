#!/usr/bin/env bash
# restore_framework.sh — self-healing, reset-resistant WORKING_FRAMEWORK.md store.
#
# The framework (binding constitution of this project) lives in redundant copies:
#   CANONICAL  : framework/WORKING_FRAMEWORK.md in the GitHub repo (survives everything)
#   HOME ZONE  : /home/z/.config/working_framework.md, /home/z/.working_framework.md
#   PROJECT    : /home/z/my-project/framework/WORKING_FRAMEWORK.md
# This script scans all zones, and if ANY single copy survives, rebuilds all the
# others. If multiple survivors disagree, the NEWEST wins and a warning is printed.
#
# Usage:   bash /home/z/my-project/scripts/restore_framework.sh
#          (also present in the repo as scripts/restore_framework.sh)
# Idempotent — safe to run at the start of every session, next to
#          restore_github_pat.sh (which heals the PAT store the same way).

set -u

PROJECT_ZONE="/home/z/my-project"
HOME_ZONE="/home/z"
FILENAME="WORKING_FRAMEWORK.md"

declare -a CANDIDATES=(
  "$HOME_ZONE/.config/working_framework.md"
  "$HOME_ZONE/.working_framework.md"
  "$PROJECT_ZONE/framework/$FILENAME"
)
# auto-detect repo clones (any repo containing framework/WORKING_FRAMEWORK.md)
for d in "$PROJECT_ZONE"/repos/*/framework; do
  [[ -f "$d/$FILENAME" ]] && CANDIDATES+=("$d/$FILENAME")
done

# collect survivors
declare -a SURVIVORS=()
for f in "${CANDIDATES[@]}"; do
  [[ -s "$f" ]] && SURVIVORS+=("$f")
done

if [[ ${#SURVIVORS[@]} -eq 0 ]]; then
  echo "ERROR: no surviving WORKING_FRAMEWORK.md copy found in any zone." >&2
  echo "Recovery floor: git clone https://github.com/MIKEAA2020/universal-consciousness-mathematical.git" >&2
  echo "then re-run this script from the repo's scripts/ directory." >&2
  exit 1
fi

# pick the NEWEST survivor (amendments are dated; newest = most current ruling)
SOURCE=""
for f in "${SURVIVORS[@]}"; do
  if [[ -z "$SOURCE" ]] || [[ "$f" -nt "$SOURCE" ]]; then
    SOURCE="$f"
  fi
done

# warn on checksum divergence among survivors
for f in "${SURVIVORS[@]}"; do
  if ! cmp -s "$SOURCE" "$f"; then
    echo "WARNING: framework copies diverge: '$f' != '$SOURCE' (newest kept; reconcile manually if intended)." >&2
  fi
done

# rebuild every mirror
mkdir -p "$HOME_ZONE/.config" "$PROJECT_ZONE/framework"
TARGETS=(
  "$HOME_ZONE/.config/working_framework.md"
  "$HOME_ZONE/.working_framework.md"
  "$PROJECT_ZONE/framework/$FILENAME"
)
RESTORED=0
for t in "${TARGETS[@]}"; do
  if ! cmp -s "$SOURCE" "$t"; then
    cp "$SOURCE" "$t"
    echo "restored: $t"
    RESTORED=$((RESTORED + 1))
  fi
done
[[ $RESTORED -eq 0 ]] && echo "OK: all framework copies already in sync (source: $SOURCE)." \
                      || echo "OK: framework restored to all zones from: $SOURCE"

# report checksums for audit
sha256sum "$SOURCE" "${TARGETS[@]}" 2>/dev/null | sed 's|/home/z/my-project|<project>|; s|/home/z|<home>|'
