#!/bin/bash
# remove_tracked_artifacts.sh
# Run from repository root. Removes tracked build artifacts (cached) and commits the deletions.
# WARNING: This will remove the matching files from the repository in the current branch.

set -euo pipefail

BRANCH=$(git rev-parse --abbrev-ref HEAD)
if [ "$BRANCH" != "try/borrarbins" ]; then
  echo "Warning: you are on branch $BRANCH. It's recommended to run this script on branch try/borrarbins. Continue? (y/N)"
  read ans
  if [ "$ans" != "y" ]; then
    echo "Aborted."
    exit 1
  fi
fi

# Patterns to untrack
PATTERNS=("*.mod" "*.o" "*.obj" ".vs" "*.suo" "*.user" "build" "bin" "*.log")

# Build a list of files/dirs tracked that match patterns
MATCHES=()
for p in "${PATTERNS[@]}"; do
  # find tracked files matching the pattern
  while IFS= read -r -d $'\0' f; do
    MATCHES+=("$f")
  done < <(git ls-files -z -- "${p}")
done

if [ ${#MATCHES[@]} -eq 0 ]; then
  echo "No tracked files matched the patterns. Nothing to do."
  exit 0
fi

echo "The following tracked paths will be removed from the repo index (kept on disk):"
printf '%s
' "${MATCHES[@]}"

echo
read -p "Proceed to remove these from the index and commit? (y/N) " confirm
if [ "$confirm" != "y" ]; then
  echo "Aborted."
  exit 1
fi

# Remove each matched path from index (preserve files on disk)
for f in "${MATCHES[@]}"; do
  git rm --cached -r --ignore-unmatch -- "$f" || true
done

git add .gitignore || true

git commit -m "Remove tracked build artifacts and add .gitignore" || true

echo "Committed removals. Push the branch with: git push origin try/borrarbins"
