#!/usr/bin/env bash
# generate.sh - regenerate the whole C library from the Fortran sources.
#
# Writes <name>.c and <name>.h for every handled file into transpiler/generated/.
# This is the deterministic, reproducible counterpart to the hand-written C in
# c_translation/src: re-run it whenever the upstream Fortran changes.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/files.sh"
OUT="${1:-$HERE/generated}"
mkdir -p "$OUT"

n=0
for row in "${FILES[@]}"; do
  IFS='|' read -r name src only extra <<< "$row"
  onlyarg=(); [ -n "$only" ] && onlyarg=(--only "$only")
  python3 "$HERE/transpile.py" "$src" --basename "$name" "${onlyarg[@]}" > "$OUT/$name.c"
  python3 "$HERE/transpile.py" "$src" --basename "$name" --header "${onlyarg[@]}" > "$OUT/$name.h"
  n=$((n+1))
done
echo "generated $n files (.c + .h) into $OUT"
