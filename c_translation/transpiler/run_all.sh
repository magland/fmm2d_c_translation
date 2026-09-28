#!/usr/bin/env bash
# run_all.sh - regression suite for the transpiler.
#
# Runs run_difftest.sh on every file the transpiler is known to handle,
# asserting bit-for-bit -O0 agreement with the Fortran original. Each row is
#   <name> <fortran-source> [comma,separated,routine,subset]
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/files.sh"

# map basename -> source path (for building cross-file reference objects)
declare -A SRCOF
for row in "${FILES[@]}"; do IFS='|' read -r n s _ _ <<< "$row"; SRCOF[$n]="$s"; done
CT="$(cd "$HERE/.." && pwd)"            # c_translation/ (sources are ../src/...)
REFDIR="$HERE/out/xrefs"; mkdir -p "$REFDIR"
build_ref() {  # basename -> echoes the .o path, building it once
  local dep="$1"
  [ -f "$REFDIR/$dep.o" ] || gfortran -O0 -g -std=legacy -fPIC -w \
      -c "$CT/${SRCOF[$dep]}" -o "$REFDIR/$dep.o" >/dev/null 2>&1
  echo "$REFDIR/$dep.o"
}

pass=0; fail=0; failed=()
for row in "${FILES[@]}"; do
  IFS='|' read -r name src only extra <<< "$row"
  xo=""
  for dep in $extra; do xo="$xo $(build_ref "$dep")"; done
  if out=$(EXTRAOBJS="$xo" bash "$HERE/run_difftest.sh" "$name" "$src" "$only" 2>&1); then
    if echo "$out" | grep -q "PASS"; then
      printf "  PASS  %s\n" "$name"; pass=$((pass+1))
    else
      printf "  FAIL  %s (no PASS line)\n" "$name"; fail=$((fail+1)); failed+=("$name")
    fi
  else
    printf "  FAIL  %s\n" "$name"; fail=$((fail+1)); failed+=("$name")
    echo "$out" | tail -3 | sed 's/^/        /'
  fi
done
echo "----------------------------------------"
echo "  $pass passed, $fail failed"
[ $fail -eq 0 ] || { echo "  failed: ${failed[*]}"; exit 1; }
