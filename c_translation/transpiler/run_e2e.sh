#!/usr/bin/env bash
# run_e2e.sh - end-to-end test of the transpiler output.
#
# Generates the whole Laplace call graph with transpile.py, compiles every
# file at -O3 with -DFMM2D_DROP_IN (so each exports the bare Fortran symbol
# name), and links the matlab-path smoke tests in ../e2e against them. Any
# symbol the C set does not provide is resolved from the Fortran libfmm2d.a,
# so this exercises the generated C as a genuine drop-in replacement.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CT="$(cd "$HERE/.." && pwd)"
ROOT="$(cd "$CT/.." && pwd)"
OUT="$HERE/out_dropin"
mkdir -p "$OUT"

CFLAGS="-std=c99 -O3 -march=native -funroll-loops -fPIC -fno-strict-aliasing \
        -fwrapv -DFMM2D_DROP_IN -Wall -Wextra -Wno-unused-parameter \
        -Wno-unused-variable -Wno-unused-but-set-variable -Wno-strict-prototypes"
FFLAGS="-O3 -march=native -funroll-loops -std=legacy -fPIC -w"
FMMLIB="$ROOT/lib-static/libfmm2d.a"
TESTOBJS="$ROOT/src/common/hkrand.o $ROOT/src/common/dlaran.o"

source "$HERE/files.sh"   # the whole generated library, built as one drop-in

echo "==> generating + compiling drop-in objects (-O3 -DFMM2D_DROP_IN)"
OBJS=()
for row in "${FILES[@]}"; do
  IFS='|' read -r name src only extra <<< "$row"
  case "$name" in mbh*) continue;; esac   # no modified-biharmonic e2e driver
  onlyarg=(); [ -n "$only" ] && onlyarg=(--only "$only")
  python3 "$HERE/transpile.py" "$src" --basename "$name" "${onlyarg[@]}" > "$OUT/$name.c"
  python3 "$HERE/transpile.py" "$src" --basename "$name" --header "${onlyarg[@]}" > "$OUT/$name.h"
  gcc $CFLAGS -I"$CT/include" -I"$OUT" -c "$OUT/$name.c" -o "$OUT/$name.o"
  OBJS+=("$OUT/$name.o")
done

run_one() {  # <label> <driver.f>
  local label="$1" drv="$2"
  echo "==> e2e $label"
  gfortran $FFLAGS "$CT/e2e/$drv" "${OBJS[@]}" $TESTOBJS "$FMMLIB" -lm \
      -o "$OUT/e2e_$label"
  "$OUT/e2e_$label"
}

run_one rfmm2d test_rfmm2d_e2e.f
run_one cfmm2d test_cfmm2d_e2e.f
run_one lfmm2d test_lfmm2d_e2e.f
run_one stfmm2d test_stfmm2d_e2e.f
echo "==> all e2e drop-in tests passed"
