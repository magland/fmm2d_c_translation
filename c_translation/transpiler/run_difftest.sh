#!/usr/bin/env bash
# run_difftest.sh <name> <fortran-source>
#
# Generates <name>.c / <name>.h with transpile.py, builds them at -O0 together
# with the existing Fortran diff-test driver (test/test_<name>.f) and the
# -O0 Fortran reference object, then runs the test. Bit-for-bit PASS means
# the transpiler output is numerically identical to the Fortran original.
set -euo pipefail

NAME="$1"
FSRC="$2"
ONLY="${3:-}"                          # optional comma-separated routine subset
ONLYARG=()
[ -n "$ONLY" ] && ONLYARG=(--only "$ONLY")

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CT="$(cd "$HERE/.." && pwd)"          # c_translation/
ROOT="$(cd "$CT/.." && pwd)"          # repo root
OUT="$HERE/out"
mkdir -p "$OUT"

CFLAGS="-std=c99 -O0 -g -fPIC -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -Wno-strict-prototypes"
FFLAGS_REF="-O0 -g -std=legacy -fPIC -w"

FMMLIB="$ROOT/lib-static/libfmm2d.a"
TESTOBJS="$ROOT/src/common/hkrand.o $ROOT/src/common/dlaran.o"
# EXTRAOBJS: extra Fortran reference objects to link (for cross-file deps not
# present in libfmm2d.a, e.g. the modified-biharmonic sources). Optional.
EXTRAOBJS="${EXTRAOBJS:-}"

echo "==> generating $NAME.c / $NAME.h"
python3 "$HERE/transpile.py" "$FSRC" --basename "$NAME" "${ONLYARG[@]}" > "$OUT/$NAME.c"
python3 "$HERE/transpile.py" "$FSRC" --basename "$NAME" --header "${ONLYARG[@]}" > "$OUT/$NAME.h"

echo "==> compiling generated C (-O0)"
gcc $CFLAGS -I"$CT/include" -I"$OUT" -c "$OUT/$NAME.c" -o "$OUT/$NAME.o"

echo "==> compiling Fortran reference (-O0)"
gfortran $FFLAGS_REF -c "$FSRC" -o "$OUT/${NAME}_ref.o"

echo "==> linking diff test"
gfortran $FFLAGS_REF "$CT/test/test_$NAME.f" "$OUT/$NAME.o" "$OUT/${NAME}_ref.o" \
    $TESTOBJS $EXTRAOBJS "$FMMLIB" -lm -o "$OUT/test_$NAME"

echo "==> running"
"$OUT/test_$NAME"
