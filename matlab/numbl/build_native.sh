#!/bin/bash
# Build a numbl-compatible native shared library that exposes the fmm2d
# mexFunction (from upstream matlab/fmm2d.c) via a small mex shim.
# Produces fmm2d.so (Linux) or fmm2d.dylib (macOS) in this directory.
# All Fortran fmm2d routines are statically linked from libfmm2d.a.
#
# Usage:
#   cd fmm2d/matlab/numbl && bash build_native.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FMM2D_SRC="${FMM2D_SRC:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
BUILD_DIR="$SCRIPT_DIR/build_native"

# Make sure both Homebrew prefixes are on PATH. CI runners often start
# with neither: /opt/homebrew/bin is the ARM brew location and
# /usr/local/bin is the Intel brew location.
export PATH="/opt/homebrew/bin:/usr/local/bin:$PATH"

# Detect platform.
case "$(uname -s)" in
  Linux*)  EXT="so";    SHARED_FLAGS="-shared -fPIC" ;;
  Darwin*) EXT="dylib"; SHARED_FLAGS="-dynamiclib" ;;
  *)       echo "Unsupported OS: $(uname -s)" >&2; exit 1 ;;
esac

# Locate gfortran. GitHub Actions macOS runners ship gcc/gfortran but
# only as versioned binaries (gfortran-11..15) — there is no plain
# `gfortran` symlink. Try the plain name first, then fall back to
# versioned ones in newest-first order.
FC=""
for cand in gfortran gfortran-15 gfortran-14 gfortran-13 gfortran-12 gfortran-11; do
  if command -v "$cand" >/dev/null 2>&1; then
    FC="$cand"
    break
  fi
done
if [ -z "$FC" ]; then
  echo "Error: no gfortran found on PATH (tried gfortran, gfortran-{11..15})." >&2
  echo "PATH=$PATH" >&2
  exit 1
fi

echo "fmm2d source: $FMM2D_SRC"
echo "Output:       fmm2d.$EXT"
echo "Fortran:      $FC ($(command -v "$FC"))"

# Step 1: build the Fortran static library if it isn't already there.
LIBFMM2D="$FMM2D_SRC/lib-static/libfmm2d.a"
if [ ! -f "$LIBFMM2D" ]; then
  echo "=== Building libfmm2d.a (Fortran) ==="
  pushd "$FMM2D_SRC" > /dev/null

  if [ "$EXT" = "dylib" ]; then
    LIBGFORTRAN_NAME="libgfortran.dylib"
  else
    LIBGFORTRAN_NAME="libgfortran.so"
  fi

  cat > make.inc <<EOF
FC=$FC
CC=gcc
FDIR=\$(shell dirname \`$FC --print-file-name $LIBGFORTRAN_NAME\`)
MFLAGS+=-L\${FDIR}
OMPFLAGS=-fopenmp
OMPLIBS=-lgomp
FFLAGS=-fPIC -O3 -funroll-loops -std=legacy -w
EOF

  make lib OMP=ON FC="$FC"
  popd > /dev/null
fi

# Step 2: compile the upstream mwrap-generated MEX source against our
# tiny mex shim header, plus the shim implementation itself.
mkdir -p "$BUILD_DIR"

SHIM_INC="-I$SCRIPT_DIR/mex_shim"
DEFS="-DMX_HAS_INTERLEAVED_COMPLEX=1 -DMWF77_UNDERSCORE1"

FMM2D_OBJ="$BUILD_DIR/fmm2d.o"
SHIM_OBJ="$BUILD_DIR/mex_shim.o"

echo "=== Compiling fmm2d.c ==="
gcc -O3 -std=c99 -fPIC -fvisibility=hidden \
    -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable \
    $SHIM_INC $DEFS \
    -c "$FMM2D_SRC/matlab/fmm2d.c" -o "$FMM2D_OBJ"

echo "=== Compiling mex_shim.cpp ==="
g++ -O3 -fPIC -fvisibility=hidden \
    -Wall -Wextra -Wno-unused-parameter \
    $SHIM_INC $DEFS \
    -c "$SCRIPT_DIR/mex_shim.cpp" -o "$SHIM_OBJ"

# Step 3: link the shared library. Use gfortran as the linker driver so
# it pulls in libgfortran/libquadmath automatically (Ubuntu CI runners
# only have libgfortran.so.5 in the default install; relying on a manual
# `-lgfortran` from gcc/g++ would fail to resolve).
#
# We hide the Fortran symbols from the dynamic export table so the only
# public ABI of fmm2d.{so,dylib} is mex_*/my_*. On Linux that's a
# version script; on macOS it's an exported_symbols_list.
echo "=== Linking fmm2d.$EXT ==="
if [ "$EXT" = "so" ]; then
  VSCRIPT="$BUILD_DIR/fmm2d.ver"
  cat > "$VSCRIPT" <<'VEOF'
{ global: mex_*; my_*; local: *; };
VEOF
  EXPORT_FLAGS="-Wl,--version-script=$VSCRIPT"
else
  EXPORT_LIST="$BUILD_DIR/fmm2d.exports"
  : > "$EXPORT_LIST"
  for sym in mex_alloc_args mex_dispatch mex_free_args mex_free_array \
             mex_get_arg mex_get_classid mex_get_error mex_get_is_complex \
             mex_get_m mex_get_n mex_make_complex_matrix mex_make_double_scalar \
             mex_make_real_matrix mex_make_string mex_make_struct \
             mex_read_complex mex_read_double_scalar mex_read_real \
             mex_read_string mex_set_arg mex_struct_set_field \
             my_malloc my_free; do
    echo "_${sym}" >> "$EXPORT_LIST"
  done
  EXPORT_FLAGS="-Wl,-exported_symbols_list,$EXPORT_LIST"
fi

"$FC" $SHARED_FLAGS \
    "$FMM2D_OBJ" "$SHIM_OBJ" \
    "$LIBFMM2D" \
    $EXPORT_FLAGS \
    -fopenmp \
    -lm -lstdc++ \
    -o "$SCRIPT_DIR/fmm2d.$EXT"

echo "=== Built fmm2d.$EXT ($(wc -c < "$SCRIPT_DIR/fmm2d.$EXT") bytes) ==="
