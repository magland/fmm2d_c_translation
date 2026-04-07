#!/bin/bash
# Build fmm2d as a native shared library for use with numbl.
# Produces fmm2d.so (Linux) or fmm2d.dylib (macOS) in this directory.
#
# Currently only the rfmm2d entry point is exposed (see rfmm2d_wrapper.c).
#
# Unlike the WASM build, the native build links against the existing
# Fortran library (top-level `make lib OMP=OFF`) rather than the C
# translation. The wrapper calls hndiv2d_ and rfmm2d_ndiv_ — those are
# the bare Fortran symbols exported by libfmm2d.a, so the wrapper code
# is identical to the WASM build.
#
# Usage:
#   cd fmm2d/matlab/numbl && bash build_native.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FMM2D_SRC="${FMM2D_SRC:-$(cd "$SCRIPT_DIR/../.." && pwd)}"

# Detect platform.
case "$(uname -s)" in
  Linux*)  EXT="so";    SHARED_FLAGS="-shared -fPIC" ;;
  Darwin*) EXT="dylib"; SHARED_FLAGS="-dynamiclib" ;;
  *)       echo "Unsupported OS: $(uname -s)" >&2; exit 1 ;;
esac

echo "fmm2d source: $FMM2D_SRC"
echo "Output:       fmm2d.$EXT"

# Step 1: build the Fortran static library if it isn't already there.
LIBFMM2D="$FMM2D_SRC/lib-static/libfmm2d.a"
if [ ! -f "$LIBFMM2D" ]; then
  echo "=== Building libfmm2d.a (Fortran) ==="
  pushd "$FMM2D_SRC" > /dev/null

  # Generate a make.inc with -fPIC and OpenMP off (we serialize for numbl
  # since the JS runtime is single-threaded anyway). Same recipe as the
  # mip-core compile.m but without the matlab MEX bits.
  if [ "$EXT" = "dylib" ]; then
    LIBGFORTRAN_NAME="libgfortran.dylib"
  else
    LIBGFORTRAN_NAME="libgfortran.so"
  fi

  cat > make.inc <<EOF
FDIR=\$(shell dirname \`gfortran --print-file-name $LIBGFORTRAN_NAME\`)
MFLAGS+=-L\${FDIR}
OMPFLAGS=
OMPLIBS=
FFLAGS=-fPIC -O3 -funroll-loops -std=legacy -w
EOF

  make lib OMP=OFF
  popd > /dev/null
fi

# Step 2: locate libgfortran for runtime linking.
LIBGFORTRAN_DIR=$(dirname "$(gfortran --print-file-name=lib${EXT/dylib/dylib}gfortran.${EXT} 2>/dev/null || gfortran --print-file-name=libgfortran.${EXT})")

# Step 3: compile the wrapper.
BUILD_DIR="$SCRIPT_DIR/build_native"
mkdir -p "$BUILD_DIR"
WRAPPER_OBJ="$BUILD_DIR/rfmm2d_wrapper.o"
echo "=== Compiling rfmm2d_wrapper.c ==="
gcc -O3 -std=c99 -fPIC -fvisibility=hidden \
    -Wall -Wextra -Wno-unused-parameter \
    -c "$SCRIPT_DIR/rfmm2d_wrapper.c" -o "$WRAPPER_OBJ"

# Step 4: link the shared library. The wrapper only references
# rfmm2d_ndiv_ and hndiv2d_, but those routines pull in (transitively)
# the rest of the rfmm2d call graph via cfmm2d_ndiv_, etc. — and ld
# resolves all of those out of libfmm2d.a in the normal way, so
# --whole-archive isn't needed. We DO want to hide the Fortran symbols
# from the dynamic export table though, so the only public ABI of
# fmm2d.so is rfmm2d_w. On Linux that's a version script; on macOS
# it's an exported_symbols_list.
echo "=== Linking fmm2d.$EXT ==="
if [ "$EXT" = "so" ]; then
  VSCRIPT="$BUILD_DIR/fmm2d.ver"
  cat > "$VSCRIPT" <<'VEOF'
{ global: rfmm2d_w; local: *; };
VEOF
  EXPORT_FLAGS="-Wl,--version-script=$VSCRIPT"
else
  EXPORT_LIST="$BUILD_DIR/fmm2d.exports"
  echo "_rfmm2d_w" > "$EXPORT_LIST"
  EXPORT_FLAGS="-Wl,-exported_symbols_list,$EXPORT_LIST"
fi

g++ $SHARED_FLAGS \
    -fvisibility=hidden \
    "$WRAPPER_OBJ" \
    "$LIBFMM2D" \
    $EXPORT_FLAGS \
    -L"$LIBGFORTRAN_DIR" -lgfortran -lquadmath -lm \
    -o "$SCRIPT_DIR/fmm2d.$EXT"

echo "=== Built fmm2d.$EXT ($(wc -c < "$SCRIPT_DIR/fmm2d.$EXT") bytes) ==="
