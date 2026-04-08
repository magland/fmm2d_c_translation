#!/bin/bash
# Build fmm2d as a native shared library for use with numbl.
# Produces fmm2d.so (Linux) or fmm2d.dylib (macOS) in this directory.
#
# Exposes the rfmm2d, cfmm2d, lfmm2d, and stfmm2d entry points (see
# rfmm2d_wrapper.c).
#
# Unlike the WASM build, the native build links against the existing
# Fortran library (top-level `make lib OMP=OFF`) rather than the C
# translation. The wrapper calls hndiv2d_, rfmm2d_ndiv_, cfmm2d_ndiv_,
# lfmm2d_ndiv_, and stfmm2d_ — those are the bare Fortran symbols
# exported by libfmm2d.a, so the wrapper code is identical to the
# WASM build.
#
# Usage:
#   cd fmm2d/matlab/numbl && bash build_native.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FMM2D_SRC="${FMM2D_SRC:-$(cd "$SCRIPT_DIR/../.." && pwd)}"

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

  # Generate a make.inc with -fPIC, OpenMP off, and the discovered FC.
  # We serialize for numbl since the JS runtime is single-threaded
  # anyway. The make.inc also pins FC so the makefile uses the
  # versioned binary we found above (instead of the default `gfortran`,
  # which may not exist).
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
OMPFLAGS=
OMPLIBS=
FFLAGS=-fPIC -O3 -funroll-loops -std=legacy -w
EOF

  make lib OMP=OFF FC="$FC"
  popd > /dev/null
fi

# Step 2: compile the wrapper.
BUILD_DIR="$SCRIPT_DIR/build_native"
mkdir -p "$BUILD_DIR"
WRAPPER_OBJ="$BUILD_DIR/rfmm2d_wrapper.o"
echo "=== Compiling rfmm2d_wrapper.c ==="
gcc -O3 -std=c99 -fPIC -fvisibility=hidden \
    -Wall -Wextra -Wno-unused-parameter \
    -c "$SCRIPT_DIR/rfmm2d_wrapper.c" -o "$WRAPPER_OBJ"

# Step 3: link the shared library. We use $FC (gfortran) as the linker
# driver, not gcc/g++, for two reasons:
#
#   1. gfortran knows where its own runtime libraries live (libgfortran,
#      libquadmath, ...) and passes the right -L paths to ld
#      automatically. On Ubuntu CI runners, libgfortran.so is only in
#      the -dev package — only libgfortran.so.5 is installed by default
#      — so a manual `-lgfortran` from gcc/g++ fails to resolve. Letting
#      gfortran do the link sidesteps that entirely.
#
#   2. The upstream fmm2d makefile already builds libfmm2d.so with
#      `gfortran -shared`, so this is the same path the existing native
#      MEX build uses for its shared library.
#
# The wrapper only references rfmm2d_ndiv_ and hndiv2d_, but those
# routines pull in (transitively) the rest of the rfmm2d call graph via
# cfmm2d_ndiv_, etc. — and ld resolves all of those out of libfmm2d.a
# in the normal way, so --whole-archive isn't needed. We DO want to
# hide the Fortran symbols from the dynamic export table though, so the
# only public ABI of fmm2d.{so,dylib} is rfmm2d_w. On Linux that's a
# version script; on macOS it's an exported_symbols_list.
echo "=== Linking fmm2d.$EXT ==="
if [ "$EXT" = "so" ]; then
  VSCRIPT="$BUILD_DIR/fmm2d.ver"
  cat > "$VSCRIPT" <<'VEOF'
{ global: rfmm2d_w; cfmm2d_w; lfmm2d_w; stfmm2d_w; local: *; };
VEOF
  EXPORT_FLAGS="-Wl,--version-script=$VSCRIPT"
else
  EXPORT_LIST="$BUILD_DIR/fmm2d.exports"
  cat > "$EXPORT_LIST" <<'EEOF'
_rfmm2d_w
_cfmm2d_w
_lfmm2d_w
_stfmm2d_w
EEOF
  EXPORT_FLAGS="-Wl,-exported_symbols_list,$EXPORT_LIST"
fi

"$FC" $SHARED_FLAGS \
    "$WRAPPER_OBJ" \
    "$LIBFMM2D" \
    $EXPORT_FLAGS \
    -lm \
    -o "$SCRIPT_DIR/fmm2d.$EXT"

echo "=== Built fmm2d.$EXT ($(wc -c < "$SCRIPT_DIR/fmm2d.$EXT") bytes) ==="
