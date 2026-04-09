#!/bin/bash
# Build a numbl-compatible WebAssembly module that exposes the fmm2d
# mexFunction (from upstream matlab/fmm2d.c) via a small mex shim.
# Produces fmm2d.wasm in this directory.
#
# The Fortran sources are not built — Fortran-to-WASM is not practical
# with current toolchains, which is why c_translation/ exists. We compile
# each .c file in c_translation/src/ with emcc using -DFMM2D_DROP_IN so
# the symbols are exported under the bare Fortran ABI names
# (rfmm2d_ndiv_, hndiv2d_, cfmm2d_, ...) and link them with the upstream
# matlab/fmm2d.c (also compiled against our mex shim).
#
# Prerequisites: emcc on PATH (Emscripten SDK)
#
# Usage:
#   cd fmm2d/matlab/numbl && bash build_wasm.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FMM2D_SRC="${FMM2D_SRC:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
CT_DIR="$FMM2D_SRC/c_translation"
BUILD_DIR="$SCRIPT_DIR/build_wasm"

if ! command -v emcc &> /dev/null; then
  echo "Error: emcc (Emscripten) not found on PATH." >&2
  echo "Install: https://emscripten.org/docs/getting_started/downloads.html" >&2
  exit 1
fi

if [ ! -d "$CT_DIR/src" ] || [ ! -d "$CT_DIR/include" ]; then
  echo "Error: c_translation/ not found at $CT_DIR" >&2
  echo "This build requires the magland/fmm2d_c_translation fork." >&2
  exit 1
fi

echo "fmm2d source:    $FMM2D_SRC"
echo "c_translation:   $CT_DIR"
echo "Build directory: $BUILD_DIR"

mkdir -p "$BUILD_DIR"

SHIM_INC="-I$SCRIPT_DIR/mex_shim"
# c_translation provides Fortran routines as `void` functions, but
# fmm2d.c declares them with MWF77_RETURN which defaults to `int`.
# Override to `void` so wasm-ld accepts the matching signatures.
DEFS="-DMX_HAS_INTERLEAVED_COMPLEX=1 -DMWF77_UNDERSCORE1 -DMWF77_RETURN=void"

# Step 1: compile each c_translation .c file with -DFMM2D_DROP_IN so it
# exports the bare Fortran ABI symbol name (e.g. rfmm2d_ndiv_, not
# rfmm2d_ndiv_c_).
CFLAGS="-O3 -std=c99 -DFMM2D_DROP_IN -msimd128 -Wno-unused-parameter -Wno-unused-variable"
INCLUDES="-I$CT_DIR/include"

OBJS=()
for src in "$CT_DIR"/src/*.c; do
  base=$(basename "$src" .c)
  obj="$BUILD_DIR/${base}.o"
  echo "  CC  $base.c"
  emcc $CFLAGS $INCLUDES -c "$src" -o "$obj"
  OBJS+=("$obj")
done

# Stubs for the Fortran entry points c_translation does not yet port
# (hfmm2d_ndiv, h2d_*, l2d_*, r2d_*, c2d_*, st2d_* direct evaluators).
# These let fmm2d.o link cleanly; calling them at runtime sets *ier=1.
echo "  CC  wasm_stubs.c"
WASM_STUBS_OBJ="$BUILD_DIR/wasm_stubs.o"
emcc -O3 -std=c99 -msimd128 -Wno-unused-parameter \
     -c "$SCRIPT_DIR/wasm_stubs.c" -o "$WASM_STUBS_OBJ"
OBJS+=("$WASM_STUBS_OBJ")

# Step 2: compile the upstream MEX source against our mex shim.
echo "  CC  fmm2d.c"
FMM2D_OBJ="$BUILD_DIR/fmm2d.o"
emcc -O3 -std=c99 -msimd128 \
     -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable \
     $SHIM_INC $DEFS \
     -c "$FMM2D_SRC/matlab/fmm2d.c" -o "$FMM2D_OBJ"
OBJS+=("$FMM2D_OBJ")

# Step 3: compile the shim implementation (C++).  SUPPORT_LONGJMP=wasm
# is also a compile-time flag — it must be passed here so the setjmp
# call inside mex_dispatch is lowered to wasm-native sjlj instead of
# the emscripten JS shim (which standalone wasm doesn't have).
echo "  CXX mex_shim.cpp"
SHIM_OBJ="$BUILD_DIR/mex_shim.o"
em++ -O3 -msimd128 \
     -Wno-unused-parameter \
     -s SUPPORT_LONGJMP=wasm \
     $SHIM_INC $DEFS \
     -c "$SCRIPT_DIR/mex_shim.cpp" -o "$SHIM_OBJ"
OBJS+=("$SHIM_OBJ")

# Step 4: link into a standalone WASM module. STANDALONE_WASM produces a
# module that can be instantiated by any wasm runtime (numbl runs it under
# its own loader, not emscripten's JS shim).
echo "=== Linking fmm2d.wasm ==="
em++ "${OBJS[@]}" \
  -O3 \
  -msimd128 \
  -s STANDALONE_WASM \
  -s SUPPORT_LONGJMP=wasm \
  --no-entry \
  -s TOTAL_MEMORY=67108864 \
  -s ALLOW_MEMORY_GROWTH=1 \
  -o "$SCRIPT_DIR/fmm2d.wasm"

echo "=== Built fmm2d.wasm ($(wc -c < "$SCRIPT_DIR/fmm2d.wasm") bytes) ==="
