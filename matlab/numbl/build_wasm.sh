#!/bin/bash
# Build the C-translation of fmm2d as a WebAssembly module for use with numbl.
# Produces fmm2d.wasm in this directory.
#
# Currently only the rfmm2d entry point is exposed (see rfmm2d_wrapper.c).
#
# The Fortran sources are not built — Fortran-to-WASM is not practical with
# current toolchains, which is why c_translation/ exists. We compile each .c
# file in c_translation/src/ with emcc using -DFMM2D_DROP_IN so the symbols
# are exported under the bare Fortran ABI names (rfmm2d_ndiv_, hndiv2d_,
# cfmm2d_, ...) and link them with the rfmm2d wrapper.
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

# Step 2: compile the rfmm2d wrapper.
WRAPPER_OBJ="$BUILD_DIR/rfmm2d_wrapper.o"
echo "  CC  rfmm2d_wrapper.c"
emcc $CFLAGS $INCLUDES -c "$SCRIPT_DIR/rfmm2d_wrapper.c" -o "$WRAPPER_OBJ"
OBJS+=("$WRAPPER_OBJ")

# Step 3: link into a standalone WASM module. STANDALONE_WASM produces a
# module that can be instantiated by any wasm runtime (numbl runs it under
# its own loader, not emscripten's JS shim).
echo "=== Linking fmm2d.wasm ==="
emcc "${OBJS[@]}" \
  -O3 \
  -msimd128 \
  -s STANDALONE_WASM \
  --no-entry \
  -s TOTAL_MEMORY=67108864 \
  -s ALLOW_MEMORY_GROWTH=1 \
  -o "$SCRIPT_DIR/fmm2d.wasm"

echo "=== Built fmm2d.wasm ($(wc -c < "$SCRIPT_DIR/fmm2d.wasm") bytes) ==="
