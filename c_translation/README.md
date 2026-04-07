# fmm2d C translation

This is a fork of [flatironinstitute/fmm2d](https://github.com/flatironinstitute/fmm2d).
Everything outside of `c_translation/` is unchanged from upstream.

## Why this exists

The end goal is a **WebAssembly build of fmm2d** so the library can run
in the browser. Compiling the existing Fortran sources to WASM is not
currently practical: an earlier proof-of-concept established that the
LFortran WASM backend (the only Fortran-to-WASM toolchain available)
fails on basic patterns used throughout fmm2d (`goto`, mixed real/complex
subroutine arguments, missing complex math intrinsics in the WASM
runtime). Those are upstream issues that would take months of LFortran
work to fix.

The pragmatic alternative is to **port the relevant subset of fmm2d to
C, then compile the C to WASM** with `clang --target=wasm32-wasi` (or
emscripten). C-to-WASM is a mature, well-tested path; the only piece
missing was the C source code itself. That's what this directory
provides.

## What's in `c_translation/`

A C99 port of the subset of fmm2d needed to support the MATLAB entry
point [`matlab/rfmm2d.m`](../matlab/rfmm2d.m) — the real-valued 2D
Laplace Fast Multipole Method. The C library is a **drop-in replacement**
for the corresponding Fortran objects: same symbol names, same calling
convention, same numerical results within machine precision. The
WASM build itself is not yet wired up, but every translated routine
has been verified bit-for-bit against the Fortran original on x86 so
that the C source is known-correct before it gets compiled to a new
target.

For instructions on continuing the translation (porting more entry
points, adding routines), see [TRANSLATION_GUIDE.md](TRANSLATION_GUIDE.md).

## Status

11 of the original Fortran source files have been translated. Every
translated routine has been verified against its Fortran original with
a bit-for-bit differential test, and the assembled C library has been
verified end-to-end against `test/laplace/test_rfmm2d.f` and a custom
matlab-path smoke test.

| # | Source | C file | Routines | Lines |
|---|---|---|---|---|
| 1 | [src/helmholtz/hndiv2d.f](../src/helmholtz/hndiv2d.f) | [src/hndiv2d.c](src/hndiv2d.c) | 1 | 39 → 51 |
| 2 | [src/common/cumsum.f](../src/common/cumsum.f) | [src/cumsum.c](src/cumsum.c) | 3 | 239 → ~150 |
| 3 | [src/common/fmmcommon2d.f](../src/common/fmmcommon2d.f) | [src/fmmcommon2d.c](src/fmmcommon2d.c) | 3 of 5 | 166 → ~120 |
| 4 | [src/laplace/l2dterms.f](../src/laplace/l2dterms.f) | [src/l2dterms.c](src/l2dterms.c) | 1 of 6 | 564 → ~80 |
| 5 | [src/laplace/cauchykernels2d.f](../src/laplace/cauchykernels2d.f) | [src/cauchykernels2d.c](src/cauchykernels2d.c) | 9 | 683 → ~430 |
| 6 | [src/laplace/laprouts2d.f](../src/laplace/laprouts2d.f) | [src/laprouts2d.c](src/laprouts2d.c) | 16 | 1283 → ~700 |
| 7 | [src/common/tree_routs2d.f](../src/common/tree_routs2d.f) | [src/tree_routs2d.c](src/tree_routs2d.c) | 7 of 8 | 711 → ~450 |
| 8 | [src/common/pts_tree2d.f](../src/common/pts_tree2d.f) | [src/pts_tree2d.c](src/pts_tree2d.c) | 6 | 1250 → ~1115 |
| 9 | [src/laplace/cfmm2d.f](../src/laplace/cfmm2d.f) | [src/cfmm2d.c](src/cfmm2d.c) | 5 | 1689 → ~1400 |
| 10 | [src/laplace/cfmm2d_ndiv.f](../src/laplace/cfmm2d_ndiv.f) | [src/cfmm2d_ndiv.c](src/cfmm2d_ndiv.c) | 1 | 435 → ~340 |
| 11 | [src/laplace/rfmm2d_ndiv.f](../src/laplace/rfmm2d_ndiv.f) | [src/rfmm2d_ndiv.c](src/rfmm2d_ndiv.c) | 1 | 222 → ~180 |

Routines listed as "N of M" mean only the routines reachable from the
`rfmm2d` call graph were translated; the rest are unused by this entry
point and were intentionally skipped.

## Layout

```
c_translation/
├── README.md                     ← you are here
├── TRANSLATION_GUIDE.md          ← detailed guide for continuing the work
├── Makefile                      ← standalone build (not coupled to parent)
├── include/                      ← C headers
│   ├── fmm2d_c.h                 ← shared definitions (FNAME, fint, fcomplex, FA2/FA3)
│   └── <name>.h                  ← one per translated file
├── src/                          ← C translations (one .c per .f)
├── test/                         ← Fortran differential-test drivers
│   └── test_<name>.f             ← compares the C translation to the Fortran reference
├── e2e/                          ← end-to-end smoke test
│   └── test_rfmm2d_e2e.f         ← exercises the matlab call path through the C drop-in
└── build/                        ← (gitignore-able) intermediate objects and binaries
```

## Building and testing

The Makefile assumes the parent fmm2d static library has already been
built (it links against `../lib-static/libfmm2d.a` for cross-file
dependencies and reference symbols).

```bash
# 1. Build the parent Fortran library once (from the repo root):
cd ..
make lib OMP=OFF
gfortran -fPIC -O3 -march=native -funroll-loops -std=legacy -w \
    -c src/common/hkrand.o src/common/dlaran.o   # test helpers
cd c_translation/

# 2. Build the C translations and the per-file diff tests:
make

# 3. Run all 11 differential tests:
make tests
# Expected output: 11 lines of "PASS: ..."

# 4. Run the end-to-end tests through the C drop-in library:
make e2e
# Expected output:
#   ==> build/end2end_rfmm2d (existing test_rfmm2d.f via C drop-in)
#       PASS (output in build/end2end_rfmm2d.log)
#   ==> build/end2end_rfmm2d_ndiv (matlab path smoke test)
#       charge potential at targets: ... rel err= 5.7E-10
#       PASS: rfmm2d_ndiv end-to-end smoke test

# 5. Build the drop-in objects (-O3 with bare Fortran symbol names):
make dropin
# Produces build/dropin/*.o that you can link in place of the
# corresponding objects from libfmm2d.a.

# 6. Clean:
make clean
```

### Differential tests

Each translated file has a Fortran 77 test driver that calls **both**
the original Fortran routine and the C translation on identical inputs
and asserts bit-for-bit equality of every output element.

Both sides of the diff test are compiled at `-O0`. This is essential:
gcc and gfortran at `-O3` don't agree on the last bit of complex
arithmetic (vectorization, fused multiply-add, and operation reordering
all differ subtly), so a sane bit-equal comparison requires both
compilers in their most predictable mode. The drop-in build (`make
dropin`) re-compiles at `-O3` for production speed; the resulting
~1 ULP differences from the Fortran library are well below the FMM
precision target.

### End-to-end tests

`make e2e` builds two integration tests against the C drop-in library:

1. **`build/end2end_rfmm2d`** — runs the existing
   [`test/laplace/test_rfmm2d.f`](../test/laplace/test_rfmm2d.f)
   (which compares FMM output against direct summation) but linked
   against the C drop-in objects instead of the Fortran ones. The
   call graph exercised is:

       test_rfmm2d.f
           ↓ Fortran (rfmm2dwrap.f, rfmm2d.f — not translated)
       cfmm2d  (C)
       cfmm2dmain  (C)
       l2dformmpc/d/cd, l2dmpmp, l2dmploc, l2dlocloc, ...  (all C)
       c2d_directcp/cg/ch/...  (C)
       pts_tree_mem, pts_tree_build, ...  (all C)

2. **`build/end2end_rfmm2d_ndiv`** — a custom driver in
   [`e2e/test_rfmm2d_e2e.f`](e2e/test_rfmm2d_e2e.f) that calls
   `rfmm2d_ndiv` directly (matching what
   [`matlab/rfmm2d.m`](../matlab/rfmm2d.m) does via the MEX interface).
   It compares the FMM output against a brute-force direct sum.

Both tests produce errors at ~1e-10 relative, which is below the
requested `eps = 1e-6` precision and consistent with the underlying
algorithm.

## Calling convention

Translated functions use the gfortran ABI: lowercase name with a
trailing underscore, all arguments passed by pointer, multidimensional
arrays in column-major layout. The header [`include/fmm2d_c.h`](include/fmm2d_c.h)
defines the helpers:

```c
typedef int32_t fint;            // Fortran default integer
typedef double _Complex fcomplex; // Fortran complex *16

#define FA2(i, j, ld1)        // 1-indexed col-major (i, j) with leading dim ld1
#define FA3(i, j, k, ld1, ld2) // 1-indexed col-major (i, j, k)

#ifdef FMM2D_DROP_IN
#define FNAME(x) x##_         // exports `cumsum_`, `cfmm2d_`, etc.
#else
#define FNAME(x) x##_c_       // exports `cumsum_c_`, `cfmm2d_c_`, etc. (diff-test mode)
```

The `FNAME` macro is the trick that lets the same `.c` file serve both
roles: in diff-test mode every translated function gets a `_c_` suffix
so it can coexist with its Fortran original in one binary; in drop-in
mode (`-DFMM2D_DROP_IN`) it gets the canonical Fortran name and replaces
the original at link time.

## What's not done

- **OpenMP parallelism.** All `c$omp` directives in the original are
  comments. The C library is sequential. Adding OpenMP (`#pragma omp`)
  is a separate pass.
- **Other entry points.** Only the routines reachable from
  `rfmm2d_ndiv` are translated. The Helmholtz, complex-Laplace,
  Stokes, biharmonic, and modified-biharmonic FMMs are not.
- **Static archive packaging.** `make dropin` produces individual `.o`
  files but does not bundle them into a `libfmm2d_c.a`. Trivial follow-up.
- **MATLAB MEX glue.** The C library exports the right symbols, but
  rebuilding [`matlab/fmm2d.c`](../matlab/fmm2d.c) against it has not
  been attempted.

For the gory details of how the translation was done — and how to
continue it for the rest of the library — see
[TRANSLATION_GUIDE.md](TRANSLATION_GUIDE.md).
