# A deterministic Fortran→C transpiler for fmm2d

This directory holds an **auto-generating** alternative to the hand-written
C in [`../src/`](../src). Instead of translating each Fortran file by hand
(or by one-shot LLM), the transpiler parses the Fortran and emits C
mechanically, following the exact conventions documented in
[`../TRANSLATION_GUIDE.md`](../TRANSLATION_GUIDE.md).

> **The transpiler implementation now lives in the standalone `fortc`
> package** (a sibling checkout at `../../../fortc`, `pip install -e`-able).
> The `transpile.py` here is a thin shim that delegates to it, so the build
> scripts below keep working unchanged. See `fortc/README.md` for the package
> docs and its general (non-fmm2d) usage.

The generated C is validated against the **same** per-file differential
tests as the hand translation: build the generated `.c` and the Fortran
original both at `-O0`, run identical inputs through both, assert
bit-for-bit equality.

## Why a deterministic transpiler

The hardest correctness constraint — bit-for-bit `-O0` agreement — is a
*structural* property: it requires preserving Fortran's exact operation
order (left-to-right associativity, `1/(z/r)` vs `r/z`, `(-1)**j`). A
parser-based emitter that never reassociates gets this **for free**, where
a human or LLM gets it only by vigilance. That, plus reproducibility (re-run
when upstream changes) and the existing diff tests as a ready-made oracle,
is why the deterministic route was chosen.

## Architecture

```
Fortran .f ──▶ fparser2 ──▶ Scope/Sym  ──▶ Emitter ──▶ C .c / .h
              (front end)   (symbol      (back end)
                            table)
```

- **Front end — fparser2.** A faithful, source-level Fortran AST (it keeps
  `2 * INT(base / 2D0 + .9999D0)` with correct associativity, labels, and
  GOTOs). No text-scraping of compiler debug output.
- **Symbol table — `Scope` / `Sym`.** Walks the `Specification_Part` to
  resolve each name's C type, intent, array shape, `parameter`/`data`
  values, and `allocatable`-ness. This is the only semantic layer; the
  front end and back end are otherwise decoupled, so the parser substrate
  could be swapped later without touching the emitter.
- **Back end — `Emitter`.** Recursive node→C dispatch encoding the
  TRANSLATION_GUIDE rules.

Every unsupported construct raises `Unsupported(node)` — coverage gaps are
**loud, never silent**.

## Translation rules currently encoded

| Concern | Rule |
|---|---|
| Symbol names | `void FNAME(name)(...)` so one source serves diff-test (`name_c_`) and drop-in (`name_`) roles |
| Args | by pointer; `intent(in)` → `const`; scalar dummy read/write → `(*x)` |
| Same-file calls | `FNAME(callee)(...)` |
| Cross-file calls | bare `callee_(...)` |
| Actual args | array/allocatable → base pointer; scalar dummy → pointer; scalar local/param → `&x`; literal/expr → materialized temp |
| Types | `integer`→`fint`, `real*8`/`double precision`→`double`, `complex*16`→`fcomplex` |
| Implicit typing | honors `implicit none`, `implicit real*8(a-h,o-z)`, and the default I-N rule; body-referenced locals discovered and typed |
| 1-D arrays | `a(i)` → `a[i-1]` (lower-bound aware) |
| 2-/3-D arrays | `FA2`/`FA3` macros (1-based dims); general column-major offset for 0-based / higher-rank |
| Complex | `complex*16`→`fcomplex`; `dcmplx(a,b)`→`(a)+(b)*I`; real-vs-complex dispatch for `abs`/`log`/`exp`/`sqrt`/`dble` |
| Power `**` | complex base→`cpow`, real→`pow`, `(-1)**n`→parity |
| Routine subset | `--only a,b` emits a subset (reachable routines) while keeping all names for same-file call resolution |
| Loops | `do i=1,n` → `for(i=1;i<=n;i++)`; `do while` → `while`; negative literal step handled |
| Control flow | `if`/`else if`/`else`, single-line `IF`, `goto`/labels (`Lnnn:`), `return`, `exit`→`break`, `cycle`→`continue` |
| Expressions | binary ops + relational/logical; operation order preserved verbatim |
| Intrinsics | `int`/`nint`,`dble`,`abs`,`sqrt`,`log`,`exp`,`atan`/`atan2`/`sin`/`cos`/`tan`/`sign`/`aint`/`anint` (+`d*`/`cd*`); `max`/`min`/`mod`; `maxval`/`minval`/`sum` over sections → generated reduction loop |
| Complex literals | `(0.0d0,1.0d0)` → `(0.0 + (1.0)*I)`; signed literals handled |
| Memory | `allocate`(multi-target) → `malloc`; `deallocate` → `free`; allocated bounds captured for later indexing |
| Automatic arrays | local arrays sized by a dummy/expr → C99 VLA |
| Cross-file protos | `extern void callee_();` emitted for every cross-file `CALL` (gcc-14 needs the declaration) |
| Declarations | locals, `parameter` (const-evaluated for array sizes), `data` inits |
| Stripped | OpenMP comment-directives (drop out as Fortran comments); `prini`/`prinf`/`prin2` logging removed; `cpu_time`/`second()`/`omp_get_wtime()` → deterministic `0.0` (diff tests don't compare `timeinfo`) |

## Usage

```bash
python3 transpile.py ../../src/common/next235.f            # -> next235.c on stdout
python3 transpile.py ../../src/common/next235.f --header   # -> next235.h on stdout

# generate + build + run the existing diff test, bit-for-bit:
./run_difftest.sh next235 ../src/common/next235.f
```

`run_difftest.sh` requires the parent library and test helpers (see
[`../README.md`](../README.md)):

```bash
cd ../.. && make lib OMP=OFF
gfortran -fPIC -O3 -march=native -funroll-loops -std=legacy -w -c src/common/hkrand.f -o src/common/hkrand.o
gfortran -fPIC -O3 -march=native -funroll-loops -std=legacy -w -c src/common/dlaran.f -o src/common/dlaran.o
```

## Status — complete

**All 36 translated Fortran files generate and pass their diff tests
bit-for-bit at `-O0`** — i.e. the transpiler reproduces the entire hand
translation deterministically. From `c_translation/` run:

```bash
make transpiler-check      # all 36 per-file diff tests, bit-for-bit
make transpiler-e2e        # 4 matlab-path end-to-end drop-in tests
make generate              # regenerate every .c/.h into transpiler/generated/
```

| Entry point | Files | Result |
|---|---|---|
| Laplace (rfmm2d/cfmm2d/lfmm2d) | next235, cumsum, hndiv2d, fmmcommon2d, l2dterms, cauchykernels2d, lapkernels2d, rlapkernels2d, laprouts2d, tree_routs2d, pts_tree2d, cfmm2d, cfmm2d_ndiv, rfmm2d_ndiv, lfmm2d_ndiv | ✅ all bit-for-bit |
| Helmholtz (hfmm2d) | cdjseval2d, hank103, dfft_threadsafe, h2dcommon, h2dterms, helmkernels2d, helmrouts2d, wideband2d, hfmm2d, hfmm2d_ndiv | ✅ all bit-for-bit |
| Biharmonic | bhndiv2d, bh2dterms, bhkernels2d, bhrouts2d, bhfmm2d | ✅ all bit-for-bit |
| Stokes (stfmm2d) | stfmm2d, stokkernels2d | ✅ all bit-for-bit |
| Modified biharmonic (mbhfmm2d) | mbhgreen2d, mbhkernels2d, mbhrouts2d, mbhfmm2d | ✅ all bit-for-bit |

**End-to-end (`run_e2e.sh`):** the whole library is generated, compiled at
`-O3 -DFMM2D_DROP_IN` (bare Fortran symbol names) and linked into the
matlab-path smoke tests in [`../e2e`](../e2e) as a drop-in replacement for the
Fortran. All pass at the algorithm's accuracy, matching the hand translation:
`rfmm2d` 5.7e-10, `cfmm2d` 2.3e-10, `lfmm2d` 1.8e-9, `stfmm2d` 1.4e-10 (eps 1e-6).

### What it took, milestone by milestone

- **M1** the pipeline (next235, cumsum, hndiv2d): scalar/array dummies,
  control flow, `goto`/labels, `parameter`/`data`, `allocate`.
- **M2** numerics: implicit typing, complex arithmetic (`complex*16`/`double
  complex`, `dcmplx`, `dimag`/`dble`/`dconjg`/`cdabs`, real-vs-complex
  dispatch), `**`, 0-based / general N-D indexing.
- **M3** tree code + operators: allocatables (bounds from `ALLOCATE`,
  multi-target), C99-VLA automatic arrays, `maxval`/`minval`/`sum` reductions,
  `max`/`min`/`mod`, `exit`/`cycle`, cross-file `extern` prototypes.
- **M4** FMM mains: trig/`nint`/`sign`, complex & signed literals,
  logging/timing stripping (`prini`-family, `cpu_time`/`second`→0).
- **M5** all remaining entry points: `dimension` statements, assumed-size
  `a(0:*)`, arithmetic `IF`, `print` stripping, `allocate(...,stat=)`,
  whole-array and subscripted (designated-init) `DATA`, **case-insensitive
  identifiers** (gfortran lowercases all symbols), **`double complex`** typing,
  **same-file call argument casts** (Fortran's typeless arg passing, e.g. a
  real array landing on a `complex*16` parameter), **identifier mangling** for
  C-keyword/libm collisions (a Fortran variable literally named `pow`), and an
  exact integer-power helper for runtime exponents.

The structural approach paid off repeatedly — every bug the diff-test oracle
caught was a concrete, fixable rule (a `double` temp where an `int` was needed;
`x**2` via `pow()` instead of gfortran's exact `x*x`; `double complex` falling
through the type map and halving array sizes), not a fuzzy "close but not
bit-equal" tail. The guide's error-prone "split `x=x+a+b` into two `+=`" rule
is satisfied automatically (fparser keeps the left-associative tree).

## Known limitations / future work

1. **`free` on every return path** for `allocate`d arrays (today heap
   allocations may leak; harmless for the diff tests, worth fixing for prod).
2. **Precise negative-exponent `**`** beyond what the helper covers (no diff
   test exercises it yet).
3. **`save`/`character`/`equivalence`/`common`** are not needed by any
   currently-translated file; add them node-by-node, gated by a diff test, if
   a future file uses them.
