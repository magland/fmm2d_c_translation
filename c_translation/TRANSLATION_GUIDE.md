# fmm2d → C translation: continuation guide

This document explains how the existing 18-file translation was built
and how to extend it. The intended reader is an AI coding agent (or a
human) who wants to translate additional Fortran sources from
[fmm2d](https://github.com/flatironinstitute/fmm2d) to C.

The work currently covers the call graphs of four entry points:
- `rfmm2d_ndiv` — real-valued Laplace FMM
- `cfmm2d_ndiv` — complex Cauchy-kernel Laplace FMM
- `lfmm2d_ndiv` — log-kernel Laplace FMM with complex densities
- `stfmm2d` / `bhfmm2d` — 2D Stokes / biharmonic FMM

Adding more entry points (Helmholtz `hfmm2d`, modified biharmonic
`mbhfmm2d`) follows the same recipe.

## Table of contents

1. [Goals and scope philosophy](#1-goals-and-scope-philosophy)
2. [The drop-in replacement trick (`FNAME` macro)](#2-the-drop-in-replacement-trick-fname-macro)
3. [Translation conventions](#3-translation-conventions)
4. [The differential test architecture](#4-the-differential-test-architecture)
5. [Build system mechanics](#5-build-system-mechanics)
6. [Recipe: translating one new file](#6-recipe-translating-one-new-file)
7. [Pitfalls (everything that broke during the first 11 files)](#7-pitfalls-everything-that-broke-during-the-first-11-files)
8. [Subagent prompt template](#8-subagent-prompt-template)
9. [Suggested next files](#9-suggested-next-files)
10. [Quick reference: where things live](#10-quick-reference-where-things-live)

---

## 1. Goals and scope philosophy

**Goal:** A C library that can replace the Fortran object files in
`lib-static/libfmm2d.a`, one file at a time, with bit-for-bit
equivalent results (modulo compiler optimization choices).

**Scope reduction:** Only translate routines actually reachable from
the entry point of interest. Skip unreachable routines (`print_tree`,
`l2dterms_far`, `geterrstr`, `ireorderf`, etc.) entirely — don't even
stub them. This is a defensible scope cut: if a future caller needs
those routines, they can be added later, but porting them now wastes
effort.

**Verification first.** Every translated routine has a Fortran-side
test driver that calls both versions on the same inputs and asserts
agreement. The test is the spec; if the test passes bit-for-bit, the
translation is provably correct (modulo coverage of the test inputs).

**No drift.** Strict 1:1. Don't refactor. Don't add error handling.
Don't add asserts. Don't add OpenMP. Don't restructure loops. Don't
"clean up" gotos. The Fortran source is the spec; the C is just a
re-encoding. Cleanup, parallelization, and refactoring all come *after*
correctness is established.

---

## 2. The drop-in replacement trick (`FNAME` macro)

The single most important piece of the translation is in
[`include/fmm2d_c.h`](include/fmm2d_c.h):

```c
#ifdef FMM2D_DROP_IN
#define FNAME(x) x##_
#else
#define FNAME(x) x##_c_
#endif
```

Every translated function is **defined** via `FNAME(...)`:

```c
void FNAME(cumsum)(const fint *n, const fint *a, fint *b) { ... }
```

In **diff-test mode** (the default), `FNAME(cumsum)` expands to
`cumsum_c_`, which can coexist with the gfortran-compiled original
`cumsum_` in the same executable. The test driver calls both and
compares.

In **drop-in mode** (compile with `-DFMM2D_DROP_IN`), `FNAME(cumsum)`
expands to `cumsum_` — the canonical Fortran symbol. The C object
file becomes a drop-in replacement for `src/common/cumsum.o` and can
be linked into any binary that previously used the Fortran version.

This is the only mechanism needed. There is no second build of the
sources, no preprocessor flag soup, no #ifdef branches inside the
translation logic. One source file, two roles.

### Same-file vs cross-file calls

When translating, the rules for *calls* (not just definitions) are:

- **Same-file calls** (one routine in `foo.c` calling another in
  `foo.c`): use `FNAME(name)`. They flip together with the
  definitions.
- **Cross-file calls** (calling a routine defined in another `.c`):
  use the **bare Fortran name** (e.g. `cumsum_`) with an `extern`
  forward declaration. **Do not** use `FNAME` here.

Why bare names for cross-file calls? Because:

- In diff-test mode, the bare name resolves to the trusted Fortran
  reference from `libfmm2d.a`, isolating the test to ONE file at a
  time. (Otherwise a bug in `cumsum.c` would corrupt the diff test
  for `tree_routs2d.c`, since `tree_routs2d` calls `cumsum`.)
- In drop-in mode, the bare name resolves to the C drop-in version of
  the dependency, since both are now exporting the same symbol.

This is the convention that makes the diff tests sound: each test
exercises exactly one new routine against a fully-trusted
reference for everything else.

---

## 3. Translation conventions

### Calling convention

The Fortran ABI used by gfortran on Linux:

- Symbol names are lowercase with a single trailing underscore.
- All arguments are passed by reference (pointer to scalar, or pointer
  to the first element of an array).
- There is no return value for `subroutine` (return value goes through
  an output pointer parameter).
- Multidimensional arrays are stored column-major.

### Type mapping

| Fortran type | C type | Notes |
|---|---|---|
| `integer` | `fint *` | `fint = int32_t`, the default Fortran integer |
| `real *8`, `double precision` | `double *` | |
| `complex *16` | `fcomplex *` | `fcomplex = double _Complex` |
| `logical` | not used in fmm2d2D path | use `fint` if encountered |

### Indexing macros

The header [`include/fmm2d_c.h`](include/fmm2d_c.h) defines:

```c
#define FA2(i, j, ld1)         (((j) - 1) * (ld1) + ((i) - 1))
#define FA3(i, j, k, ld1, ld2) ((((k) - 1) * (ld2) + ((j) - 1)) * (ld1) + ((i) - 1))
#define FA4(i, j, k, l, ld1, ld2, ld3) ...
```

These take **1-based** Fortran subscripts and return the linear C
offset for column-major storage. The convention is to keep loop
variables 1-based in the C code so that the visual structure mirrors
the Fortran source — this makes off-by-one bugs much easier to spot in
review.

```fortran
do j = 1, nt
   do i = 1, ns
      pot(ii, j) = pot(ii, j) + charge(ii, i) * something
   enddo
enddo
```
becomes
```c
for (j = 1; j <= nt_v; j++) {
    for (i = 1; i <= ns_v; i++) {
        pot[FA2(ii, j, nd_v)] += charge[FA2(ii, i, nd_v)] * something;
    }
}
```

For arrays with **0-based** Fortran subscripts (e.g. `mpole(nd, 0:nterms)`,
`carray(0:ldc, 0:ldc)`, `boxsize(0:nlevels)`), the FA macros don't fit
because they assume 1-based indices. Define a small per-file helper:

```c
/* mpole(ii, n) where ii is 1-based and n is 0-based, leading dim nd */
#define MIDX(ii, n, nd) ((n) * (nd) + ((ii) - 1))

/* carray(l, m) where both are 0-based, leading dim is ldc+1 */
#define CIDX(l, m, ldc) ((m) * ((ldc) + 1) + (l))

/* laddr(k, ilev) where k is 1-based and ilev is 0-based, leading dim 2 */
#define LADDR(k, ilev) ((ilev) * 2 + ((k) - 1))
```

For 1D 0-based arrays like `boxsize(0:nlevels)`, just use direct
indexing: `boxsize[ilev]`.

### Slice arguments

Many Fortran calls pass a slice into the middle of an array:

```fortran
call l2dformmpc(nd, rscales(ilev), sourcesort(1, istart), npts,
                chargesort(1, istart), centers(1, ibox), nterms(ilev),
                rmlexp(iaddr(1, ibox)))
```

In C, these become pointer arithmetic:

```c
l2dformmpc_(&nd, &rscales[ilev], &sourcesort[FA2(1, istart, 2)], &npts,
            &chargesort[FA2(1, istart, nd_v)],
            &centers[FA2(1, ibox, 2)], &nterms[ilev],
            (fcomplex *)&rmlexp[iaddr[FA2(1, ibox, 2)] - 1]);
```

Note the `iaddr[FA2(1, ibox, 2)] - 1`: the value stored in `iaddr` is
a 1-based Fortran index into `rmlexp`, so we subtract 1 to get the C
0-based offset. The cast to `fcomplex *` is because `rmlexp` is
declared as `real *8` in Fortran but holds complex multipole/local
expansion coefficients packed as pairs of doubles — the laprouts2d
routines reinterpret the slice as `complex *16`.

### Constants and intrinsics

| Fortran | C |
|---|---|
| `1.0d0`, `0.5d-2` | `1.0`, `0.5e-2` |
| `dcmplx(a, b)` | `(a) + (b) * I` |
| `dble(z)` | `creal(z)` |
| `imag(z)` | `cimag(z)` |
| `dabs(x)` | `fabs(x)` |
| `cdabs(z)` | `cabs(z)` |
| `dsqrt(x)` | `sqrt(x)` |
| `dlog(x)` | `log(x)` |
| `(0.0d0, 1.0d0)` (eye) | `I` from `<complex.h>` |
| `1/z` (complex) | `1.0 / z` |
| `(-1)**j` | `((j) % 2 == 0 ? 1 : -1)` (don't use `pow`) |

### `implicit real *8 (a-h, o-z)`

A few older Fortran files (notably `sort_pts_to_children`,
`pts_tree_sort`) use Fortran's implicit typing rule: variables whose
names start with letters in `a-h` or `o-z` are implicitly `real *8`
(double); variables starting with `i-n` are implicitly `integer`. When
translating such a file, you must read every undeclared variable and
infer its type from its name.

### Allocations

Fortran:
```fortran
integer, allocatable :: isum(:)
allocate(isum(nbloc))
... use ...
! implicit deallocate at subroutine exit
```

C:
```c
fint *isum = (fint *)malloc(nbloc_v * sizeof(fint));
... use ...
free(isum);
```

For arrays with 0-based Fortran indices: `allocate(arr(0:n))` →
`malloc((n + 1) * sizeof(...))`.

Each Fortran routine has exactly one return point in this codebase, so
all `free` calls go in a single cleanup block at the end of the C
function. Be careful to free *every* malloc; the diff tests don't
catch leaks but a long-running process would.

### `cpu_time` and timing

Many Fortran routines call `cpu_time` and `omp_get_wtime` to measure
section timings, then store them in a `timeinfo(8)` output array.
These are **non-deterministic** — they depend on wall clock — and
must not be in the C translation. The convention adopted is:

- Don't call `cpu_time` or `omp_get_wtime` at all.
- Set `time1 = time2 = 0.0`.
- Compute `timeinfo[i] = time2 - time1`, which stores 0.0.
- The diff test driver explicitly does **not** compare `timeinfo`.

### `prinf`, `prin2`, `prini` (logging)

All wrapped in `if (ifprint .ge. 1)` guards in the Fortran. We
hardcode `ifprint = 0`, so these are dead code. Strip them entirely
from the C version.

### OpenMP pragmas

All `c$omp`, `c$OMP`, `C$OMP$`, and `c$` lines are Fortran comments
to the non-OpenMP build, and the parent fmm2d makefile is invoked
with `OMP=OFF` for this work. Strip them entirely from the C
translation. Do **not** insert `#pragma omp` directives.

### Goto preservation

Translate `goto 1233` and `1233 continue` literally to C `goto skip;`
and `skip: ;`. Don't try to restructure into structured control
flow — the goto is part of the spec and a 1:1 translation is the only
defensible thing.

---

## 4. The differential test architecture

Each translated file has a Fortran 77 fixed-form test driver in
[`test/test_<name>.f`](test/). The driver:

1. Constructs deterministic random inputs (via `hkrand`, seeded
   once with `hkrand(1234)`).
2. For each test case, prepares two output buffers (one for the
   Fortran reference, one for the C translation).
3. Calls the Fortran reference (`<name>`, no underscore — gfortran
   adds it) with one buffer.
4. Calls the C translation (`<name>_c`) with the other buffer.
5. Compares the buffers element-by-element. For integer outputs:
   exact equality. For floating-point outputs: max absolute
   difference, with the goal being **exactly zero**.
6. Tracks `nfail`, prints one line per case, `stop 1` if anything
   failed.

### Why -O0 on both sides

The diff tests build **both** the C translation and the corresponding
Fortran original at `-O0`. This is non-negotiable for bit-for-bit
diff testing. At `-O3`:

- gcc does loop vectorization and may use FMA instructions
  (`vfmadd...`) to fuse `a*b + c` into one rounding step.
- gfortran does the same but with subtly different choices.
- C99 `_Complex` division goes through different code paths than
  gfortran's complex division.

The result: at -O3 the C and Fortran versions produce results that
differ by 1-3 ULPs in the last bit. This is harmless for the FMM
(precision target ~1e-6) but defeats the diff test (tolerance 0).

At -O0, both compilers do exactly what the source code requests, in
the order it requests it. The C translation, written 1:1 against the
Fortran source, then produces bit-identical results.

Importantly: the -O0 constraint only applies to the **diff test
build**. The drop-in build (`make dropin`) re-compiles at -O3 with
`-DFMM2D_DROP_IN` for production speed.

### How -O0 references are produced

The Makefile has explicit per-file rules:

```make
FORT_SRC_cumsum = $(ROOT)/src/common/cumsum.f
FFLAGS_REF      = -O0 -g -std=legacy -fPIC -w

$(BUILD)/cumsum_ref.o: $(FORT_SRC_cumsum) | $(BUILD)
	$(FC) $(FFLAGS_REF) -c $< -o $@
```

The diff test for cumsum links:

```
test/test_cumsum.f
    + build/cumsum.o          (the C translation, also at -O0)
    + build/cumsum_ref.o      (the Fortran original recompiled at -O0)
    + libfmm2d.a              (everything else, at -O3)
```

The `_ref.o` is named on the command line **before** `libfmm2d.a`,
so the linker uses our deterministic copy of the routine under test
instead of the -O3 copy in the .a. Cross-file dependencies (e.g.
`hkrand_`, or other Fortran routines that the routine under test
happens to call) get resolved against `libfmm2d.a` at -O3, which is
fine because they are not under test.

### Critical: static pattern rule for diff tests

GNU Make rejects an ordinary pattern rule for the test binaries with
"Avoiding implicit rule recursion" because the rule chains through
intermediate `%.o` files. The Makefile uses a **static pattern rule**
restricted to the explicit list of test binaries:

```make
$(TEST_BINS): $(BUILD)/test_%: $(TESTDIR)/test_%.f $(BUILD)/%.o $(BUILD)/%_ref.o $(FMMLIB) | $(BUILD)
	$(FC) $(FFLAGS_REF) $< $(BUILD)/$*.o $(BUILD)/$*_ref.o \
	    $(TESTOBJS) $(FMMLIB) -lm -o $@
```

If you remove the `$(TEST_BINS):` prefix and convert it to a regular
pattern rule, Make will silently fall back to its built-in `%: %.f`
rule and you will get cryptic "no rule" errors. Don't do that.

---

## 5. Build system mechanics

### Targets in `Makefile`

- `all` (default) — builds every `.c` in `src/` to a `.o` at -O0,
  every `_ref.o` (the Fortran reference at -O0), and every diff
  test binary.
- `tests` — runs every diff test binary in sequence and prints
  PASS/FAIL.
- `clean` — `rm -rf build/`.
- `dropin` — builds every `.c` at `-O3 -DFMM2D_DROP_IN`, producing
  `build/dropin/*.o` files that export the bare Fortran symbol names.
- `e2e` — builds and runs two end-to-end integration tests against
  the drop-in objects.

### Per-file mappings the Makefile needs

When you add a new translated file `c_translation/src/foo.c`, you need
to also tell the Makefile where the original Fortran lives so it can
build the -O0 reference. Add a line to the `FORT_SRC_*` block:

```make
FORT_SRC_foo = $(ROOT)/src/<subdir>/foo.f
```

And add an explicit reference-build rule:

```make
$(BUILD)/foo_ref.o: $(FORT_SRC_foo) | $(BUILD)
	$(FC) $(FFLAGS_REF) -c $< -o $@
```

(This part could be macro-generated with `$(eval)` but I tried that
and ran into multi-line `define` issues with the recipe getting
dropped. Eleven explicit rules turned out to be simpler.)

### Static lib dependency

The Makefile assumes `../lib-static/libfmm2d.a` exists. Build it from
the repo root with `make lib OMP=OFF`. The diff tests also need
`../src/common/hkrand.o` and `../src/common/dlaran.o` (built when you
run any of the parent makefile's `test/*` targets, or manually).

---

## 6. Recipe: translating one new file

This is the playbook used for files 2 through 11.

### Step 1: Pick the next file from the call graph

For new entry points, first map the call graph. Either by hand
inspection or with an Explore subagent. Identify which Fortran files
contain reachable code, and which routines within each file are
reachable. Translate **only the reachable subset**.

For the existing 11 files, the order was determined by dependency:
leaves first (cumsum, fmmcommon2d, l2dterms, cauchykernels2d), then
operators (laprouts2d, tree_routs2d), then tree builder (pts_tree2d),
then the FMM main (cfmm2d), then the wrappers (cfmm2d_ndiv,
rfmm2d_ndiv).

Bottom-up order means each new file's dependencies are already
trusted, but technically it's not required: cross-file calls go to
the bare Fortran name, so a new file can be translated even before
its dependencies are. The bottom-up order is just for psychological
comfort and easier debugging.

### Step 2: Read the Fortran source carefully

Read every routine you intend to translate, in full. Note:

- The set of subroutines in the file. Which are reachable?
- The argument signatures (especially which args are output).
- Multidimensional array shapes (1-based or 0-based for each dim).
- Calls to other routines (same-file or cross-file).
- Allocations.
- Use of `goto`, complex arithmetic, intrinsics, implicit typing.

### Step 3: Write the C translation

Follow [section 3](#3-translation-conventions). The biggest sources
of bugs are:

- **`itree(iptr(N) + offset)` indexing.** This pattern appears
  throughout the tree code. Translate as `itree[iptr[N-1] + offset - 1]`.
  Read every such access twice; off-by-ones here are silent and hard
  to debug.
- **Multi-term Fortran additions.** `pot = pot + a + b` is
  left-associative in Fortran: `(pot + a) + b`. Writing `pot += a + b`
  in C is right-associative: `pot + (a + b)`. **These give different
  last-bit results.** Always split into two `+=` statements.
- **`1/(z/r)` vs `r/z`.** Different operations, different rounding.
  Translate Fortran's `1/(z/rscale)` literally as `1.0 / (z / rscale_v)`,
  not as `rscale_v / z`.
- **`(-1)**j` for integer j.** Use `(j % 2 == 0 ? 1 : -1)`. Don't use
  `pow(-1, j)` (returns float, may round).

### Step 4: Write the header

Standard pattern (see [`include/cumsum.h`](include/cumsum.h) for an
example):

```c
#ifndef FMM2D_FOO_H
#define FMM2D_FOO_H
#include "fmm2d_c.h"

void FNAME(routine_name)(<args>);

#endif
```

Don't redeclare cross-file dependencies; their headers should be
included by callers separately.

### Step 5: Write the diff test driver

Pattern (see [`test/test_cumsum.f`](test/test_cumsum.f) for the
simplest example, [`test/test_laprouts2d.f`](test/test_laprouts2d.f)
for floating-point):

1. Declare `external` for both `<name>` and `<name>_c` for each
   translated routine, plus `external hkrand`.
2. Allocate input arrays with realistic sizes.
3. Initialize `hkrand` once with `dummy = hkrand(1234)`.
4. Fill inputs with `hkrand(0)` calls.
5. For each test case, prepare two output buffers (Fortran and C).
6. Zero the buffers if the routine *increments* its outputs.
7. Call the Fortran reference into the first buffer.
8. Call the C translation into the second buffer.
9. Compare element-by-element. For integers: `errcount = errcount + 1`
   on any mismatch. For floats: track `errmax` as max abs difference.
10. Use a `report(name, errcount, nfail)` or `report(name, errmax,
    nfail)` helper to print one line per test.
11. `stop 1` at the end if `nfail > 0`.

#### Test buffer leading dimensions

If the routine takes a leading dimension `nd` as an argument and your
test sweeps multiple `nd` values, you must allocate fresh buffers per
sweep iteration with the correct shape. **Do not** declare a fixed
`buf(ndmax, n)` and pass it for both nd=1 and nd=3 — when nd=1 the
routine writes at stride 1 while you read back at stride ndmax,
producing junk. The fix is `allocate(buf(nd, n))` inside the sweep
loop.

#### Fortran 72-column limit

Fortran 77 fixed-form has a 72-column line limit; anything past that
gets silently truncated. Long variable declarations like
`real *8, allocatable :: foo_a(:, :), foo_b(:, :)` can easily exceed
72 columns. When this happens you get cryptic "Function 'foo_b' has no
implicit type" errors. Split into two declarations.

```bash
# Quick check for over-length lines:
awk 'length($0) > 72 {print NR": "length($0)}' test/test_<name>.f
```

### Step 6: Update the Makefile

Add entries for the new file:

1. Add to `FORT_SRC_*` block.
2. Add an explicit `_ref.o` rule.

### Step 7: Build and run

```bash
cd c_translation
make             # build the new C object, ref object, and test driver
build/test_<name>
```

If the test fails:

- **For integer outputs**, the bug is almost certainly an off-by-one
  in array indexing. Bisect by adding `print` statements (or `printf`
  in the C version) to narrow it down.
- **For floating-point outputs**, the bug is almost certainly:
  - A multi-term addition that wasn't split (see [section 7](#7-pitfalls-everything-that-broke-during-the-first-11-files)).
  - A test buffer leading-dimension mismatch.
  - An algorithmic typo (`-zinv*zinv` vs `zinv*zinv`).
- **For build errors**, usually a missing `extern` for a cross-file
  call, or a header guard typo.

When the test prints `PASS: ...`, you're done with that file.

### Step 8: Move to the next file

That's it. The whole loop is mechanical once you've internalized
sections 2 and 3.

---

## 7. Pitfalls (everything that broke during the first 11 files)

Documented here so you don't have to rediscover them.

### 7.1 Fortran left-to-right addition

**Symptom:** `c2d_directcdg`, `c2d_directcdh`, and `c2d_directcdp`
diff tests failed at -O0 with errors at the 1e-12 to 1e-9 level.

**Cause:** The Fortran source was

```fortran
pot(ii,j) = pot(ii,j) + rtmp*charge(ii,i)
   1                  + zinv*dipstr(ii,i)
```

which evaluates left-to-right as `((pot + rtmp*charge) + zinv*dipstr)`.
The first C translation wrote it as

```c
pot[idx] += rtmp * charge[idx] + zinv * dipstr[idx];
```

which parses as `pot + (rtmp*charge + zinv*dipstr)` — a different
parenthesization, different last bit.

**Fix:** Split into two statements:

```c
pot[idx] += rtmp * charge[idx];
pot[idx] += zinv * dipstr[idx];
```

**Generalization:** Any multi-term Fortran assignment `x = x + a + b + c`
must become a sequence of `+=` statements in the same order:

```c
x += a;
x += b;
x += c;
```

This is now baked into the agent prompt template.

### 7.2 Test buffer leading-dimension mismatch

**Symptom:** Same `c2d_direct*` failures as above, but compounded.

**Cause:** The first test driver declared

```fortran
parameter (ndmax = 3)
complex *16 pot_f(ndmax, ntmax)
```

and called the routine with `nd = 1`. The routine treated `pot_f` as
having stride 1 in the first dim and wrote contiguous elements 0..nt-1.
The test then read `pot_f(1, j)` for `j = 1..nt` at offsets `0, 3, 6,
...`, getting different values.

**Fix:** Allocate per-iteration with the actual leading dim:

```fortran
do ind = 1, 2
   nd = nds(ind)
   allocate(pot_f(nd, nt), pot_c(nd, nt))
   ...
   deallocate(pot_f, pot_c)
enddo
```

### 7.3 -O3 numerical drift

**Symptom:** All cauchykernels diff tests reported `errmax` ~1e-13
even after fixing 7.1 and 7.2.

**Cause:** gcc -O3 and gfortran -O3 don't agree bit-for-bit on
complex arithmetic. Confirmed by an isolated repro: a single tiny
loop that computes the same complex multiply-add gives identical
results at -O0 but diverges at -O3.

**Fix:** Build both sides at -O0 for the diff test. The Makefile
infrastructure is set up to do this via the `_ref.o` rules and
`CFLAGS = $(CFLAGS_TEST)`.

### 7.4 `cpu_time` non-determinism

**Symptom:** The `timeinfo(8)` output array of cfmm2d differed between
runs (and obviously between Fortran and C translations).

**Cause:** `cpu_time(time1)` measures wall-clock time. Different
between runs.

**Fix:** Strip all `cpu_time` calls from the C translation. Set
`time1 = time2 = 0.0`. Set `timeinfo[i] = 0.0`. The diff test driver
explicitly does not compare `timeinfo`.

### 7.5 GNU Make implicit rule recursion

**Symptom:** "No rule to make target 'build/test_cauchykernels2d',
needed by 'all'." Even though there was a pattern rule:

```make
$(BUILD)/test_%: $(TESTDIR)/test_%.f $(BUILD)/%.o $(BUILD)/%_ref.o ...
```

**Cause:** GNU Make rejects pattern rules that chain through
intermediate `%.o` files when both rules are pattern rules. It logs
"Avoiding implicit rule recursion."

**Fix:** Convert to a static pattern rule restricted to the explicit
list of test binaries:

```make
$(TEST_BINS): $(BUILD)/test_%: $(TESTDIR)/test_%.f $(BUILD)/%.o $(BUILD)/%_ref.o ...
```

The `$(TEST_BINS):` prefix tells Make "these specific targets follow
this pattern" rather than "any target matching this pattern."

### 7.6 Fortran 72-column limit

**Symptom:** "Function 'reorg_centers_c' has no implicit type" — but
the variable IS declared.

**Cause:** The declaration line was 74 columns long. Fortran 77
fixed-form silently truncates at column 72.

**Fix:** Split long declarations across two lines.

### 7.7 `ar rcs` overwrites .a contents

**Symptom:** `make test/rfmm2d` from the repo root suddenly couldn't
find symbols in `lib-static/libfmm2d.a`.

**Cause:** I had run `ar rcs libfmm2d.a foo.o` thinking it would
*add* foo.o to the archive. It actually *replaces* the archive
contents with just foo.o.

**Fix:** Don't do this. Always rebuild via `make lib` if you need to
modify the .a.

### 7.8 Stray character at start of file (auto-generated)

**Symptom:** A subagent wrote a header file that began with `1/*`
instead of `/*`. Compile failed with "expected identifier or '('
before numeric constant."

**Cause:** Likely a bash heredoc / sed escape interaction.

**Fix:** Read the file and fix the first line. (This happened once in
the entire 11-file translation; if you see weird first-character
errors, check the literal bytes of the file.)

### 7.9 `iaddr` is the most error-prone array

**Symptom:** None observed — but there's a high pre-test risk.

`iaddr(2, nboxes)` is a Fortran-1-based pointer into the workspace
array `rmlexp`. The C translation needs to:

1. Access `iaddr` itself with `FA2(k, ibox, 2)`.
2. Subtract 1 from the value retrieved (since the value is a 1-based
   Fortran index).
3. Cast the resulting pointer to `fcomplex *` for the laprouts2d
   routines.

The full pattern, from `cfmm2d.c`:

```c
(fcomplex *)&rmlexp[iaddr[FA2(1, ibox, 2)] - 1]
```

Read this twice every time you write it.

### 7.10 Cross-file calls to functions translated in the *same* commit

**Note for parallel translation:** if you're translating multiple
files in parallel (e.g., dispatching subagents to do
`tree_routs2d` and `pts_tree2d` simultaneously), this works fine for
the *translation* (each agent sees only the original Fortran source)
but NOT for the *diff test*. Each diff test isolates one file by
calling all cross-file dependencies via the bare Fortran name (which
resolves to libfmm2d.a). So as long as the Fortran library exists,
you can translate and test files independently and in any order.

---

## 8. Subagent prompt template

This is the template that worked for 11 files. Copy and adapt for new
files. The most important sections are the "CRITICAL TRANSLATION
RULES" and the "Pitfalls" section embedded in the prompt.

````
You are translating a single Fortran source file from the fmm2d 2D
Fast Multipole Method library to C, as part of a multi-file porting
effort. STRICT 1:1 translation — same algorithm, same control flow,
same operation order. Your job is research + writing code, not
running the build.

# Repo root
`/home/magland/src/fmm2d`

# Source to translate
`<full path to .f file>`

# Routines to translate (translate ALL of these — they are reachable)
- routine_a(args) — short description
- routine_b(args) — short description
- ...

(Skip these unreachable routines: <list>.)

# Required output files
- `c_translation/src/<name>.c`
- `c_translation/include/<name>.h`
- `c_translation/test/test_<name>.f`

# READ THESE FIRST (mandatory pattern)
- `/home/magland/src/fmm2d/c_translation/include/fmm2d_c.h` — defines FNAME, fint, fcomplex, FA2, FA3
- `/home/magland/src/fmm2d/c_translation/src/<closest analog>.c` — pattern reference
- `/home/magland/src/fmm2d/c_translation/test/test_<closest analog>.f` — test pattern
- `<full path to .f file>` — the FULL source. Read every line.

# CRITICAL TRANSLATION RULES

## Rule 1: Symbol naming
Define each routine via `void FNAME(routine_name)(...)`. All args by
pointer. FNAME defaults to `name##_c_` (diff tests) or `name##_` with
`-DFMM2D_DROP_IN`.

## Rule 2: Cross-file calls
For routines in OTHER source files, use the bare Fortran symbol name
(`other_routine_`) with an `extern` forward declaration. Get exact
signatures from `c_translation/include/<other>.h`. Do NOT use
`FNAME(...)` for cross-file calls.

## Rule 3: Same-file calls
For routines in THE SAME source file, use `FNAME(other_routine)`.

## Rule 4: Floating-point operation order MUST match Fortran exactly
The diff tests build BOTH sides at -O0 for bit-for-bit equality.
This means the C source code must do the same operations in the same
order as the Fortran source.

**Multi-term additions:** Fortran `x = x + a + b` is left-to-right:
`(x + a) + b`. The C compound `x += a + b` parses as `x + (a + b)` —
a different rounding! ALWAYS split multi-term Fortran assignments
into separate `+=` statements in the SAME left-to-right order.

```c
// Fortran: x = x + a + b
x += a;     // not x += a + b;
x += b;
```

## Rule 5: Types
- `integer` → `fint *`
- `real *8` / `double precision` → `double *`
- `complex *16` → `fcomplex *`

## Rule 6: Column-major arrays
Access multidim arrays via FA2/FA3 macros, keeping loop variables
1-based to mirror Fortran. For 0-based dims (`mpole(nd, 0:nterms)`,
`carray(0:ldc, 0:ldc)`, `laddr(2, 0:nlevels)`), define per-file helper
macros.

## Rule 7: Slice arguments
Fortran `arr(i, j)` passed as a slice → C `&arr[FA2(i, j, ld)]`.
Pattern: `array(first_idx, second_idx)` becomes `&array[<linear offset>]`.

## Rule 8: NO OpenMP, NO error handling, NO logging
Strip all `c$omp`, `c$OMP`, `c$` lines. Strip `prinf`/`prin2`/`prini`
calls. Don't add validation. Don't add error checks. Don't add assertions.

## Rule 9: NO cpu_time
Replace `time1 = cpu_time(...)` with `time1 = 0.0;`. Same for
`omp_get_wtime` and `second()`. The diff test driver does NOT
compare `timeinfo`.

## Rule 10: Allocations
`allocate(x(n))` → `malloc(n_v * sizeof(...))`. `allocate(x(0:n))` →
`malloc((n_v + 1) * sizeof(...))`. Free everything before return.

## Rule 11: Goto preservation
Translate `goto N` as C `goto label_N;` and `N continue` as
`label_N: ;`. Don't restructure.

# Header file

Standard pattern with header guard, `#include "fmm2d_c.h"`, and one
prototype per routine via `FNAME(...)`.

# Diff test driver

Modeled on `c_translation/test/test_<closest analog>.f`. For each
routine:
- Allocate two output buffers at the actual per-call leading dimension
  (NOT a fixed maximum).
- Generate randomized inputs via `hkrand` (initialize with `dummy =
  hkrand(1234)`, then `hkrand(0)` per element).
- Call the Fortran reference and the C version on identical inputs.
- Compare element-by-element. Tolerance 1e-15 for FP, exact equality
  for integers. The goal is errmax = 0.
- Track `nfail`, print one line per case, `stop 1` on any failure.

# When done
Report (1) the routines translated, (2) any non-trivial decisions,
(3) the test cases. Do NOT build or run anything yourself.
````

After the agent returns, the user (or you) runs:

```bash
cd c_translation
make 2>&1 | tail -30           # build
build/test_<name>               # run the diff test
```

If it fails, debug with the pitfalls in [section 7](#7-pitfalls-everything-that-broke-during-the-first-11-files)
as your first hypothesis pool. Most failures are one of the eleven
listed there.

---

## 9. Suggested next files

To add another entry point to the C library, the workflow is:

1. **Map the call graph** of the new entry point. Use an Explore
   subagent or grep manually. Identify:
   - Which Fortran source files contain reachable code.
   - Which routines within each file are reachable.
   - Which routines are *already* translated in `c_translation/src/`
     (these are free dependencies).
   - Which Fortran files would need to be translated.

2. **Translate bottom-up.** Leaves first, wrappers last.

3. **For each new file**, follow the recipe in [section 6](#6-recipe-translating-one-new-file).

### Remaining candidates and their dependencies

The 2D fmm2d library has these top-level entry points still to do:

| Entry point | Description | Estimated new files |
|---|---|---|
| `hfmm2d` | Helmholtz FMM | 8-10 files (hfmm2d.f, hfmm2d_ndiv.f, helmrouts2d.f, h2dterms.f, h2dcommon.f, helmkernels2d.f, wideband2d.f, hfmm2d_mps.f, hank103.f) — substantial |
| `mbhfmm2d` | Modified biharmonic FMM | 4-6 files (mbhfmm2d.f, mbhrouts2d.f, mbhkernels2d.f, mbhgreen2d.f) |

The Helmholtz FMM is the biggest single addition. It introduces
Hankel function evaluation (`hank103.f`) which is a standalone
~1000-line numerical-special-functions file.

### What's reusable

The following are already translated and will be free dependencies
for any new entry point:

- All of `cumsum.f`, `fmmcommon2d.f`, `tree_routs2d.f`, `pts_tree2d.f`
  (the tree management infrastructure)
- `l2dterms.f` (Laplace term-count)
- `cauchykernels2d.f`, `laprouts2d.f` (Laplace kernels and operators)
- `cfmm2d.f`, `cfmm2d_ndiv.f`, `rfmm2d_ndiv.f`, `lfmm2d_ndiv.f`
  (Laplace pipeline — real, complex Cauchy, complex log)
- `bhndiv2d.f`, `bh2dterms.f`, `bhkernels2d.f`, `bhrouts2d.f`,
  `bhfmm2d.f` (the biharmonic pipeline)
- `stfmm2d.f` (Stokes wrapper around biharmonic)

For Helmholtz, you'd need to translate the `helm*` and `hank103.f`
files but the tree code is shared.

For modified-biharmonic, the tree code is free; only the
kernel-specific files need translation.

---

## 10. Quick reference: where things live

### In `c_translation/`

- [`Makefile`](Makefile) — standalone build, diff tests, drop-in build, e2e
- [`include/fmm2d_c.h`](include/fmm2d_c.h) — `FNAME`, `fint`, `fcomplex`, `FA2`/`FA3`/`FA4` macros
- [`include/<name>.h`](include/) — one header per translated file
- [`src/<name>.c`](src/) — the 11 translated source files
- [`test/test_<name>.f`](test/) — one diff test per translated file
- [`e2e/test_rfmm2d_e2e.f`](e2e/test_rfmm2d_e2e.f) — matlab-path smoke test
- `build/` — gitignore-able output

### In the parent fmm2d repo

- [`src/laplace/`](../src/laplace/) — Laplace FMM Fortran sources (the main target)
- [`src/common/`](../src/common/) — tree code, prefix sums, FFTs, randoms (mostly translated)
- [`src/helmholtz/`](../src/helmholtz/) — Helmholtz FMM (1 file translated: hndiv2d)
- [`src/biharmonic/`](../src/biharmonic/), [`src/stokes/`](../src/stokes/), [`src/modified-biharmonic/`](../src/modified-biharmonic/) — other entry points (untranslated)
- [`test/laplace/test_rfmm2d.f`](../test/laplace/test_rfmm2d.f) — the existing test that the e2e build runs through the C drop-in
- [`matlab/rfmm2d.m`](../matlab/rfmm2d.m) — the MATLAB entry point this work was scoped to support
- [`lib-static/libfmm2d.a`](../lib-static/) — the parent Fortran library (built by `make lib OMP=OFF` from the repo root). The diff tests link against it for the trusted reference.

### Key things the diff tests need from the parent

- `../lib-static/libfmm2d.a` — the Fortran library at -O3
- `../src/common/hkrand.o` — random number helper (built by parent makefile's test targets)
- `../src/common/dlaran.o` — random number helper

If you cleaned the parent repo and these are missing, rebuild with:

```bash
cd ..
make lib OMP=OFF
gfortran -fPIC -O3 -march=native -funroll-loops -std=legacy -w -c src/common/hkrand.f -o src/common/hkrand.o
gfortran -fPIC -O3 -march=native -funroll-loops -std=legacy -w -c src/common/dlaran.f -o src/common/dlaran.o
cd c_translation
```

---

## Appendix: complete diff test status

Last verified: all 11 diff tests pass at -O0 with errmax = 0.

```
$ make tests
==> build/test_cauchykernels2d
 [ ok ] c2d_directcp   nd=1 ns= 50 nt= 40  errmax=  0.000E+00
 ... 18 cases ...
 PASS: all cauchykernels2d cases match
==> build/test_cfmm2d
 [ ok ] cfmm2d  ifc=1 ifd=0 ifpgh=1 ifpght=1  errmax= 0.000E+00
 ... 9 cases + l2dmpalloc ...
 PASS: all cfmm2d cases match
==> build/test_cfmm2d_ndiv
 [ ok ] cfmm2d_ndiv ifc=1 ifd=0 ifpgh=1 ifpght=1  errmax= 0.000E+00
 ... 9 cases ...
 PASS: all cfmm2d_ndiv cases match
==> build/test_cumsum               PASS (12 cases)
==> build/test_fmmcommon2d          PASS (11 cases)
==> build/test_hndiv2d              PASS (12 cases)
==> build/test_l2dterms             PASS (8 cases)
==> build/test_laprouts2d           PASS (16 routines)
==> build/test_pts_tree2d           PASS (6 routines)
==> build/test_rfmm2d_ndiv          PASS (9 cases)
==> build/test_tree_routs2d         PASS (7 routines)

$ make e2e
==> build/end2end_rfmm2d (existing test_rfmm2d.f via C drop-in)
    PASS (output in build/end2end_rfmm2d.log)
==> build/end2end_rfmm2d_ndiv (matlab path smoke test)
 charge potential at targets: eps= 1.000E-06 max abs err= 5.623E-09 rel err= 5.700E-10
 PASS: rfmm2d_ndiv end-to-end smoke test
```

If you start translating new files and any of these regress, the new
file is the prime suspect — the existing 11 are stable.
