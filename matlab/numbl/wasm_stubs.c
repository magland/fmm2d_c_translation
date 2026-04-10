/*
 * wasm_stubs.c
 *
 * Stub implementations of the Fortran-ABI fmm2d entry points that are not
 * (yet) ported in c_translation/.  Linked into the wasm build only — the
 * native build links against the full Fortran libfmm2d.a, which has real
 * implementations.
 *
 * Each stub calls mexErrMsgTxt to raise a clear "not implemented" error
 * that propagates through mex_dispatch as a RuntimeError in JS. This
 * ensures callers get an immediate, obvious failure rather than silently
 * receiving zero results with ier=1.
 *
 * Signatures match the MWF77_* declarations at the top of matlab/fmm2d.c
 * exactly so wasm-ld accepts the resulting object alongside fmm2d.o.
 */

#include <complex.h>
#include <stdint.h>

typedef _Complex double dcomplex;

/* mexErrMsgTxt is provided by mex_shim.cpp; it longjmps back to
 * mex_dispatch which returns rc=1 and surfaces the message to JS. */
extern void mexErrMsgTxt(const char *msg);

#define STUB_NOT_IMPLEMENTED(name) \
    mexErrMsgTxt("fmm2d wasm: " name " is not yet ported to C/wasm")

/* ── Helmholtz stubs removed: real implementations now in
      c_translation/src/hfmm2d_ndiv.c and helmkernels2d.c ────────── */

/* ── FFT routines removed: real implementations now in
      c_translation/src/dfft_threadsafe.c ───────────────────────── */

/* ── Laplace direct evaluators (l2d_direct*, r2d_direct*) removed:
      real implementations now in c_translation/src/lapkernels2d.c
      and c_translation/src/rlapkernels2d.c ─────────────────────── */

/* ── Stokes direct evaluators (st2ddirect*) removed: real
      implementations now in c_translation/src/stokkernels2d.c ──── */
