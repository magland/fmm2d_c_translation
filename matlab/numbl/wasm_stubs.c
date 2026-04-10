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

/* ── Helmholtz FMM driver ─────────────────────────────────────────── */

void hfmm2d_ndiv_(int *nd, double *eps, dcomplex *zk,
                  int *ifcharge, double *sources,
                  int *ns, dcomplex *charge,
                  int *ifdipole, dcomplex *dipstr, double *dipvec,
                  int *iper, int *ifpgh,
                  dcomplex *pot, dcomplex *grad, dcomplex *hess,
                  int *nt, double *targ, int *ifpghtarg,
                  dcomplex *pottarg, dcomplex *gradtarg, dcomplex *hesstarg,
                  int *ndiv, int *idivflag, int *ifnear,
                  double *timeinfo, int *ier) {
    STUB_NOT_IMPLEMENTED("hfmm2d_ndiv");
}

/* ── Helmholtz direct evaluators (h2d_direct*) ────────────────────── */

void h2d_directcp_(int *nd, dcomplex *zk, double *sources,
                   int *ns, dcomplex *charge, double *targ,
                   int *nt, dcomplex *pot, double *thresh) {
    STUB_NOT_IMPLEMENTED("h2d_directcp");
}
void h2d_directdp_(int *nd, dcomplex *zk, double *sources,
                   int *ns, dcomplex *dipstr, double *dipvec,
                   double *targ, int *nt, dcomplex *pot, double *thresh) {
    STUB_NOT_IMPLEMENTED("h2d_directdp");
}
void h2d_directcdp_(int *nd, dcomplex *zk, double *sources,
                    int *ns, dcomplex *charge, dcomplex *dipstr,
                    double *dipvec, double *targ, int *nt,
                    dcomplex *pot, double *thresh) {
    STUB_NOT_IMPLEMENTED("h2d_directcdp");
}
void h2d_directcg_(int *nd, dcomplex *zk, double *sources,
                   int *ns, dcomplex *charge, double *targ, int *nt,
                   dcomplex *pot, dcomplex *grad, double *thresh) {
    STUB_NOT_IMPLEMENTED("h2d_directcg");
}
void h2d_directdg_(int *nd, dcomplex *zk, double *sources,
                   int *ns, dcomplex *dipstr, double *dipvec,
                   double *targ, int *nt,
                   dcomplex *pot, dcomplex *grad, double *thresh) {
    STUB_NOT_IMPLEMENTED("h2d_directdg");
}
void h2d_directcdg_(int *nd, dcomplex *zk, double *sources,
                    int *ns, dcomplex *charge, dcomplex *dipstr,
                    double *dipvec, double *targ, int *nt,
                    dcomplex *pot, dcomplex *grad, double *thresh) {
    STUB_NOT_IMPLEMENTED("h2d_directcdg");
}
void h2d_directch_(int *nd, dcomplex *zk, double *sources,
                   int *ns, dcomplex *charge, double *targ, int *nt,
                   dcomplex *pot, dcomplex *grad, dcomplex *hess,
                   double *thresh) {
    STUB_NOT_IMPLEMENTED("h2d_directch");
}
void h2d_directdh_(int *nd, dcomplex *zk, double *sources,
                   int *ns, dcomplex *dipstr, double *dipvec,
                   double *targ, int *nt,
                   dcomplex *pot, dcomplex *grad, dcomplex *hess,
                   double *thresh) {
    STUB_NOT_IMPLEMENTED("h2d_directdh");
}
void h2d_directcdh_(int *nd, dcomplex *zk, double *sources,
                    int *ns, dcomplex *charge, dcomplex *dipstr,
                    double *dipvec, double *targ, int *nt,
                    dcomplex *pot, dcomplex *grad, dcomplex *hess,
                    double *thresh) {
    STUB_NOT_IMPLEMENTED("h2d_directcdh");
}

/* ── Laplace direct evaluators (l2d_direct*) — complex out ────────── */

void l2d_directcp_(int *nd, double *sources, int *ns,
                   dcomplex *charge, double *targ, int *nt,
                   dcomplex *pot, double *thresh) {
    STUB_NOT_IMPLEMENTED("l2d_directcp");
}
void l2d_directdp_(int *nd, double *sources, int *ns,
                   dcomplex *dipstr, double *dipvec,
                   double *targ, int *nt, dcomplex *pot, double *thresh) {
    STUB_NOT_IMPLEMENTED("l2d_directdp");
}
void l2d_directcdp_(int *nd, double *sources, int *ns,
                    dcomplex *charge, dcomplex *dipstr, double *dipvec,
                    double *targ, int *nt, dcomplex *pot, double *thresh) {
    STUB_NOT_IMPLEMENTED("l2d_directcdp");
}
void l2d_directcg_(int *nd, double *sources, int *ns,
                   dcomplex *charge, double *targ, int *nt,
                   dcomplex *pot, dcomplex *grad, double *thresh) {
    STUB_NOT_IMPLEMENTED("l2d_directcg");
}
void l2d_directdg_(int *nd, double *sources, int *ns,
                   dcomplex *dipstr, double *dipvec, double *targ, int *nt,
                   dcomplex *pot, dcomplex *grad, double *thresh) {
    STUB_NOT_IMPLEMENTED("l2d_directdg");
}
void l2d_directcdg_(int *nd, double *sources, int *ns,
                    dcomplex *charge, dcomplex *dipstr, double *dipvec,
                    double *targ, int *nt,
                    dcomplex *pot, dcomplex *grad, double *thresh) {
    STUB_NOT_IMPLEMENTED("l2d_directcdg");
}
void l2d_directch_(int *nd, double *sources, int *ns,
                   dcomplex *charge, double *targ, int *nt,
                   dcomplex *pot, dcomplex *grad, dcomplex *hess,
                   double *thresh) {
    STUB_NOT_IMPLEMENTED("l2d_directch");
}
void l2d_directdh_(int *nd, double *sources, int *ns,
                   dcomplex *dipstr, double *dipvec, double *targ, int *nt,
                   dcomplex *pot, dcomplex *grad, dcomplex *hess,
                   double *thresh) {
    STUB_NOT_IMPLEMENTED("l2d_directdh");
}
void l2d_directcdh_(int *nd, double *sources, int *ns,
                    dcomplex *charge, dcomplex *dipstr, double *dipvec,
                    double *targ, int *nt,
                    dcomplex *pot, dcomplex *grad, dcomplex *hess,
                    double *thresh) {
    STUB_NOT_IMPLEMENTED("l2d_directcdh");
}

/* ── Real direct evaluators (r2d_direct*) — real out ──────────────── */

void r2d_directcp_(int *nd, double *sources, int *ns,
                   double *charge, double *targ, int *nt,
                   double *pot, double *thresh) {
    STUB_NOT_IMPLEMENTED("r2d_directcp");
}
void r2d_directdp_(int *nd, double *sources, int *ns,
                   double *dipstr, double *dipvec,
                   double *targ, int *nt, double *pot, double *thresh) {
    STUB_NOT_IMPLEMENTED("r2d_directdp");
}
void r2d_directcdp_(int *nd, double *sources, int *ns,
                    double *charge, double *dipstr, double *dipvec,
                    double *targ, int *nt, double *pot, double *thresh) {
    STUB_NOT_IMPLEMENTED("r2d_directcdp");
}
void r2d_directcg_(int *nd, double *sources, int *ns,
                   double *charge, double *targ, int *nt,
                   double *pot, double *grad, double *thresh) {
    STUB_NOT_IMPLEMENTED("r2d_directcg");
}
void r2d_directdg_(int *nd, double *sources, int *ns,
                   double *dipstr, double *dipvec, double *targ, int *nt,
                   double *pot, double *grad, double *thresh) {
    STUB_NOT_IMPLEMENTED("r2d_directdg");
}
void r2d_directcdg_(int *nd, double *sources, int *ns,
                    double *charge, double *dipstr, double *dipvec,
                    double *targ, int *nt,
                    double *pot, double *grad, double *thresh) {
    STUB_NOT_IMPLEMENTED("r2d_directcdg");
}
void r2d_directch_(int *nd, double *sources, int *ns,
                   double *charge, double *targ, int *nt,
                   double *pot, double *grad, double *hess, double *thresh) {
    STUB_NOT_IMPLEMENTED("r2d_directch");
}
void r2d_directdh_(int *nd, double *sources, int *ns,
                   double *dipstr, double *dipvec, double *targ, int *nt,
                   double *pot, double *grad, double *hess, double *thresh) {
    STUB_NOT_IMPLEMENTED("r2d_directdh");
}
void r2d_directcdh_(int *nd, double *sources, int *ns,
                    double *charge, double *dipstr, double *dipvec,
                    double *targ, int *nt,
                    double *pot, double *grad, double *hess, double *thresh) {
    STUB_NOT_IMPLEMENTED("r2d_directcdh");
}

/* ── Stokes direct evaluators (st2d_direct*) ──────────────────────── */

void st2ddirectstokg_(int *nd, double *sources, double *stoklet,
                      int *ns, double *targ, int *nt,
                      double *pot, double *pre, double *grad, double *thresh) {
    STUB_NOT_IMPLEMENTED("st2ddirectstokg");
}
void st2ddirectstokstrsg_(int *nd, double *sources,
                          int *ifstoklet, double *stoklet,
                          int *ifstrslet, double *strslet, double *strsvec,
                          int *ns, double *targ, int *nt,
                          double *pot, double *pre, double *grad,
                          double *thresh) {
    STUB_NOT_IMPLEMENTED("st2ddirectstokstrsg");
}
