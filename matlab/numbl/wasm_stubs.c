/*
 * wasm_stubs.c
 *
 * Stub implementations of the Fortran-ABI fmm2d entry points that are not
 * (yet) ported in c_translation/.  Linked into the wasm build only — the
 * native build links against the full Fortran libfmm2d.a, which has real
 * implementations.
 *
 * Each stub sets `ier` (the trailing int* arg) to a nonzero error code
 * when present, so any caller that actually invokes one via fmm2d.c
 * gets a clean runtime failure rather than crashing on a missing symbol.
 * The wasm backend therefore cannot run the parts of rfmm2dTest /
 * lfmm2dTest / cfmm2dTest / hfmm2dTest / stfmm2dTest that exercise
 * the direct evaluators or the Helmholtz FMM until c_translation grows
 * the missing entry points.  The plain FMM driver paths
 * (rfmm2d_ndiv, lfmm2d_ndiv, cfmm2d_ndiv, stfmm2d) work fine in wasm
 * because c_translation provides them.
 *
 * Signatures match the MWF77_* declarations at the top of matlab/fmm2d.c
 * exactly so wasm-ld accepts the resulting object alongside fmm2d.o.
 */

#include <complex.h>
#include <stdint.h>

typedef _Complex double dcomplex;

/* Set *ier (always the LAST int* arg of these stubs that take one) to a
 * nonzero value.  fmm2d.c's stubs propagate it back through plhs and the
 * MATLAB caller will see ier != 0.  Stubs that don't take an ier just
 * abort silently (they're never invoked through the test path). */
#define STUB_IER(p) do { if (p) *(p) = 1; } while (0)

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
    (void)nd; (void)eps; (void)zk; (void)ifcharge; (void)sources;
    (void)ns; (void)charge; (void)ifdipole; (void)dipstr; (void)dipvec;
    (void)iper; (void)ifpgh; (void)pot; (void)grad; (void)hess;
    (void)nt; (void)targ; (void)ifpghtarg;
    (void)pottarg; (void)gradtarg; (void)hesstarg;
    (void)ndiv; (void)idivflag; (void)ifnear; (void)timeinfo;
    STUB_IER(ier);
}

/* ── Helmholtz direct evaluators (h2d_direct*) ────────────────────── */

void h2d_directcp_(int *nd, dcomplex *zk, double *sources,
                   int *ns, dcomplex *charge, double *targ,
                   int *nt, dcomplex *pot, double *thresh) {
    (void)nd; (void)zk; (void)sources; (void)ns; (void)charge;
    (void)targ; (void)nt; (void)pot; (void)thresh;
}
void h2d_directdp_(int *nd, dcomplex *zk, double *sources,
                   int *ns, dcomplex *dipstr, double *dipvec,
                   double *targ, int *nt, dcomplex *pot, double *thresh) {
    (void)nd; (void)zk; (void)sources; (void)ns; (void)dipstr;
    (void)dipvec; (void)targ; (void)nt; (void)pot; (void)thresh;
}
void h2d_directcdp_(int *nd, dcomplex *zk, double *sources,
                    int *ns, dcomplex *charge, dcomplex *dipstr,
                    double *dipvec, double *targ, int *nt,
                    dcomplex *pot, double *thresh) {
    (void)nd; (void)zk; (void)sources; (void)ns; (void)charge;
    (void)dipstr; (void)dipvec; (void)targ; (void)nt; (void)pot; (void)thresh;
}
void h2d_directcg_(int *nd, dcomplex *zk, double *sources,
                   int *ns, dcomplex *charge, double *targ, int *nt,
                   dcomplex *pot, dcomplex *grad, double *thresh) {
    (void)nd; (void)zk; (void)sources; (void)ns; (void)charge;
    (void)targ; (void)nt; (void)pot; (void)grad; (void)thresh;
}
void h2d_directdg_(int *nd, dcomplex *zk, double *sources,
                   int *ns, dcomplex *dipstr, double *dipvec,
                   double *targ, int *nt,
                   dcomplex *pot, dcomplex *grad, double *thresh) {
    (void)nd; (void)zk; (void)sources; (void)ns; (void)dipstr;
    (void)dipvec; (void)targ; (void)nt; (void)pot; (void)grad; (void)thresh;
}
void h2d_directcdg_(int *nd, dcomplex *zk, double *sources,
                    int *ns, dcomplex *charge, dcomplex *dipstr,
                    double *dipvec, double *targ, int *nt,
                    dcomplex *pot, dcomplex *grad, double *thresh) {
    (void)nd; (void)zk; (void)sources; (void)ns; (void)charge;
    (void)dipstr; (void)dipvec; (void)targ; (void)nt;
    (void)pot; (void)grad; (void)thresh;
}
void h2d_directch_(int *nd, dcomplex *zk, double *sources,
                   int *ns, dcomplex *charge, double *targ, int *nt,
                   dcomplex *pot, dcomplex *grad, dcomplex *hess,
                   double *thresh) {
    (void)nd; (void)zk; (void)sources; (void)ns; (void)charge;
    (void)targ; (void)nt; (void)pot; (void)grad; (void)hess; (void)thresh;
}
void h2d_directdh_(int *nd, dcomplex *zk, double *sources,
                   int *ns, dcomplex *dipstr, double *dipvec,
                   double *targ, int *nt,
                   dcomplex *pot, dcomplex *grad, dcomplex *hess,
                   double *thresh) {
    (void)nd; (void)zk; (void)sources; (void)ns; (void)dipstr;
    (void)dipvec; (void)targ; (void)nt;
    (void)pot; (void)grad; (void)hess; (void)thresh;
}
void h2d_directcdh_(int *nd, dcomplex *zk, double *sources,
                    int *ns, dcomplex *charge, dcomplex *dipstr,
                    double *dipvec, double *targ, int *nt,
                    dcomplex *pot, dcomplex *grad, dcomplex *hess,
                    double *thresh) {
    (void)nd; (void)zk; (void)sources; (void)ns; (void)charge;
    (void)dipstr; (void)dipvec; (void)targ; (void)nt;
    (void)pot; (void)grad; (void)hess; (void)thresh;
}

/* ── Laplace direct evaluators (l2d_direct*) — complex out ────────── */

void l2d_directcp_(int *nd, double *sources, int *ns,
                   dcomplex *charge, double *targ, int *nt,
                   dcomplex *pot, double *thresh) {
    (void)nd; (void)sources; (void)ns; (void)charge;
    (void)targ; (void)nt; (void)pot; (void)thresh;
}
void l2d_directdp_(int *nd, double *sources, int *ns,
                   dcomplex *dipstr, double *dipvec,
                   double *targ, int *nt, dcomplex *pot, double *thresh) {
    (void)nd; (void)sources; (void)ns; (void)dipstr; (void)dipvec;
    (void)targ; (void)nt; (void)pot; (void)thresh;
}
void l2d_directcdp_(int *nd, double *sources, int *ns,
                    dcomplex *charge, dcomplex *dipstr, double *dipvec,
                    double *targ, int *nt, dcomplex *pot, double *thresh) {
    (void)nd; (void)sources; (void)ns; (void)charge; (void)dipstr;
    (void)dipvec; (void)targ; (void)nt; (void)pot; (void)thresh;
}
void l2d_directcg_(int *nd, double *sources, int *ns,
                   dcomplex *charge, double *targ, int *nt,
                   dcomplex *pot, dcomplex *grad, double *thresh) {
    (void)nd; (void)sources; (void)ns; (void)charge;
    (void)targ; (void)nt; (void)pot; (void)grad; (void)thresh;
}
void l2d_directdg_(int *nd, double *sources, int *ns,
                   dcomplex *dipstr, double *dipvec, double *targ, int *nt,
                   dcomplex *pot, dcomplex *grad, double *thresh) {
    (void)nd; (void)sources; (void)ns; (void)dipstr; (void)dipvec;
    (void)targ; (void)nt; (void)pot; (void)grad; (void)thresh;
}
void l2d_directcdg_(int *nd, double *sources, int *ns,
                    dcomplex *charge, dcomplex *dipstr, double *dipvec,
                    double *targ, int *nt,
                    dcomplex *pot, dcomplex *grad, double *thresh) {
    (void)nd; (void)sources; (void)ns; (void)charge; (void)dipstr;
    (void)dipvec; (void)targ; (void)nt;
    (void)pot; (void)grad; (void)thresh;
}
void l2d_directch_(int *nd, double *sources, int *ns,
                   dcomplex *charge, double *targ, int *nt,
                   dcomplex *pot, dcomplex *grad, dcomplex *hess,
                   double *thresh) {
    (void)nd; (void)sources; (void)ns; (void)charge;
    (void)targ; (void)nt; (void)pot; (void)grad; (void)hess; (void)thresh;
}
void l2d_directdh_(int *nd, double *sources, int *ns,
                   dcomplex *dipstr, double *dipvec, double *targ, int *nt,
                   dcomplex *pot, dcomplex *grad, dcomplex *hess,
                   double *thresh) {
    (void)nd; (void)sources; (void)ns; (void)dipstr; (void)dipvec;
    (void)targ; (void)nt;
    (void)pot; (void)grad; (void)hess; (void)thresh;
}
void l2d_directcdh_(int *nd, double *sources, int *ns,
                    dcomplex *charge, dcomplex *dipstr, double *dipvec,
                    double *targ, int *nt,
                    dcomplex *pot, dcomplex *grad, dcomplex *hess,
                    double *thresh) {
    (void)nd; (void)sources; (void)ns; (void)charge; (void)dipstr;
    (void)dipvec; (void)targ; (void)nt;
    (void)pot; (void)grad; (void)hess; (void)thresh;
}

/* ── Real direct evaluators (r2d_direct*) — real out ──────────────── */

void r2d_directcp_(int *nd, double *sources, int *ns,
                   double *charge, double *targ, int *nt,
                   double *pot, double *thresh) {
    (void)nd; (void)sources; (void)ns; (void)charge;
    (void)targ; (void)nt; (void)pot; (void)thresh;
}
void r2d_directdp_(int *nd, double *sources, int *ns,
                   double *dipstr, double *dipvec,
                   double *targ, int *nt, double *pot, double *thresh) {
    (void)nd; (void)sources; (void)ns; (void)dipstr; (void)dipvec;
    (void)targ; (void)nt; (void)pot; (void)thresh;
}
void r2d_directcdp_(int *nd, double *sources, int *ns,
                    double *charge, double *dipstr, double *dipvec,
                    double *targ, int *nt, double *pot, double *thresh) {
    (void)nd; (void)sources; (void)ns; (void)charge; (void)dipstr;
    (void)dipvec; (void)targ; (void)nt; (void)pot; (void)thresh;
}
void r2d_directcg_(int *nd, double *sources, int *ns,
                   double *charge, double *targ, int *nt,
                   double *pot, double *grad, double *thresh) {
    (void)nd; (void)sources; (void)ns; (void)charge;
    (void)targ; (void)nt; (void)pot; (void)grad; (void)thresh;
}
void r2d_directdg_(int *nd, double *sources, int *ns,
                   double *dipstr, double *dipvec, double *targ, int *nt,
                   double *pot, double *grad, double *thresh) {
    (void)nd; (void)sources; (void)ns; (void)dipstr; (void)dipvec;
    (void)targ; (void)nt; (void)pot; (void)grad; (void)thresh;
}
void r2d_directcdg_(int *nd, double *sources, int *ns,
                    double *charge, double *dipstr, double *dipvec,
                    double *targ, int *nt,
                    double *pot, double *grad, double *thresh) {
    (void)nd; (void)sources; (void)ns; (void)charge; (void)dipstr;
    (void)dipvec; (void)targ; (void)nt;
    (void)pot; (void)grad; (void)thresh;
}
void r2d_directch_(int *nd, double *sources, int *ns,
                   double *charge, double *targ, int *nt,
                   double *pot, double *grad, double *hess, double *thresh) {
    (void)nd; (void)sources; (void)ns; (void)charge;
    (void)targ; (void)nt; (void)pot; (void)grad; (void)hess; (void)thresh;
}
void r2d_directdh_(int *nd, double *sources, int *ns,
                   double *dipstr, double *dipvec, double *targ, int *nt,
                   double *pot, double *grad, double *hess, double *thresh) {
    (void)nd; (void)sources; (void)ns; (void)dipstr; (void)dipvec;
    (void)targ; (void)nt;
    (void)pot; (void)grad; (void)hess; (void)thresh;
}
void r2d_directcdh_(int *nd, double *sources, int *ns,
                    double *charge, double *dipstr, double *dipvec,
                    double *targ, int *nt,
                    double *pot, double *grad, double *hess, double *thresh) {
    (void)nd; (void)sources; (void)ns; (void)charge; (void)dipstr;
    (void)dipvec; (void)targ; (void)nt;
    (void)pot; (void)grad; (void)hess; (void)thresh;
}

/* ── Stokes direct evaluators (st2d_direct*) ──────────────────────── */

void st2ddirectstokg_(int *nd, double *sources, double *stoklet,
                      int *ns, double *targ, int *nt,
                      double *pot, double *pre, double *grad, double *thresh) {
    (void)nd; (void)sources; (void)stoklet; (void)ns;
    (void)targ; (void)nt; (void)pot; (void)pre; (void)grad; (void)thresh;
}
void st2ddirectstokstrsg_(int *nd, double *sources,
                          int *ifstoklet, double *stoklet,
                          int *ifstrslet, double *strslet, double *strsvec,
                          int *ns, double *targ, int *nt,
                          double *pot, double *pre, double *grad,
                          double *thresh) {
    (void)nd; (void)sources; (void)ifstoklet; (void)stoklet;
    (void)ifstrslet; (void)strslet; (void)strsvec; (void)ns;
    (void)targ; (void)nt; (void)pot; (void)pre; (void)grad; (void)thresh;
}
