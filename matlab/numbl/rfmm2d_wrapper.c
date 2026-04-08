/*
 * rfmm2d_wrapper.c
 *
 * C wrappers around the Fortran-ABI fmm2d entry points used by numbl,
 * for both the WebAssembly build and the native shared-library build.
 *
 * Currently exports:
 *   rfmm2d_w  - real Laplace FMM (calls rfmm2d_ndiv_)
 *   cfmm2d_w  - complex Cauchy-kernel Laplace FMM (calls cfmm2d_ndiv_)
 *   lfmm2d_w  - complex log-kernel Laplace FMM (calls lfmm2d_ndiv_)
 *   stfmm2d_w - Stokes FMM (calls stfmm2d_)
 *
 * The first three entry points call hndiv2d_ to pick a default
 * subdivision criterion before calling the appropriate _ndiv_ Fortran
 * routine. stfmm2d_w calls stfmm2d_ directly because the Stokes driver
 * does its own tree management via bhfmm2d/bhndiv2d internally.
 * All callees use the gfortran ABI (lowercase + trailing underscore,
 * arguments passed by pointer, fortran integers are int32).
 *
 * Linkage:
 *   WASM    - linked against c_translation/src .c files compiled with
 *             -DFMM2D_DROP_IN, which exports the bare Fortran names.
 *   native  - linked against the top-level libfmm2d.a built by
 *             gfortran, which exports the real Fortran symbols.
 * In both cases the externs declared here resolve to the same names.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#define EXPORT(name) __attribute__((export_name(#name)))
#else
#define EXPORT(name) __attribute__((visibility("default")))
#endif

/* hndiv2d_(eps, ns, nt, ifcharge, ifdipole, ifpgh, ifpghtarg, ndiv, idivflag) */
extern void hndiv2d_(const double *eps,
                     const int32_t *ns, const int32_t *nt,
                     const int32_t *ifcharge, const int32_t *ifdipole,
                     const int32_t *ifpgh, const int32_t *ifpghtarg,
                     int32_t *ndiv, int32_t *idivflag);

/* rfmm2d_ndiv_(nd, eps, ns, sources, ifcharge, charge, ifdipole, dipstr,
 *              dipvec, iper, ifpgh, pot, grad, hess, nt, targ, ifpghtarg,
 *              pottarg, gradtarg, hesstarg, ndiv, idivflag, ifnear,
 *              timeinfo, ier) */
extern void rfmm2d_ndiv_(const int32_t *nd, const double *eps,
                         const int32_t *ns, const double *sources,
                         const int32_t *ifcharge, const double *charge,
                         const int32_t *ifdipole, const double *dipstr,
                         const double *dipvec,
                         int32_t *iper, const int32_t *ifpgh,
                         double *pot, double *grad, double *hess,
                         const int32_t *nt, const double *targ,
                         const int32_t *ifpghtarg,
                         double *pottarg, double *gradtarg,
                         double *hesstarg,
                         const int32_t *ndiv, const int32_t *idivflag,
                         const int32_t *ifnear, double *timeinfo,
                         int32_t *ier);

/*
 * rfmm2d_w
 *
 * Single entry point that the JS shim calls. All array buffers are
 * caller-owned (allocated on the JS / WASM heap by the shim). Output
 * buffers must be sized by the caller; the JS shim sizes them based
 * on nd, ns, nt, ifpgh, ifpghtarg.
 *
 * Buffer shapes (column-major, matching the Fortran callee):
 *   sources : (2, ns)
 *   charge  : (nd, ns)        ignored if ifcharge == 0
 *   dipstr  : (nd, ns)        ignored if ifdipole == 0
 *   dipvec  : (2*nd, ns)      ignored if ifdipole == 0
 *   targ    : (2, nt)         ignored if nt == 0
 *   pot     : (nd, ns)        written if ifpgh >= 1
 *   grad    : (2*nd, ns)      written if ifpgh >= 2
 *   hess    : (3*nd, ns)      written if ifpgh >= 3
 *   pottarg : (nd, nt)        written if ifpghtarg >= 1
 *   gradtarg: (2*nd, nt)      written if ifpghtarg >= 2
 *   hesstarg: (3*nd, nt)      written if ifpghtarg >= 3
 */
EXPORT(rfmm2d_w)
int32_t rfmm2d_w(int32_t nd, double eps, int32_t ns, double *sources,
                 int32_t ifcharge, double *charge,
                 int32_t ifdipole, double *dipstr, double *dipvec,
                 int32_t ifpgh, double *pot, double *grad, double *hess,
                 int32_t nt, double *targ, int32_t ifpghtarg,
                 double *pottarg, double *gradtarg, double *hesstarg)
{
    int32_t ndiv = 20;
    int32_t idivflag = 0;
    int32_t iper = 1;
    int32_t ifnear = 1;
    double timeinfo[8] = {0};
    int32_t ier = 0;

    hndiv2d_(&eps, &ns, &nt, &ifcharge, &ifdipole,
             &ifpgh, &ifpghtarg, &ndiv, &idivflag);

    rfmm2d_ndiv_(&nd, &eps, &ns, sources,
                 &ifcharge, charge, &ifdipole, dipstr, dipvec,
                 &iper, &ifpgh, pot, grad, hess,
                 &nt, targ, &ifpghtarg, pottarg, gradtarg, hesstarg,
                 &ndiv, &idivflag, &ifnear, timeinfo, &ier);

    return ier;
}

/* cfmm2d_ndiv_(nd, eps, ns, sources, ifcharge, charge, ifdipole, dipstr,
 *              iper, ifpgh, pot, grad, hess, nt, targ, ifpghtarg,
 *              pottarg, gradtarg, hesstarg, ndiv, idivflag, ifnear,
 *              timeinfo, ier)
 *
 * Note: cfmm2d_ndiv has NO dipvec argument (dipoles are encoded fully in
 * the complex dipstr) and pot/grad/hess are each (nd, *) complex
 * (i.e. one complex per source/target — grad is d/dz, hess is d^2/dz^2).
 */
extern void cfmm2d_ndiv_(const int32_t *nd, const double *eps,
                         const int32_t *ns, const double *sources,
                         const int32_t *ifcharge, const double *charge,
                         const int32_t *ifdipole, const double *dipstr,
                         int32_t *iper, const int32_t *ifpgh,
                         double *pot, double *grad, double *hess,
                         const int32_t *nt, const double *targ,
                         const int32_t *ifpghtarg,
                         double *pottarg, double *gradtarg,
                         double *hesstarg,
                         const int32_t *ndiv, const int32_t *idivflag,
                         const int32_t *ifnear, double *timeinfo,
                         int32_t *ier);

/*
 * cfmm2d_w
 *
 * Same calling pattern as rfmm2d_w but for the complex-valued Cauchy
 * kernel. Note the buffer-shape differences vs rfmm2d:
 *   - charge / dipstr are complex (each element is a pair of doubles)
 *   - there is NO dipvec
 *   - pot / grad / hess have shape (nd, ns) complex (NOT (2*nd,ns) /
 *     (3*nd,ns) like the real case): grad holds d/dz, hess d^2/dz^2.
 *
 * All buffers are passed as `double *` for ABI simplicity. The JS shim
 * is responsible for laying out complex arrays as interleaved
 * (re, im, re, im, ...) doubles, matching gfortran's `complex *16` layout.
 *
 * Buffer shapes (column-major, doubles per element):
 *   sources : (2, ns)
 *   charge  : (2*nd, ns)      ignored if ifcharge == 0
 *   dipstr  : (2*nd, ns)      ignored if ifdipole == 0
 *   targ    : (2, nt)         ignored if nt == 0
 *   pot     : (2*nd, ns)      written if ifpgh >= 1
 *   grad    : (2*nd, ns)      written if ifpgh >= 2
 *   hess    : (2*nd, ns)      written if ifpgh >= 3
 *   pottarg : (2*nd, nt)      written if ifpghtarg >= 1
 *   gradtarg: (2*nd, nt)      written if ifpghtarg >= 2
 *   hesstarg: (2*nd, nt)      written if ifpghtarg >= 3
 */
EXPORT(cfmm2d_w)
int32_t cfmm2d_w(int32_t nd, double eps, int32_t ns, double *sources,
                 int32_t ifcharge, double *charge,
                 int32_t ifdipole, double *dipstr,
                 int32_t ifpgh, double *pot, double *grad, double *hess,
                 int32_t nt, double *targ, int32_t ifpghtarg,
                 double *pottarg, double *gradtarg, double *hesstarg)
{
    int32_t ndiv = 20;
    int32_t idivflag = 0;
    int32_t iper = 1;
    int32_t ifnear = 1;
    double timeinfo[8] = {0};
    int32_t ier = 0;

    hndiv2d_(&eps, &ns, &nt, &ifcharge, &ifdipole,
             &ifpgh, &ifpghtarg, &ndiv, &idivflag);

    cfmm2d_ndiv_(&nd, &eps, &ns, sources,
                 &ifcharge, charge, &ifdipole, dipstr,
                 &iper, &ifpgh, pot, grad, hess,
                 &nt, targ, &ifpghtarg, pottarg, gradtarg, hesstarg,
                 &ndiv, &idivflag, &ifnear, timeinfo, &ier);

    return ier;
}

/* lfmm2d_ndiv_(nd, eps, ns, sources, ifcharge, charge, ifdipole, dipstr,
 *              dipvec, iper, ifpgh, pot, grad, hess, nt, targ, ifpghtarg,
 *              pottarg, gradtarg, hesstarg, ndiv, idivflag, ifnear,
 *              timeinfo, ier)
 *
 * lfmm2d uses real dipvec (like rfmm2d) but complex charge/dipstr/pot/grad/hess.
 * grad shape is (nd, 2, *), hess shape is (nd, 3, *), all complex.
 */
extern void lfmm2d_ndiv_(const int32_t *nd, const double *eps,
                         const int32_t *ns, const double *sources,
                         const int32_t *ifcharge, const double *charge,
                         const int32_t *ifdipole, const double *dipstr,
                         const double *dipvec,
                         int32_t *iper, const int32_t *ifpgh,
                         double *pot, double *grad, double *hess,
                         const int32_t *nt, const double *targ,
                         const int32_t *ifpghtarg,
                         double *pottarg, double *gradtarg,
                         double *hesstarg,
                         const int32_t *ndiv, const int32_t *idivflag,
                         const int32_t *ifnear, double *timeinfo,
                         int32_t *ier);

/*
 * lfmm2d_w
 *
 * Real-valued log-kernel Laplace FMM with complex densities. Buffer
 * shapes mirror rfmm2d but with complex (2x storage) for the density
 * and output fields:
 *   sources : (2, ns)
 *   charge  : (2*nd, ns)      ignored if ifcharge == 0
 *   dipstr  : (2*nd, ns)      ignored if ifdipole == 0
 *   dipvec  : (2*nd, ns)      ignored if ifdipole == 0  (REAL — same as rfmm2d)
 *   targ    : (2, nt)         ignored if nt == 0
 *   pot     : (2*nd, ns)      written if ifpgh >= 1     (complex)
 *   grad    : (4*nd, ns)      written if ifpgh >= 2     (nd*2 complex)
 *   hess    : (6*nd, ns)      written if ifpgh >= 3     (nd*3 complex)
 *   pottarg : (2*nd, nt)      written if ifpghtarg >= 1
 *   gradtarg: (4*nd, nt)      written if ifpghtarg >= 2
 *   hesstarg: (6*nd, nt)      written if ifpghtarg >= 3
 */
EXPORT(lfmm2d_w)
int32_t lfmm2d_w(int32_t nd, double eps, int32_t ns, double *sources,
                 int32_t ifcharge, double *charge,
                 int32_t ifdipole, double *dipstr, double *dipvec,
                 int32_t ifpgh, double *pot, double *grad, double *hess,
                 int32_t nt, double *targ, int32_t ifpghtarg,
                 double *pottarg, double *gradtarg, double *hesstarg)
{
    int32_t ndiv = 20;
    int32_t idivflag = 0;
    int32_t iper = 1;
    int32_t ifnear = 1;
    double timeinfo[8] = {0};
    int32_t ier = 0;

    hndiv2d_(&eps, &ns, &nt, &ifcharge, &ifdipole,
             &ifpgh, &ifpghtarg, &ndiv, &idivflag);

    lfmm2d_ndiv_(&nd, &eps, &ns, sources,
                 &ifcharge, charge, &ifdipole, dipstr, dipvec,
                 &iper, &ifpgh, pot, grad, hess,
                 &nt, targ, &ifpghtarg, pottarg, gradtarg, hesstarg,
                 &ndiv, &idivflag, &ifnear, timeinfo, &ier);

    return ier;
}

/* stfmm2d_(nd, eps, nsource, source, ifstoklet, stoklet, ifstrslet,
 *          strslet, strsvec, ifppreg, pot, pre, grad, ntarg, targ,
 *          ifppregtarg, pottarg, pretarg, gradtarg, ier)
 *
 * stfmm2d is a real-valued Stokes FMM. All inputs and outputs are
 * REAL doubles (no complex). Internally it builds complex biharmonic
 * arrays and calls bhfmm2d.
 *
 * Buffer shapes (column-major doubles):
 *   source  : (2, ns)
 *   stoklet : (nd, 2, ns)         ignored if ifstoklet == 0
 *   strslet : (nd, 2, ns)         ignored if ifstrslet == 0
 *   strsvec : (nd, 2, ns)         ignored if ifstrslet == 0
 *   targ    : (2, nt)             ignored if nt == 0
 *   pot     : (nd, 2, ns)         written if ifppreg >= 1  (velocity)
 *   pre     : (nd, ns)            written if ifppreg >= 2  (pressure)
 *   grad    : (nd, 2, 2, ns)      written if ifppreg >= 3  (velocity grad)
 *   pottarg : (nd, 2, nt)         written if ifppregtarg >= 1
 *   pretarg : (nd, nt)            written if ifppregtarg >= 2
 *   gradtarg: (nd, 2, 2, nt)      written if ifppregtarg >= 3
 */
extern void stfmm2d_(const int32_t *nd, const double *eps,
                     const int32_t *nsource, const double *source,
                     const int32_t *ifstoklet, const double *stoklet,
                     const int32_t *ifstrslet, const double *strslet,
                     const double *strsvec,
                     const int32_t *ifppreg, double *pot, double *pre,
                     double *grad,
                     const int32_t *ntarg, const double *targ,
                     const int32_t *ifppregtarg, double *pottarg,
                     double *pretarg, double *gradtarg, int32_t *ier);

EXPORT(stfmm2d_w)
int32_t stfmm2d_w(int32_t nd, double eps, int32_t ns, double *source,
                  int32_t ifstoklet, double *stoklet,
                  int32_t ifstrslet, double *strslet, double *strsvec,
                  int32_t ifppreg, double *pot, double *pre, double *grad,
                  int32_t nt, double *targ, int32_t ifppregtarg,
                  double *pottarg, double *pretarg, double *gradtarg)
{
    int32_t ier = 0;

    stfmm2d_(&nd, &eps, &ns, source, &ifstoklet, stoklet,
             &ifstrslet, strslet, strsvec,
             &ifppreg, pot, pre, grad,
             &nt, targ, &ifppregtarg, pottarg, pretarg, gradtarg, &ier);

    return ier;
}

/*
 * Heap helpers used by the WASM JS shim to allocate/free buffers in the
 * WASM linear memory. (Native FFI uses koffi-style direct buffer passing
 * and does not need these.)
 */
#ifdef __EMSCRIPTEN__
EXPORT(my_malloc)
void *my_malloc(int32_t size) { return malloc((size_t)size); }

EXPORT(my_free)
void my_free(void *ptr) { free(ptr); }
#endif
