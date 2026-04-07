/*
 * rfmm2d_wrapper.c
 *
 * C wrapper around the Fortran-ABI rfmm2d entry point, used for both the
 * WebAssembly build and the native shared-library build that numbl loads.
 *
 * Exports a single function rfmm2d_w(...) that:
 *   1. Calls hndiv2d_ to pick a default subdivision criterion (ndiv).
 *   2. Calls rfmm2d_ndiv_ to actually run the FMM.
 * Both callees use the gfortran ABI (lowercase + trailing underscore,
 * arguments passed by pointer, fortran integers are int32).
 *
 * Linkage:
 *   WASM    - linked against c_translation/src .c files compiled with
 *             -DFMM2D_DROP_IN, which exports the bare Fortran names.
 *   native  - linked against the top-level libfmm2d.a built by
 *             gfortran, which exports the real Fortran symbols.
 * In both cases the externs declared here resolve to the same names.
 *
 * Only the rfmm2d entry point is currently supported.
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
