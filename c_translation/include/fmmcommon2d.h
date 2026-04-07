/*
 * fmmcommon2d.h - C translation of selected routines from
 * src/common/fmmcommon2d.f
 *
 * Only the routines reachable from rfmm2d_ndiv are translated:
 *   dreorderf, dreorderi, init_carray
 *
 * The other routines in the original Fortran file (ireorderf,
 * geterrstr) are intentionally not ported.
 */

#ifndef FMM2D_FMMCOMMON2D_H
#define FMM2D_FMMCOMMON2D_H

#include "fmm2d_c.h"

/*
 * void dreorderf(ndim, n, arr, arrsort, iarr)
 *
 * Gather: arrsort(:,i) = arr(:,iarr(i)), for i = 1..n.
 * arr and arrsort are dimensioned (ndim, n); iarr has length n.
 */
void FNAME(dreorderf)(const fint *ndim, const fint *n,
                      const double *arr, double *arrsort,
                      const fint *iarr);

/*
 * void dreorderi(ndim, n, arr, arrsort, iarr)
 *
 * Scatter (inverse of dreorderf): arrsort(:,iarr(i)) = arr(:,i),
 * for i = 1..n. arr and arrsort are dimensioned (ndim, n);
 * iarr has length n.
 */
void FNAME(dreorderi)(const fint *ndim, const fint *n,
                      const double *arr, double *arrsort,
                      const fint *iarr);

/*
 * void init_carray(carray, ldc)
 *
 * Fills a binomial coefficient table. The Fortran source declares
 *   real *8 carray(0:ldc, 0:ldc)
 * i.e. 0-based indexing in both dimensions, with leading dimension
 * (ldc+1). The caller must therefore allocate (ldc+1)*(ldc+1) doubles.
 */
void FNAME(init_carray)(double *carray, const fint *ldc);

#endif /* FMM2D_FMMCOMMON2D_H */
