/*
 * fmmcommon2d.c - C translation of selected routines from
 * src/common/fmmcommon2d.f
 *
 * Translated routines:
 *
 *   subroutine dreorderf(ndim, n, arr, arrsort, iarr)
 *   subroutine dreorderi(ndim, n, arr, arrsort, iarr)
 *   subroutine init_carray(carray, ldc)
 *
 * The other routines in the original Fortran file (ireorderf,
 * geterrstr) are intentionally not ported because they are not
 * reachable from rfmm2d_ndiv.
 *
 * OpenMP pragmas (C$OMP) in the Fortran source are comments and are
 * ignored by this translation; the resulting code runs sequentially
 * but produces identical numerical output.
 */

#include "fmmcommon2d.h"

/*
 * dreorderf: gather arr columns into arrsort according to iarr.
 *
 * Fortran:
 *   double precision arr(ndim,n), arrsort(ndim,n)
 *   integer iarr(n)
 *   do i = 1,n
 *      do idim = 1,ndim
 *         arrsort(idim,i) = arr(idim, iarr(i))
 *      enddo
 *   enddo
 */
void FNAME(dreorderf)(const fint *ndim, const fint *n,
                      const double *arr, double *arrsort,
                      const fint *iarr)
{
    fint i, idim;
    fint ndim_v = *ndim;
    fint n_v = *n;

    for (i = 1; i <= n_v; i++) {
        for (idim = 1; idim <= ndim_v; idim++) {
            arrsort[FA2(idim, i, ndim_v)] =
                arr[FA2(idim, iarr[i - 1], ndim_v)];
        }
    }

    return;
}

/*
 * dreorderi: scatter arr columns into arrsort according to iarr
 * (inverse of dreorderf).
 *
 * Fortran:
 *   double precision arr(ndim,1), arrsort(ndim,1)
 *   integer iarr(1)
 *   do i = 1,n
 *      do idim = 1,ndim
 *         arrsort(idim, iarr(i)) = arr(idim, i)
 *      enddo
 *   enddo
 */
void FNAME(dreorderi)(const fint *ndim, const fint *n,
                      const double *arr, double *arrsort,
                      const fint *iarr)
{
    fint i, idim;
    fint ndim_v = *ndim;
    fint n_v = *n;

    for (i = 1; i <= n_v; i++) {
        for (idim = 1; idim <= ndim_v; idim++) {
            arrsort[FA2(idim, iarr[i - 1], ndim_v)] =
                arr[FA2(idim, i, ndim_v)];
        }
    }

    return;
}

/*
 * init_carray: build the binomial coefficient table.
 *
 * The Fortran source declares
 *
 *     real *8 carray(0:ldc, 0:ldc)
 *
 * i.e. both subscripts are 0-based, and the leading dimension is
 * (ldc+1). The Fortran access carray(l, m) therefore corresponds to
 * the C offset
 *
 *     carray[m * (ldc+1) + l]
 *
 * with 0 <= l, m <= ldc. Note that we do NOT use the FA2 macro here
 * because FA2 assumes 1-based indexing; the off-by-one differs.
 *
 * Fortran:
 *   do l = 0, ldc
 *      carray(l, 0) = 1.0d0
 *   enddo
 *   do m = 1, ldc
 *      carray(m, m) = 1.0d0
 *      do l = m+1, ldc
 *         carray(l, m) = carray(l-1, m) + carray(l-1, m-1)
 *      enddo
 *   enddo
 */
void FNAME(init_carray)(double *carray, const fint *ldc)
{
    fint l, m;
    fint ldc_v = *ldc;
    fint ld = ldc_v + 1; /* leading dimension of carray(0:ldc, 0:ldc) */

    for (l = 0; l <= ldc_v; l++) {
        carray[0 * ld + l] = 1.0;
    }
    for (m = 1; m <= ldc_v; m++) {
        carray[m * ld + m] = 1.0;
        for (l = m + 1; l <= ldc_v; l++) {
            carray[m * ld + l] =
                carray[m * ld + (l - 1)] +
                carray[(m - 1) * ld + (l - 1)];
        }
    }

    return;
}
