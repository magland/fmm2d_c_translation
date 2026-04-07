/*
 * l2dterms.c - C translation of subroutine l2dterms in
 * src/laplace/l2dterms.f
 *
 * Original Fortran:
 *
 *   subroutine l2dterms(eps, nterms, ier)
 *
 * Determines the number of terms in a Laplace 2D multipole expansion
 * needed to reach the requested precision eps. The method examines the
 * decay of rho^n / r^(n+1) for worst-case source rho = sqrt(2)/2 and
 * worst-case target r = 1.5.
 *
 * The Fortran routine declares z1, z2, jfun, hfun as complex *16 even
 * though they only ever take real values; we mirror that with fcomplex
 * to keep the translation byte-for-byte faithful (cpow / cabs operate
 * on the complex values, with zero imaginary parts).
 *
 * Other routines in the original file (l2dterms_far, l2dterms_list2,
 * l2dterms_list2w, l2dterms_list2e, l2dterms_list2ew, l2dterms_eval)
 * are not reachable from rfmm2d_ndiv and are intentionally not ported.
 */

#include "l2dterms.h"

void FNAME(l2dterms)(const double *eps, fint *nterms, fint *ier)
{
    /* complex *16 z1, z2, jfun(0:200), hfun(0:200) */
    fcomplex z1, z2;
    fcomplex jfun[201];
    fcomplex hfun[201];
    double xtemp1;
    fint i, j;
    fint ntmax;

    *ier = 0;

    ntmax = 100;

    /* z1 = 1.5d0 */
    z1 = 1.5 + 0.0 * I;
    /* do i = 0,ntmax ; hfun(i) = 1.0d0/(z1**(i+1)) ; enddo */
    for (i = 0; i <= ntmax; i++) {
        hfun[i] = 1.0 / cpow(z1, (double)(i + 1));
    }
    /* ccc      call prin2(' hfun is *',hfun,2*ntmax+2) */

    /* z2 = dsqrt(2d0)/2.d0 */
    z2 = (sqrt(2.0) / 2.0) + 0.0 * I;
    /* do i = 0,ntmax ; jfun(i) = z2**i ; enddo */
    for (i = 0; i <= ntmax; i++) {
        jfun[i] = cpow(z2, (double)i);
    }

    /* xtemp1 = cdabs(jfun(0)*hfun(0))  -- computed but immediately
       overwritten in the loop without being used; preserved as in the
       original Fortran. */
    xtemp1 = cabs(jfun[0] * hfun[0]);
    *nterms = 1;
    for (j = 2; j <= ntmax; j++) {
        xtemp1 = cabs(jfun[j] * hfun[j]);
        if (xtemp1 < *eps) {
            *nterms = j;
            return;
        }
    }
    return;
}
