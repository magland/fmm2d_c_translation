/*
 * bh2dterms.c - C translation of subroutine bh2dterms in
 * src/biharmonic/bh2dterms.f
 *
 * Original Fortran:
 *
 *   subroutine bh2dterms(eps, nterms, ier)
 *
 * Determines the number of terms in a biharmonic 2D multipole
 * expansion needed to reach the requested precision eps. The method
 * examines the decay of rho^n / r^(n+1) for worst-case source
 * rho = sqrt(2)/2 and worst-case target r = 1.5.
 *
 * The Fortran routine declares zk, z1, z2, z3, ztmp, ht0, ht1, ht2,
 * jfun, hfun as complex *16, but only z1, z2, jfun, hfun are actually
 * used. We preserve the unused locals as (void) casts so the
 * translation stays 1:1 with the Fortran source without tripping
 * -Wunused warnings.
 *
 * The Fortran expression `z1**(i+1)` with complex base and integer
 * exponent is evaluated by gfortran as a left-to-right running product
 * of i+1 multiplications. We match that exactly with an explicit
 * running product rather than cpow, so the floating-point results are
 * bit-identical to the Fortran reference.
 *
 * Other routines in the original file (bh2dterms_far, bh2dterms_list2)
 * are not reachable from bhfmm2d and are intentionally not ported.
 */

#include "bh2dterms.h"

void FNAME(bh2dterms)(const double *eps, fint *nterms, fint *ier)
{
    /* complex *16 zk, z1, z2, z3, jfun(0:2000), ht0,
     *             ht1, ht2, ztmp,
     *             hfun(0:2000) */
    fcomplex zk, z1, z2, z3, ztmp, ht0, ht1, ht2;
    fcomplex jfun[2001];
    fcomplex hfun[2001];
    fcomplex cum;
    double xtemp1;
    fint i, j;
    fint ntmax;

    /* zk, z3, ztmp, ht0, ht1, ht2 are declared in the Fortran source
     * but never read from or written to. Silence unused-variable
     * warnings while preserving the declarations. */
    zk   = 0.0 + 0.0 * I;
    z3   = 0.0 + 0.0 * I;
    ztmp = 0.0 + 0.0 * I;
    ht0  = 0.0 + 0.0 * I;
    ht1  = 0.0 + 0.0 * I;
    ht2  = 0.0 + 0.0 * I;
    (void)zk;
    (void)z3;
    (void)ztmp;
    (void)ht0;
    (void)ht1;
    (void)ht2;

    *ier = 0;

    ntmax = 1000;

    /* z1 = 1.5d0 */
    z1 = 1.5 + 0.0 * I;
    /* do i = 0,ntmax ; hfun(i) = 1.0d0/(z1**(i+1)) ; enddo
     *
     * Fortran evaluates z1**(i+1) as a running product of i+1 copies
     * of z1 multiplied left-to-right. We mirror that exactly with a
     * cumulative product seeded at 1: after iteration i, cum holds
     * z1**(i+1). */
    cum = 1.0 + 0.0 * I;
    for (i = 0; i <= ntmax; i++) {
        cum = cum * z1;
        hfun[i] = 1.0 / cum;
    }
    /* ccc      call prin2(' hfun is *',hfun,2*ntmax+2) */

    /* z2 = dsqrt(2d0)/2.d0 */
    z2 = (sqrt(2.0) / 2.0) + 0.0 * I;
    /* do i = 0,ntmax ; jfun(i) = z2**i ; enddo
     *
     * Running product again: after iteration i, cum = z2**i. The i=0
     * case is the empty product, which is 1 (matching Fortran
     * z2**0 = 1). */
    cum = 1.0 + 0.0 * I;
    for (i = 0; i <= ntmax; i++) {
        jfun[i] = cum;
        cum = cum * z2;
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
