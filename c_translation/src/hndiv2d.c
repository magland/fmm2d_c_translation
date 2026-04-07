/*
 * hndiv2d.c - C translation of src/helmholtz/hndiv2d.f
 *
 * Original Fortran:
 *
 *   subroutine hndiv2d(eps,ns,nt,ifcharge,ifdipole,ifpgh,
 *                      ifpghtarg,ndiv,idivflag)
 *
 * Estimates a subdivision criterion (ndiv = points-per-box) from the
 * requested precision eps. idivflag is always 0.
 */

#include "hndiv2d.h"

void FNAME(hndiv2d)(
    const double *eps,
    const fint *ns, const fint *nt,
    const fint *ifcharge, const fint *ifdipole,
    const fint *ifpgh, const fint *ifpghtarg,
    fint *ndiv, fint *idivflag)
{
    (void)ifcharge;
    (void)ifdipole;
    (void)ifpgh;
    (void)ifpghtarg;

    *idivflag = 0;

    if (*eps >= 0.5e-0) {
        *ndiv = 3;
    } else if (*eps >= 0.5e-1) {
        *ndiv = 5;
    } else if (*eps >= 0.5e-2) {
        *ndiv = 8;
    } else if (*eps >= 0.5e-3) {
        *ndiv = 10;
    } else if (*eps >= 0.5e-6) {
        *ndiv = 15;
    } else if (*eps >= 0.5e-9) {
        *ndiv = 20;
    } else if (*eps >= 0.5e-12) {
        *ndiv = 25;
    } else if (*eps >= 0.5e-15) {
        *ndiv = 45;
    } else {
        *ndiv = *ns + *nt;
    }
}
