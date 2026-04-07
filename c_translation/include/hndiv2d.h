/*
 * hndiv2d.h - C translation of src/helmholtz/hndiv2d.f
 */

#ifndef FMM2D_HNDIV2D_H
#define FMM2D_HNDIV2D_H

#include "fmm2d_c.h"

/*
 * void hndiv2d(eps, ns, nt, ifcharge, ifdipole, ifpgh, ifpghtarg,
 *              ndiv, idivflag)
 *
 * Sets ndiv (subdivision criterion, points-per-box) and idivflag based
 * on the requested precision eps. The output ndiv depends only on eps;
 * idivflag is always 0 in this routine. ns,nt are only used as the
 * fallback when eps is below 0.5e-15.
 */
void FNAME(hndiv2d)(
    const double *eps,
    const fint *ns, const fint *nt,
    const fint *ifcharge, const fint *ifdipole,
    const fint *ifpgh, const fint *ifpghtarg,
    fint *ndiv, fint *idivflag);

#endif
