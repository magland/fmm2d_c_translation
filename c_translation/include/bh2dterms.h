/*
 * bh2dterms.h - C translation of subroutine bh2dterms in
 * src/biharmonic/bh2dterms.f
 */

#ifndef FMM2D_BH2DTERMS_H
#define FMM2D_BH2DTERMS_H

#include "fmm2d_c.h"

/*
 * void bh2dterms(eps, nterms, ier)
 *
 * Determines the number of terms `nterms` needed in a biharmonic 2D
 * multipole expansion to reach the requested precision `eps`. The
 * estimate is based on examining the decay of rho^n / r^(n+1) for a
 * worst-case source rho = sqrt(2)/2 and worst-case target r = 1.5.
 * `ier` is always set to 0.
 *
 * The sibling routines bh2dterms_far and bh2dterms_list2 in the same
 * Fortran source file are not reachable from bhfmm2d and are
 * intentionally not ported (mirroring the l2dterms.c precedent).
 */
void FNAME(bh2dterms)(const double *eps, fint *nterms, fint *ier);

#endif /* FMM2D_BH2DTERMS_H */
