/*
 * l2dterms.h - C translation of subroutine l2dterms in
 * src/laplace/l2dterms.f
 */

#ifndef FMM2D_L2DTERMS_H
#define FMM2D_L2DTERMS_H

#include "fmm2d_c.h"

/*
 * void l2dterms(eps, nterms, ier)
 *
 * Determines the number of terms `nterms` needed in a Laplace 2D
 * multipole expansion to reach the requested precision `eps`. The
 * estimate is based on examining the decay of rho^n / r^(n+1) for a
 * worst-case source rho = sqrt(2)/2 and worst-case target r = 1.5.
 * `ier` is always set to 0.
 */
void FNAME(l2dterms)(const double *eps, fint *nterms, fint *ier);

#endif /* FMM2D_L2DTERMS_H */
