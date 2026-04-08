/*
 * stfmm2d.h - C translation of src/stokes/stfmm2d.f
 *
 * 2D Stokes FMM driver. Thin wrapper around bhfmm2d:
 *   1. Allocate complex *16 (nd, 2, nsource) charge and (nd, 3, nsource)
 *      dip work buffers.
 *   2. Encode the Stokeslet and stresslet vectors into the biharmonic
 *      complex charge / dip arrays.
 *   3. Call bhfmm2d.
 *   4. Decode the biharmonic potl/gradl outputs back into the Stokes
 *      pot (velocity), pre (pressure), and grad (velocity gradient).
 *
 * Translated 1:1 from the Fortran reference. Currently only ifstrslet
 * == 1 (type-I stresslet) is supported, matching the Fortran source's
 * "TODO" status.
 */

#ifndef FMM2D_STFMM2D_H
#define FMM2D_STFMM2D_H

#include "fmm2d_c.h"

void FNAME(stfmm2d)(const fint *nd, const double *eps,
                    const fint *nsource, const double *source,
                    const fint *ifstoklet, const double *stoklet,
                    const fint *ifstrslet, const double *strslet,
                    const double *strsvec,
                    const fint *ifppreg, double *pot, double *pre,
                    double *grad,
                    const fint *ntarg, const double *targ,
                    const fint *ifppregtarg, double *pottarg,
                    double *pretarg, double *gradtarg, fint *ier);

#endif /* FMM2D_STFMM2D_H */
