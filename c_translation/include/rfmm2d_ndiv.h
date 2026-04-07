/*
 * rfmm2d_ndiv.h - C translation of src/laplace/rfmm2d_ndiv.f
 *
 * Top-level user-facing 2D Laplace FMM driver with REAL-valued charges
 * and dipole strengths, plus externally supplied subdivision criterion
 * (ndiv, idivflag) and pass-through ifnear, timeinfo, ier.
 *
 * This routine is a thin wrapper around cfmm2d_ndiv: it allocates
 * complex *16 work buffers, converts the real inputs (charge, dipstr,
 * dipvec) to complex (charge1, dipstr1 = -dipstr * (dipvec_x + I*dipvec_y)),
 * calls cfmm2d_ndiv, and unpacks the complex outputs back to real.
 *
 * The gradient/hessian unpacking follows the d/dz, d^2/dz^2 convention
 * of the Cauchy FMM:
 *   grad(:,1,:) =  Re(grad1),  grad(:,2,:) = -Im(grad1)
 *   hess(:,1,:) =  Re(hess1),  hess(:,2,:) = -Im(hess1),
 *   hess(:,3,:) = -hess(:,1,:)
 *
 * Translated 1:1 from the Fortran reference (same control flow, same
 * allocations, same operation order). OpenMP, error/print logging, and
 * timing are stripped.
 *
 * Cross-file dependency (cfmm2d_ndiv) is called via the bare Fortran
 * symbol name so the diff test isolates this file.
 */

#ifndef FMM2D_RFMM2D_NDIV_H
#define FMM2D_RFMM2D_NDIV_H

#include "fmm2d_c.h"

void FNAME(rfmm2d_ndiv)(const fint *nd, const double *eps,
                        const fint *ns, const double *sources,
                        const fint *ifcharge, const double *charge,
                        const fint *ifdipole, const double *dipstr,
                        const double *dipvec,
                        fint *iper, const fint *ifpgh,
                        double *pot, double *grad, double *hess,
                        const fint *nt, const double *targ,
                        const fint *ifpghtarg,
                        double *pottarg, double *gradtarg,
                        double *hesstarg,
                        const fint *ndiv, const fint *idivflag,
                        const fint *ifnear, double *timeinfo, fint *ier);

#endif /* FMM2D_RFMM2D_NDIV_H */
