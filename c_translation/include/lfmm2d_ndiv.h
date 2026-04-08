/*
 * lfmm2d_ndiv.h - C translation of src/laplace/lfmm2d_ndiv.f
 *
 * Top-level user-facing 2D Laplace FMM driver with COMPLEX-valued
 * charges and dipole strengths (and a real dipvec), plus externally
 * supplied subdivision criterion (ndiv, idivflag) and pass-through
 * ifnear, timeinfo, ier.
 *
 * This routine is a thin wrapper around cfmm2d_ndiv with nd doubled:
 * it splits each complex charge into a (re, im) pair, packs each pair
 * as two complex *16 entries (with zero imaginary part for the charge
 * pair, and ztmp = -(dipvec(:,1) + I*dipvec(:,2)) coupling for the
 * dipstr pair), calls cfmm2d_ndiv with nd2 = 2*nd, and recombines the
 * (2,nd,*) complex outputs into the original complex (nd,*),
 * complex (nd,2,*), complex (nd,3,*) shapes.
 *
 * The gradient/hessian unpacking follows the d/dz, d^2/dz^2 convention
 * of the Cauchy FMM (same as rfmm2d_ndiv):
 *   grad(:,1,:) =  Re(grad1)+i*Re(grad1_im)
 *   grad(:,2,:) = -Im(grad1)-i*Im(grad1_im)
 *   hess(:,1,:) =  Re(hess1)+i*Re(hess1_im)
 *   hess(:,2,:) = -Im(hess1)-i*Im(hess1_im)
 *   hess(:,3,:) = -hess(:,1,:)
 *
 * Translated 1:1 from the Fortran reference (same control flow, same
 * allocations, same operation order). OpenMP, error/print logging, and
 * timing are stripped.
 *
 * Cross-file dependency (cfmm2d_ndiv) is called via the bare Fortran
 * symbol name so the diff test isolates this file.
 */

#ifndef FMM2D_LFMM2D_NDIV_H
#define FMM2D_LFMM2D_NDIV_H

#include "fmm2d_c.h"

void FNAME(lfmm2d_ndiv)(const fint *nd, const double *eps,
                        const fint *ns, const double *sources,
                        const fint *ifcharge, const fcomplex *charge,
                        const fint *ifdipole, const fcomplex *dipstr,
                        const double *dipvec,
                        fint *iper, const fint *ifpgh,
                        fcomplex *pot, fcomplex *grad, fcomplex *hess,
                        const fint *nt, const double *targ,
                        const fint *ifpghtarg,
                        fcomplex *pottarg, fcomplex *gradtarg,
                        fcomplex *hesstarg,
                        const fint *ndiv, const fint *idivflag,
                        const fint *ifnear, double *timeinfo, fint *ier);

#endif /* FMM2D_LFMM2D_NDIV_H */
