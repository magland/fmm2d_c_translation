/*
 * cfmm2d_ndiv.h - C translation of src/laplace/cfmm2d_ndiv.f
 *
 * Top-level user-facing 2D Laplace (Cauchy form) FMM driver with
 * externally supplied subdivision criterion (ndiv, idivflag) and
 * pass-through ifnear, timeinfo, ier. Translated 1:1 from the Fortran
 * reference (same control flow, same allocations, same operation
 * order). OpenMP, error/print logging, and timing (cpu_time /
 * omp_get_wtime / second) are stripped.
 *
 * Cross-file dependencies (pts_tree_*, dreorder*, l2dterms, cfmm2dmain,
 * l2dmpalloc) are called via bare Fortran symbol names so the diff
 * test isolates this file from its dependencies.
 */

#ifndef FMM2D_CFMM2D_NDIV_H
#define FMM2D_CFMM2D_NDIV_H

#include "fmm2d_c.h"

void FNAME(cfmm2d_ndiv)(const fint *nd, const double *eps,
                        const fint *ns, const double *sources,
                        const fint *ifcharge, const fcomplex *charge,
                        const fint *ifdipole, const fcomplex *dipstr,
                        fint *iper, const fint *ifpgh,
                        fcomplex *pot, fcomplex *grad, fcomplex *hess,
                        const fint *nt, const double *targ,
                        const fint *ifpghtarg,
                        fcomplex *pottarg, fcomplex *gradtarg,
                        fcomplex *hesstarg,
                        const fint *ndiv, const fint *idivflag,
                        const fint *ifnear, double *timeinfo, fint *ier);

#endif /* FMM2D_CFMM2D_NDIV_H */
