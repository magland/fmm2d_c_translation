/*
 * cfmm2d.h - C translation of src/laplace/cfmm2d.f
 *
 * Top-level user-facing 2D Laplace (Cauchy form) FMM driver and the
 * five subroutines from the Fortran source. Translated 1:1 from the
 * Fortran reference (same control flow, same allocations, same
 * operation order). OpenMP, error/print logging, and timing
 * (cpu_time / omp_get_wtime / second) are stripped.
 *
 * Cross-file dependencies (lndiv2d, pts_tree_*, dreorder*, init_carray,
 * l2dterms, computemnlists, computelists, l2dmpzero, l2dformmp*,
 * l2dformta*, l2dmpeval*, l2dtaeval*, l2dmpmp, l2dmploc, l2dlocloc,
 * c2d_direct*) are called via bare Fortran symbol names so the diff
 * test isolates this file from its dependencies; same-file calls
 * (cfmm2d -> cfmm2dmain / l2dmpalloc, cfmm2dmain -> cfmm2dexpc_direct
 * / cfmm2dpart_direct) use FNAME() dispatch.
 */

#ifndef FMM2D_CFMM2D_H
#define FMM2D_CFMM2D_H

#include "fmm2d_c.h"

void FNAME(cfmm2d)(const fint *nd, const double *eps,
                   const fint *ns, const double *sources,
                   const fint *ifcharge, const fcomplex *charge,
                   const fint *ifdipole, const fcomplex *dipstr,
                   fint *iper, const fint *ifpgh,
                   fcomplex *pot, fcomplex *grad, fcomplex *hess,
                   const fint *nt, const double *targ,
                   const fint *ifpghtarg,
                   fcomplex *pottarg, fcomplex *gradtarg,
                   fcomplex *hesstarg, fint *ier);

void FNAME(cfmm2dmain)(const fint *nd, const double *eps,
                       const fint *nsource, const double *sourcesort,
                       const fint *ifcharge, const fcomplex *chargesort,
                       const fint *ifdipole, const fcomplex *dipstrsort,
                       const fint *ntarget, const double *targetsort,
                       const fint *nexpc, const double *expcsort,
                       const fint *iaddr, double *rmlexp,
                       fcomplex *mptemp, const fint *lmptmp,
                       const fint *itree, const fint *ltree,
                       const fint *iptr, const fint *ndiv,
                       const fint *nlevels, const fint *nboxes,
                       const fint *iper, const double *boxsize,
                       const double *rscales, const double *centers,
                       const fint *laddr,
                       const fint *isrcse, const fint *itargse,
                       const fint *iexpcse, const fint *nterms,
                       const fint *ntj,
                       const fint *ifpgh, fcomplex *pot,
                       fcomplex *grad, fcomplex *hess,
                       const fint *ifpghtarg, fcomplex *pottarg,
                       fcomplex *gradtarg, fcomplex *hesstarg,
                       fcomplex *jsort, double *scjsort,
                       const fint *ifnear, double *timeinfo, fint *ier);

void FNAME(cfmm2dexpc_direct)(const fint *nd,
                              const fint *istart, const fint *iend,
                              const fint *jstart, const fint *jend,
                              const double *rscales, const fint *nlevels,
                              const double *source,
                              const fint *ifcharge, const fcomplex *charge,
                              const fint *ifdipole, const fcomplex *dipstr,
                              const double *targ, fcomplex *jexps,
                              const double *scj, const fint *ntj);

void FNAME(cfmm2dpart_direct)(const fint *nd,
                              const fint *istart, const fint *iend,
                              const fint *jstart, const fint *jend,
                              const double *source,
                              const fint *ifcharge, const fcomplex *charge,
                              const fint *ifdipole, const fcomplex *dipstr,
                              const double *targ, const fint *ifpgh,
                              fcomplex *pot, fcomplex *grad, fcomplex *hess,
                              const double *thresh);

void FNAME(l2dmpalloc)(const fint *nd, const fint *laddr, fint *iaddr,
                       const fint *nlevels, fint *lmptot,
                       const fint *nterms);

#endif /* FMM2D_CFMM2D_H */
