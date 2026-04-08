/*
 * bhfmm2d.h - C translation of src/biharmonic/bhfmm2d.f
 *
 * Top-level user-facing 2D biharmonic (Stokes) FMM driver and the
 * supporting subroutines:
 *   bhfmm2d            - top-level driver (allocates tree, calls main)
 *   bhfmm2dmain        - main 8-step FMM engine
 *   bhfmm2dpart_direct - direct particle-to-particle dispatcher
 *   bh2dmpalloc        - lay out workspace for mpole/local expansions
 *
 * Translated 1:1 from the Fortran reference: same control flow, same
 * allocations, same operation order. OpenMP, logging, and timing are
 * stripped. The biharmonic routines support ifpgh up to 2 only
 * (hessian is not implemented in the Fortran source either).
 */

#ifndef FMM2D_BHFMM2D_H
#define FMM2D_BHFMM2D_H

#include "fmm2d_c.h"

void FNAME(bhfmm2d)(const fint *nd, const double *eps,
                    const fint *ns, const double *sources,
                    const fint *ifcharge, const fcomplex *charge,
                    const fint *ifdipole, const fcomplex *dip,
                    fint *iper, const fint *ifpgh,
                    fcomplex *pot, fcomplex *grad, fcomplex *hess,
                    const fint *nt, const double *targ,
                    const fint *ifpghtarg,
                    fcomplex *pottarg, fcomplex *gradtarg,
                    fcomplex *hesstarg, fint *ier);

void FNAME(bhfmm2dmain)(const fint *nd, const double *eps,
                        const fint *nsource, const double *sourcesort,
                        const fint *ifcharge, const fcomplex *chargesort,
                        const fint *ifdipole, const fcomplex *dipsort,
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
                        fcomplex *jsort, double *scjsort);

void FNAME(bhfmm2dpart_direct)(const fint *nd,
                               const fint *istart, const fint *iend,
                               const fint *jstart, const fint *jend,
                               const double *source,
                               const fint *ifcharge, const fcomplex *charge,
                               const fint *ifdipole, const fcomplex *dip,
                               const double *targ, const fint *ifpgh,
                               fcomplex *pot, fcomplex *grad, fcomplex *hess,
                               const double *thresh);

void FNAME(bh2dmpalloc)(const fint *nd, const fint *laddr, fint *iaddr,
                        const fint *nlevels, fint *lmptot,
                        const fint *nterms);

#endif /* FMM2D_BHFMM2D_H */
