/*
 * laprouts2d.h - C translation of src/laplace/laprouts2d.f
 *
 * Multipole and local expansion form/eval/translate routines for the
 * 2D Laplace FMM (Cauchy form). The expansion convention follows the
 * Fortran library: log(|z|) for the n=0 multipole/local term, and
 * complex powers (z/rscale)^n / (rscale/z)^n for the higher order
 * terms. All routines INCREMENT their output buffers; they do not
 * overwrite them. The same is true for the M2M / L2L / M2L
 * translation routines.
 */

#ifndef FMM2D_LAPROUTS2D_H
#define FMM2D_LAPROUTS2D_H

#include "fmm2d_c.h"

void FNAME(l2dformmpc)(const fint *nd, const double *rscale, const double *source,
                       const fint *ns, const fcomplex *charge, const double *center,
                       const fint *nterms, fcomplex *mpole);

void FNAME(l2dformmpd)(const fint *nd, const double *rscale, const double *source,
                       const fint *ns, const fcomplex *dipstr, const double *center,
                       const fint *nterms, fcomplex *mpole);

void FNAME(l2dformmpcd)(const fint *nd, const double *rscale, const double *source,
                        const fint *ns, const fcomplex *charge, const fcomplex *dipstr,
                        const double *center, const fint *nterms, fcomplex *mpole);

void FNAME(l2dmpevalp)(const fint *nd, const double *rscale, const double *center,
                       const fcomplex *mpole, const fint *nterms,
                       const double *ztarg, const fint *ntarg, fcomplex *pot);

void FNAME(l2dmpevalg)(const fint *nd, const double *rscale, const double *center,
                       const fcomplex *mpole, const fint *nterms,
                       const double *ztarg, const fint *ntarg,
                       fcomplex *pot, fcomplex *grad);

void FNAME(l2dmpevalh)(const fint *nd, const double *rscale, const double *center,
                       const fcomplex *mpole, const fint *nterms,
                       const double *ztarg, const fint *ntarg,
                       fcomplex *pot, fcomplex *grad, fcomplex *hess);

void FNAME(l2dformtac)(const fint *nd, const double *rscale, const double *source,
                       const fint *ns, const fcomplex *charge, const double *center,
                       const fint *nterms, fcomplex *local);

void FNAME(l2dformtad)(const fint *nd, const double *rscale, const double *source,
                       const fint *ns, const fcomplex *dipstr, const double *center,
                       const fint *nterms, fcomplex *local);

void FNAME(l2dformtacd)(const fint *nd, const double *rscale, const double *source,
                        const fint *ns, const fcomplex *charge, const fcomplex *dipstr,
                        const double *center, const fint *nterms, fcomplex *local);

void FNAME(l2dtaevalp)(const fint *nd, const double *rscale, const double *center,
                       const fcomplex *local, const fint *nterms,
                       const double *ztarg, const fint *ntarg, fcomplex *pot);

void FNAME(l2dtaevalg)(const fint *nd, const double *rscale, const double *center,
                       const fcomplex *local, const fint *nterms,
                       const double *ztarg, const fint *ntarg,
                       fcomplex *pot, fcomplex *grad);

void FNAME(l2dtaevalh)(const fint *nd, const double *rscale, const double *center,
                       const fcomplex *local, const fint *nterms,
                       const double *ztarg, const fint *ntarg,
                       fcomplex *pot, fcomplex *grad, fcomplex *hess);

void FNAME(l2dmpmp)(const fint *nd,
                    const double *rscale1, const double *center1,
                    const fcomplex *hexp1, const fint *nterms1,
                    const double *rscale2, const double *center2,
                    fcomplex *hexp2, const fint *nterms2,
                    const double *carray, const fint *ldc);

void FNAME(l2dlocloc)(const fint *nd,
                      const double *rscale1, const double *center1,
                      const fcomplex *jexp1, const fint *nterms1,
                      const double *rscale2, const double *center2,
                      fcomplex *jexp2, const fint *nterms2,
                      const double *carray, const fint *ldc);

void FNAME(l2dmploc)(const fint *nd,
                     const double *rscale1, const double *center1,
                     const fcomplex *hexp1, const fint *nterms1,
                     const double *rscale2, const double *center2,
                     fcomplex *jexp2, const fint *nterms2,
                     const double *carray, const fint *ldc);

void FNAME(l2dmpzero)(const fint *nd, fcomplex *mpole, const fint *nterms);

#endif /* FMM2D_LAPROUTS2D_H */
