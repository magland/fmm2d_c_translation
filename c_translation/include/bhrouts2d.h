/*
 * bhrouts2d.h - C translation of src/biharmonic/bhrouts2d.f
 *
 * Multipole and local expansion form / eval / translate routines for
 * the 2D biharmonic (Stokes) FMM. Each expansion has FIVE sets of
 * coefficients per density slot, so the coefficient buffer has shape
 *   mpole(nd, 5, 0:nterms)
 * in Fortran column-major storage. All routines INCREMENT their
 * output buffers; they never overwrite them. The translation
 * preserves the Fortran left-to-right floating-point operation order
 * exactly so the C and Fortran objects can be compared bit-for-bit at
 * -O0.
 */

#ifndef FMM2D_BHROUTS2D_H
#define FMM2D_BHROUTS2D_H

#include "fmm2d_c.h"

void FNAME(bh2dmpevalp)(const fint *nd, const double *rscale,
                        const double *center, const fcomplex *mpole,
                        const fint *nterms, const double *ztarg,
                        const fint *ntarg, fcomplex *vel);

void FNAME(bh2dmpevalg)(const fint *nd, const double *rscale,
                        const double *center, const fcomplex *mpole,
                        const fint *nterms, const double *ztarg,
                        const fint *ntarg, fcomplex *vel, fcomplex *grad);

void FNAME(bh2dtaevalp)(const fint *nd, const double *rscale,
                        const double *center, const fcomplex *mpole,
                        const fint *nterms, const double *ztarg,
                        const fint *ntarg, fcomplex *vel);

void FNAME(bh2dtaevalg)(const fint *nd, const double *rscale,
                        const double *center, const fcomplex *mpole,
                        const fint *nterms, const double *ztarg,
                        const fint *ntarg, fcomplex *vel, fcomplex *grad);

void FNAME(bh2dformmpd)(const fint *nd, const double *rscale,
                        const double *sources, const fint *ns,
                        const fcomplex *dip, const double *center,
                        const fint *nterms, fcomplex *mpole);

void FNAME(bh2dformmpc)(const fint *nd, const double *rscale,
                        const double *sources, const fint *ns,
                        const fcomplex *c1, const double *center,
                        const fint *nterms, fcomplex *mpole);

void FNAME(bh2dformmpcd)(const fint *nd, const double *rscale,
                         const double *sources, const fint *ns,
                         const fcomplex *c1, const fcomplex *dip,
                         const double *center, const fint *nterms,
                         fcomplex *mpole);

void FNAME(bh2dformtac)(const fint *nd, const double *rscale,
                        const double *sources, const fint *ns,
                        const fcomplex *c1, const double *center,
                        const fint *nterms, fcomplex *mpole);

void FNAME(bh2dformtad)(const fint *nd, const double *rscale,
                        const double *sources, const fint *ns,
                        const fcomplex *dip, const double *center,
                        const fint *nterms, fcomplex *mpole);

void FNAME(bh2dformtacd)(const fint *nd, const double *rscale,
                         const double *sources, const fint *ns,
                         const fcomplex *c1, const fcomplex *dip,
                         const double *center, const fint *nterms,
                         fcomplex *mpole);

void FNAME(bh2dlocloc)(const fint *nd,
                       const double *rscale1, const double *c1,
                       const fcomplex *hexp, const fint *nterms1,
                       const double *rscale2, const double *c2,
                       fcomplex *jexp, const fint *nterms2,
                       const double *carray, const fint *ldc);

void FNAME(bh2dmpmp)(const fint *nd,
                     const double *rscale1, const double *c1,
                     const fcomplex *hexp, const fint *nterms1,
                     const double *rscale2, const double *c2,
                     fcomplex *jexp, const fint *nterms2,
                     const double *carray, const fint *ldc);

void FNAME(bh2dmploc)(const fint *nd,
                      const double *rscale1, const double *c1,
                      const fcomplex *hexp, const fint *nterms1,
                      const double *rscale2, const double *c2,
                      fcomplex *jexp, const fint *nterms2,
                      const double *carray, const fint *ldc);

void FNAME(bh2dmpzero)(const fint *nd, fcomplex *mpole, const fint *nterms);

#endif /* FMM2D_BHROUTS2D_H */
