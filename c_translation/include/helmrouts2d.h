/*
 * helmrouts2d.h - C translation of src/helmholtz/helmrouts2d.f
 */

#ifndef FMM2D_HELMROUTS2D_H
#define FMM2D_HELMROUTS2D_H

#include "fmm2d_c.h"

void FNAME(h2dformmpc)(const fint *nd, const fcomplex *zk,
    const double *rscale, const double *source, const fint *ns,
    const fcomplex *charge, const double *center, const fint *nterms,
    fcomplex *mpole);

void FNAME(h2dformmpd)(const fint *nd, const fcomplex *zk,
    const double *rscale, const double *source, const fint *ns,
    const fcomplex *dipstr, const double *dipvec,
    const double *center, const fint *nterms, fcomplex *mpole);

void FNAME(h2dformmpcd)(const fint *nd, const fcomplex *zk,
    const double *rscale, const double *source, const fint *ns,
    const fcomplex *charge, const fcomplex *dipstr,
    const double *dipvec, const double *center, const fint *nterms,
    fcomplex *mpole);

void FNAME(h2dmpevalp)(const fint *nd, const fcomplex *zk,
    const double *rscale, const double *center, const fcomplex *mpole,
    const fint *nterms, const double *ztarg, const fint *ntarg,
    fcomplex *pot1);

void FNAME(h2dmpevalg)(const fint *nd, const fcomplex *zk,
    const double *rscale, const double *center, const fcomplex *mpole,
    const fint *nterms, const double *ztarg, const fint *ntarg,
    fcomplex *pot1, fcomplex *grad1);

void FNAME(h2dmpevalh)(const fint *nd, const fcomplex *zk,
    const double *rscale, const double *center, const fcomplex *mpole,
    const fint *nterms, const double *ztarg, const fint *ntarg,
    fcomplex *pot1, fcomplex *grad1, fcomplex *hess1);

void FNAME(h2dformtac)(const fint *nd, const fcomplex *zk,
    const double *rscale, const double *source, const fint *ns,
    const fcomplex *charge, const double *center, const fint *nterms,
    fcomplex *local);

void FNAME(h2dformtad)(const fint *nd, const fcomplex *zk,
    const double *rscale, const double *source, const fint *ns,
    const fcomplex *dipstr, const double *dipvec,
    const double *center, const fint *nterms, fcomplex *local);

void FNAME(h2dformtacd)(const fint *nd, const fcomplex *zk,
    const double *rscale, const double *source, const fint *ns,
    const fcomplex *charge, const fcomplex *dipstr,
    const double *dipvec, const double *center, const fint *nterms,
    fcomplex *local);

void FNAME(h2dtaevalp)(const fint *nd, const fcomplex *zk,
    const double *rscale, const double *center, const fcomplex *local,
    const fint *nterms, const double *ztarg, const fint *ntarg,
    fcomplex *pot1);

void FNAME(h2dtaevalg)(const fint *nd, const fcomplex *zk,
    const double *rscale, const double *center, const fcomplex *local,
    const fint *nterms, const double *ztarg, const fint *ntarg,
    fcomplex *pot1, fcomplex *grad1);

void FNAME(h2dtaevalh)(const fint *nd, const fcomplex *zk,
    const double *rscale, const double *center, const fcomplex *local,
    const fint *nterms, const double *ztarg, const fint *ntarg,
    fcomplex *pot1, fcomplex *grad1, fcomplex *hess1);

void FNAME(h2dmpmp)(const fint *nd, const fcomplex *zk,
    const double *rscale1, const double *center1, const fcomplex *hexp1,
    const fint *nterms1, const double *rscale2, const double *center2,
    fcomplex *hexp2, const fint *nterms2);

void FNAME(h2dlocloc)(const fint *nd, const fcomplex *zk,
    const double *rscale1, const double *center1, const fcomplex *jexp1,
    const fint *nterms1, const double *rscale2, const double *center2,
    fcomplex *jexp2, const fint *nterms2);

void FNAME(h2dmploc)(const fint *nd, const fcomplex *zk,
    const double *rscale1, const double *center1, const fcomplex *hexp,
    const fint *nterms1, const double *rscale2, const double *center2,
    fcomplex *jexp, const fint *nterms2);

#endif
