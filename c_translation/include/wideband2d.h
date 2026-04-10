/*
 * wideband2d.h - C translation of src/helmholtz/wideband2d.f
 */

#ifndef FMM2D_WIDEBAND2D_H
#define FMM2D_WIDEBAND2D_H

#include "fmm2d_c.h"

void FNAME(h2dloclochf)(const fint *nd, const fcomplex *zk,
    const double *rscale1, const double *center1, const fcomplex *sig,
    const fint *nterms1, const fint *nsig, const double *rscale2,
    const double *center2, fcomplex *hexp2, const fint *nterms2,
    const fcomplex *transvec, fcomplex *wsave);

void FNAME(h2dmpmphf)(const fint *nd, const fcomplex *zk,
    const double *rscale1, const double *center1, const fcomplex *hexp1,
    const fint *nterms1, const double *rscale2, const double *center2,
    fcomplex *sig2, const fint *nterms2, const fint *nsig,
    fcomplex *wsave, const fcomplex *transvec);

void FNAME(h2d_mptosig)(const fint *nd, const fint *nterms1,
    const fint *nsig, const fcomplex *hexp, fcomplex *sig,
    fcomplex *wsave);

void FNAME(h2d_mkmpshift)(const fcomplex *zk, const double *center1,
    const fint *nterms1, const double *center2, const fint *nterms2,
    const fint *nsig, fcomplex *wsave, fcomplex *transvec);

void FNAME(h2d_diagtrans)(const fint *nd, const fint *nsig,
    const fcomplex *sig, const fcomplex *transvec, fcomplex *sig2);

void FNAME(h2dmplochf)(const fint *nd, const fcomplex *zk,
    const double *rscale1, const double *center1, const fcomplex *sig,
    const fint *nterms1, const double *rscale2, const double *center2,
    fcomplex *sig2, const fint *nterms2, const fint *nsig,
    fcomplex *wsave, const fcomplex *tvec);

void FNAME(h2d_sig2exp)(const fint *nd, const fint *nsig,
    fcomplex *sig, fcomplex *wsave, const fint *nterms,
    fcomplex *expans);

void FNAME(h2d_mkm2ltrans)(const fcomplex *zk, const double *center1,
    const fint *nterms1, const double *center2, const fint *nterms2,
    const fint *nsig, fcomplex *wsave, fcomplex *transvec);

#endif
