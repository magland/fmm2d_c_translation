/*
 * hfmm2d_ndiv.h - C translation of src/helmholtz/hfmm2d_ndiv.f
 */

#ifndef FMM2D_HFMM2D_NDIV_H
#define FMM2D_HFMM2D_NDIV_H

#include "fmm2d_c.h"

void FNAME(hfmm2d_ndiv)(const fint *nd, const double *eps,
    const fcomplex *zk, const fint *ns, const double *sources,
    const fint *ifcharge, const fcomplex *charge,
    const fint *ifdipole, const fcomplex *dipstr, const double *dipvec,
    fint *iper, const fint *ifpgh, fcomplex *pot, fcomplex *grad,
    fcomplex *hess, const fint *nt, const double *targ,
    const fint *ifpghtarg, fcomplex *pottarg, fcomplex *gradtarg,
    fcomplex *hesstarg, const fint *ndiv, const fint *idivflag,
    const fint *ifnear, double *timeinfo, fint *ier);

#endif
