/*
 * hfmm2d.h - C translation of src/helmholtz/hfmm2d.f
 */

#ifndef FMM2D_HFMM2D_H
#define FMM2D_HFMM2D_H

#include "fmm2d_c.h"

void FNAME(hfmm2d)(const fint *nd, const double *eps, const fcomplex *zk,
    const fint *ns, const double *sources,
    const fint *ifcharge, const fcomplex *charge,
    const fint *ifdipole, const fcomplex *dipstr, const double *dipvec,
    fint *iper, const fint *ifpgh, fcomplex *pot, fcomplex *grad,
    fcomplex *hess, const fint *nt, const double *targ,
    const fint *ifpghtarg, fcomplex *pottarg, fcomplex *gradtarg,
    fcomplex *hesstarg, fint *ier);

void FNAME(h2dmpalloc)(const fint *nd, const fint *laddr, fint *iaddr,
    const fint *nlevels, fint *lmptot, const fint *nterms);

#endif
