/*
 * helmkernels2d.h - C translation of src/helmholtz/helmkernels2d.f
 */

#ifndef FMM2D_HELMKERNELS2D_H
#define FMM2D_HELMKERNELS2D_H

#include "fmm2d_c.h"

void FNAME(h2d_directcp)(const fint *nd, const fcomplex *wavek,
    const double *sources, const fint *ns, const fcomplex *charge,
    const double *targ, const fint *nt, fcomplex *pot,
    const double *thresh);

void FNAME(h2d_directcg)(const fint *nd, const fcomplex *wavek,
    const double *sources, const fint *ns, const fcomplex *charge,
    const double *targ, const fint *nt, fcomplex *pot,
    fcomplex *grad, const double *thresh);

void FNAME(h2d_directch)(const fint *nd, const fcomplex *wavek,
    const double *sources, const fint *ns, const fcomplex *charge,
    const double *targ, const fint *nt, fcomplex *pot,
    fcomplex *grad, fcomplex *hess, const double *thresh);

void FNAME(h2d_directdp)(const fint *nd, const fcomplex *wavek,
    const double *sources, const fint *ns, const fcomplex *dipstr,
    const double *dipvec, const double *targ, const fint *nt,
    fcomplex *pot, const double *thresh);

void FNAME(h2d_directdg)(const fint *nd, const fcomplex *wavek,
    const double *sources, const fint *ns, const fcomplex *dipstr,
    const double *dipvec, const double *targ, const fint *nt,
    fcomplex *pot, fcomplex *grad, const double *thresh);

void FNAME(h2d_directdh)(const fint *nd, const fcomplex *wavek,
    const double *sources, const fint *ns, const fcomplex *dipstr,
    const double *dipvec, const double *targ, const fint *nt,
    fcomplex *pot, fcomplex *grad, fcomplex *hess,
    const double *thresh);

void FNAME(h2d_directcdp)(const fint *nd, const fcomplex *wavek,
    const double *sources, const fint *ns, const fcomplex *charge,
    const fcomplex *dipstr, const double *dipvec,
    const double *targ, const fint *nt, fcomplex *pot,
    const double *thresh);

void FNAME(h2d_directcdg)(const fint *nd, const fcomplex *wavek,
    const double *sources, const fint *ns, const fcomplex *charge,
    const fcomplex *dipstr, const double *dipvec,
    const double *targ, const fint *nt, fcomplex *pot,
    fcomplex *grad, const double *thresh);

void FNAME(h2d_directcdh)(const fint *nd, const fcomplex *wavek,
    const double *sources, const fint *ns, const fcomplex *charge,
    const fcomplex *dipstr, const double *dipvec,
    const double *targ, const fint *nt, fcomplex *pot,
    fcomplex *grad, fcomplex *hess, const double *thresh);

#endif
