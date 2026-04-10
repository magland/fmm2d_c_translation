#ifndef MBHGREEN2D_H
#define MBHGREEN2D_H

#include "fmm2d_c.h"

void FNAME(modbhgreen_all)(
    const double *beta, const double *zx, const double *zy,
    const fint *ifpot, double *pot,
    const fint *ifgrad, double *grad,
    const fint *ifhess, double *hess,
    const fint *ifder3, double *der3,
    const fint *ifder4, double *der4,
    const fint *ifder5, double *der5);

void FNAME(modbhgreen)(
    const double *beta, const double *zx, const double *zy,
    const fint *ifpot, double *pot,
    const fint *ifgrad, double *grad,
    const fint *ifhess, double *hess);

void FNAME(modbhgreend1)(
    const double *beta, const double *zx, const double *zy,
    const fint *ifpot, double *pot,
    const fint *ifgrad, double *grad,
    const fint *ifhess, double *hess,
    const double *dir1);

void FNAME(modbhgreend2)(
    const double *beta, const double *zx, const double *zy,
    const fint *ifpot, double *pot,
    const fint *ifgrad, double *grad,
    const fint *ifhess, double *hess,
    const double *dir1, const double *dir2);

void FNAME(difflogbk)(
    const double *x, const double *beta,
    const fint *if0, double *g0,
    const fint *if1, double *g1,
    const fint *if2, double *g2,
    const fint *if3, double *g3);

void FNAME(diffslogbk)(
    const double *x, const double *beta,
    const double *rscale, double *diffs,
    const fint *n);

void FNAME(diffslogbk_fast)(
    const double *x, const double *beta,
    const double *rscale, double *diffs,
    const fint *ifders, double *ders, double *kvec,
    const fint *n);

void FNAME(diffszkik)(
    const double *x, const double *beta,
    const double *rscale, double *diffs,
    const fint *n);

void FNAME(diffszkik_fast)(
    const double *x, const double *beta,
    const double *rscale, double *diffs,
    const fint *ifders, double *ders, double *ivec,
    const fint *n);

#endif /* MBHGREEN2D_H */
