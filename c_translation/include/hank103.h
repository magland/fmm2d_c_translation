/*
 * hank103.h - C translation of src/common/hank103.f
 */

#ifndef FMM2D_HANK103_H
#define FMM2D_HANK103_H

#include "fmm2d_c.h"

void FNAME(hank103)(const fcomplex *z, fcomplex *h0, fcomplex *h1,
                    const fint *ifexpon);

void FNAME(hank103u)(const fcomplex *z, fint *ier, fcomplex *h0,
                     fcomplex *h1, const fint *ifexpon);

void FNAME(hank103r)(const fcomplex *z, fint *ier, fcomplex *h0,
                     fcomplex *h1, const fint *ifexpon);

void FNAME(hank103p)(const fcomplex *p, const fint *m, const fcomplex *z,
                     fcomplex *f);

void FNAME(hank103l)(const fcomplex *z, fcomplex *h0, fcomplex *h1,
                     const fint *ifexpon);

void FNAME(hank103a)(const fcomplex *z, fcomplex *h0, fcomplex *h1,
                     const fint *ifexpon);

void FNAME(hanks103)(const fcomplex *z, fcomplex *hanks, const fint *n,
                     const fint *ifexpon);

void FNAME(hanks104)(const fcomplex *z, fcomplex *hanks, const fint *n,
                     const fint *ifexpon);

#endif
