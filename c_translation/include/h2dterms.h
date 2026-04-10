/*
 * h2dterms.h - C translation of src/helmholtz/h2dterms.f
 */

#ifndef FMM2D_H2DTERMS_H
#define FMM2D_H2DTERMS_H

#include "fmm2d_c.h"

void FNAME(h2dterms)(const double *bsize, const fcomplex *zk,
                     const double *eps, fint *nterms, fint *ier);

#endif
