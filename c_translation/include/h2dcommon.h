/*
 * h2dcommon.h - C translation of src/helmholtz/h2dcommon.f
 */

#ifndef FMM2D_H2DCOMMON_H
#define FMM2D_H2DCOMMON_H

#include "fmm2d_c.h"

void FNAME(h2cart2polar)(const double *zat, double *r, double *theta);

void FNAME(h2dall)(const fint *nterms, const fcomplex *z,
                   const double *rscale, fcomplex *hvec,
                   const fint *ifder, fcomplex *hder);

void FNAME(h2dmpzero)(const fint *nd, fcomplex *mpole, const fint *nterms);

void FNAME(h2dsigzero)(const fint *nd, fcomplex *sig, const fint *nsig);

#endif
