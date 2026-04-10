/*
 * cdjseval2d.h - C translation of src/common/cdjseval2d.f
 */

#ifndef FMM2D_CDJSEVAL2D_H
#define FMM2D_CDJSEVAL2D_H

#include "fmm2d_c.h"

void FNAME(jbessel2d)(const fint *nterms, const fcomplex *z,
                      const double *rscale, fcomplex *fjs,
                      const fint *ifder, fcomplex *fjder);

#endif
