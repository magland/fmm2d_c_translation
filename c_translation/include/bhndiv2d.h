/*
 * bhndiv2d.h - C translation of src/biharmonic/bhndiv2d.f
 *
 * Subdivision criterion estimator for the biharmonic 2D FMM.
 * Same eps thresholds as hndiv2d but separate routine for the
 * biharmonic call graph.
 */

#ifndef FMM2D_BHNDIV2D_H
#define FMM2D_BHNDIV2D_H

#include "fmm2d_c.h"

void FNAME(bhndiv2d)(
    const double *eps,
    const fint *ns, const fint *nt,
    const fint *ifcharge, const fint *ifdipole,
    const fint *ifpgh, const fint *ifpghtarg,
    fint *ndiv, fint *idivflag);

#endif /* FMM2D_BHNDIV2D_H */
