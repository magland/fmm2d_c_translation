/*
 * stokkernels2d.h - C translation of src/stokes/stokkernels2d.f
 *
 * Direct evaluation kernels for the 2D Stokes FMM.
 */

#ifndef FMM2D_STOKKERNELS2D_H
#define FMM2D_STOKKERNELS2D_H

#include "fmm2d_c.h"

void FNAME(st2ddirectstokg)(const fint *nd, const double *sources,
                            const double *stoklet, const fint *ns,
                            const double *targ, const fint *nt,
                            double *pot, double *pre, double *grad,
                            const double *thresh);

void FNAME(st2ddirectstokstrsg)(const fint *nd, const double *sources,
                                const fint *ifstoklet, const double *stoklet,
                                const fint *istress, const double *strslet,
                                const double *strsvec,
                                const fint *ns, const double *targ, const fint *nt,
                                double *pot, double *pre, double *grad,
                                const double *thresh);

#endif /* FMM2D_STOKKERNELS2D_H */
