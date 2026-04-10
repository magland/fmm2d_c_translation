/*
 * rlapkernels2d.h - C translation of src/laplace/rlapkernels2d.f
 *
 * Direct evaluation kernels for the 2D Laplace FMM with real-valued
 * charge/dipstr/dipvec/pot/grad/hess.
 */

#ifndef FMM2D_RLAPKERNELS2D_H
#define FMM2D_RLAPKERNELS2D_H

#include "fmm2d_c.h"

void FNAME(r2d_directcp)(const fint *nd, const double *sources, const fint *ns,
                         const double *charge, const double *targ, const fint *nt,
                         double *pot, const double *thresh);

void FNAME(r2d_directcg)(const fint *nd, const double *sources, const fint *ns,
                         const double *charge, const double *targ, const fint *nt,
                         double *pot, double *grad, const double *thresh);

void FNAME(r2d_directch)(const fint *nd, const double *sources, const fint *ns,
                         const double *charge, const double *targ, const fint *nt,
                         double *pot, double *grad, double *hess,
                         const double *thresh);

void FNAME(r2d_directdp)(const fint *nd, const double *sources, const fint *ns,
                         const double *dipstr, const double *dipvec,
                         const double *targ, const fint *nt,
                         double *pot, const double *thresh);

void FNAME(r2d_directdg)(const fint *nd, const double *sources, const fint *ns,
                         const double *dipstr, const double *dipvec,
                         const double *targ, const fint *nt,
                         double *pot, double *grad, const double *thresh);

void FNAME(r2d_directdh)(const fint *nd, const double *sources, const fint *ns,
                         const double *dipstr, const double *dipvec,
                         const double *targ, const fint *nt,
                         double *pot, double *grad, double *hess,
                         const double *thresh);

void FNAME(r2d_directcdp)(const fint *nd, const double *sources, const fint *ns,
                          const double *charge, const double *dipstr,
                          const double *dipvec,
                          const double *targ, const fint *nt,
                          double *pot, const double *thresh);

void FNAME(r2d_directcdg)(const fint *nd, const double *sources, const fint *ns,
                          const double *charge, const double *dipstr,
                          const double *dipvec,
                          const double *targ, const fint *nt,
                          double *pot, double *grad, const double *thresh);

void FNAME(r2d_directcdh)(const fint *nd, const double *sources, const fint *ns,
                          const double *charge, const double *dipstr,
                          const double *dipvec,
                          const double *targ, const fint *nt,
                          double *pot, double *grad, double *hess,
                          const double *thresh);

#endif /* FMM2D_RLAPKERNELS2D_H */
