/*
 * lapkernels2d.h - C translation of src/laplace/lapkernels2d.f
 *
 * Direct evaluation kernels for the 2D Laplace FMM with complex-valued
 * charge/dipstr/pot/grad/hess and real-valued dipvec.
 */

#ifndef FMM2D_LAPKERNELS2D_H
#define FMM2D_LAPKERNELS2D_H

#include "fmm2d_c.h"

void FNAME(l2d_directcp)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *charge, const double *targ, const fint *nt,
                         fcomplex *pot, const double *thresh);

void FNAME(l2d_directcg)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *charge, const double *targ, const fint *nt,
                         fcomplex *pot, fcomplex *grad, const double *thresh);

void FNAME(l2d_directch)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *charge, const double *targ, const fint *nt,
                         fcomplex *pot, fcomplex *grad, fcomplex *hess,
                         const double *thresh);

void FNAME(l2d_directdp)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *dipstr, const double *dipvec,
                         const double *targ, const fint *nt,
                         fcomplex *pot, const double *thresh);

void FNAME(l2d_directdg)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *dipstr, const double *dipvec,
                         const double *targ, const fint *nt,
                         fcomplex *pot, fcomplex *grad, const double *thresh);

void FNAME(l2d_directdh)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *dipstr, const double *dipvec,
                         const double *targ, const fint *nt,
                         fcomplex *pot, fcomplex *grad, fcomplex *hess,
                         const double *thresh);

void FNAME(l2d_directcdp)(const fint *nd, const double *sources, const fint *ns,
                          const fcomplex *charge, const fcomplex *dipstr,
                          const double *dipvec,
                          const double *targ, const fint *nt,
                          fcomplex *pot, const double *thresh);

void FNAME(l2d_directcdg)(const fint *nd, const double *sources, const fint *ns,
                          const fcomplex *charge, const fcomplex *dipstr,
                          const double *dipvec,
                          const double *targ, const fint *nt,
                          fcomplex *pot, fcomplex *grad, const double *thresh);

void FNAME(l2d_directcdh)(const fint *nd, const double *sources, const fint *ns,
                          const fcomplex *charge, const fcomplex *dipstr,
                          const double *dipvec,
                          const double *targ, const fint *nt,
                          fcomplex *pot, fcomplex *grad, fcomplex *hess,
                          const double *thresh);

#endif /* FMM2D_LAPKERNELS2D_H */
