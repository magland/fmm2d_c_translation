/*
 * cauchykernels2d.h - C translation of src/laplace/cauchykernels2d.f
 *
 * Direct (kernel) evaluation routines for the 2D Laplace FMM, using
 * the unscaled log response and complex 1/z dipole formulation.
 *
 * Each routine INCREMENTS its output arrays (pot, grad, hess); it
 * does not overwrite them. Near-source contributions are skipped
 * based on a `thresh` distance check (the precise comparison varies
 * between routines and is preserved exactly from the Fortran source).
 */

#ifndef FMM2D_CAUCHYKERNELS2D_H
#define FMM2D_CAUCHYKERNELS2D_H

#include "fmm2d_c.h"

void FNAME(c2d_directcp)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *charge, const double *targ, const fint *nt,
                         fcomplex *pot, const double *thresh);

void FNAME(c2d_directcg)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *charge, const double *targ, const fint *nt,
                         fcomplex *pot, fcomplex *grad, const double *thresh);

void FNAME(c2d_directch)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *charge, const double *targ, const fint *nt,
                         fcomplex *pot, fcomplex *grad, fcomplex *hess,
                         const double *thresh);

void FNAME(c2d_directdp)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *dipstr, const double *targ, const fint *nt,
                         fcomplex *pot, const double *thresh);

void FNAME(c2d_directdg)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *dipstr, const double *targ, const fint *nt,
                         fcomplex *pot, fcomplex *grad, const double *thresh);

void FNAME(c2d_directdh)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *dipstr, const double *targ, const fint *nt,
                         fcomplex *pot, fcomplex *grad, fcomplex *hess,
                         const double *thresh);

void FNAME(c2d_directcdp)(const fint *nd, const double *sources, const fint *ns,
                          const fcomplex *charge, const fcomplex *dipstr,
                          const double *targ, const fint *nt,
                          fcomplex *pot, const double *thresh);

void FNAME(c2d_directcdg)(const fint *nd, const double *sources, const fint *ns,
                          const fcomplex *charge, const fcomplex *dipstr,
                          const double *targ, const fint *nt,
                          fcomplex *pot, fcomplex *grad, const double *thresh);

void FNAME(c2d_directcdh)(const fint *nd, const double *sources, const fint *ns,
                          const fcomplex *charge, const fcomplex *dipstr,
                          const double *targ, const fint *nt,
                          fcomplex *pot, fcomplex *grad, fcomplex *hess,
                          const double *thresh);

#endif /* FMM2D_CAUCHYKERNELS2D_H */
