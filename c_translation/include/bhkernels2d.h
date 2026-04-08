/*
 * bhkernels2d.h - C translation of src/biharmonic/bhkernels2d.f
 *
 * Direct (kernel) evaluation routines for the 2D biharmonic FMM.
 *
 * The complex velocity at zt = xt + i*yt due to a charge with complex
 * components (c1, c2) at zs is
 *
 *     vel = 2*c1*log|zt-zs| + c2*(zt-zs)/conj(zt-zs)
 *
 * and due to a dipole with complex components (d1, d2, d3) is
 *
 *     vel = d1/(zt-zs) + d2*(zt-zs)/conj(zt-zs)^2
 *           + d3/conj(zt-zs)
 *
 * Each routine INCREMENTS its output arrays (vel, grad); it does not
 * overwrite them. Near-source contributions are skipped based on
 * abs(zt-zs) .le. thresh (preserved exactly from the Fortran source).
 *
 * Array shapes (column-major):
 *   sources(2, ns)
 *   targ(2, nt)
 *   charges(nd, 2, ns)
 *   dippar(nd, 3, ns)
 *   vel(nd, nt)
 *   grad(nd, 3, nt)
 */

#ifndef FMM2D_BHKERNELS2D_H
#define FMM2D_BHKERNELS2D_H

#include "fmm2d_c.h"

void FNAME(bh2d_directcp)(const fint *nd, const double *sources,
                          const fint *ns, const fcomplex *charges,
                          const double *targ, const fint *nt,
                          fcomplex *vel, const double *thresh);

void FNAME(bh2d_directcg)(const fint *nd, const double *sources,
                          const fint *ns, const fcomplex *charges,
                          const double *targ, const fint *nt,
                          fcomplex *vel, fcomplex *grad,
                          const double *thresh);

void FNAME(bh2d_directdp)(const fint *nd, const double *sources,
                          const fint *ns, const fcomplex *dippar,
                          const double *targ, const fint *nt,
                          fcomplex *vel, const double *thresh);

void FNAME(bh2d_directdg)(const fint *nd, const double *sources,
                          const fint *ns, const fcomplex *dippar,
                          const double *targ, const fint *nt,
                          fcomplex *vel, fcomplex *grad,
                          const double *thresh);

void FNAME(bh2d_directcdp)(const fint *nd, const double *sources,
                           const fint *ns, const fcomplex *charges,
                           const fcomplex *dippar,
                           const double *targ, const fint *nt,
                           fcomplex *vel, const double *thresh);

void FNAME(bh2d_directcdg)(const fint *nd, const double *sources,
                           const fint *ns, const fcomplex *charges,
                           const fcomplex *dippar,
                           const double *targ, const fint *nt,
                           fcomplex *vel, fcomplex *grad,
                           const double *thresh);

#endif /* FMM2D_BHKERNELS2D_H */
