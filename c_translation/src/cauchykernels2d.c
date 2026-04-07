/*
 * cauchykernels2d.c - C translation of src/laplace/cauchykernels2d.f
 *
 * Direct evaluation kernels for the 2D Laplace FMM. Each routine
 * INCREMENTS its output arrays (pot, grad, hess) at every target due
 * to charges and/or dipoles at every source. The unscaled log
 * response (log|z|) is used; gradients and Hessians are derivatives
 * d/dz, d^2/dz^2 of the complex potential.
 *
 * The threshold check that skips near-source contributions varies
 * between routines: some compare rr (squared distance) against
 * thresh*thresh with .le., some with .lt., and the dipole-only
 * routines compare cabs(z) against thresh directly. These are
 * preserved exactly from the Fortran source.
 *
 * OpenMP pragmas (cf2py / c$omp) in the Fortran source are comments
 * and are ignored by this translation; the resulting code runs
 * sequentially but produces identical numerical output (operations
 * proceed in the same order as the serial Fortran loops).
 */

#include "cauchykernels2d.h"

/*
 * c2d_directcp: charges -> potential
 *
 * Threshold: rr .le. thresh2 (squared)
 *
 * pot(ii,j) += log(|z|) * charge(ii,i)
 */
void FNAME(c2d_directcp)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *charge, const double *targ, const fint *nt,
                         fcomplex *pot, const double *thresh)
{
    fint i, j, ii;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    double xdiff, ydiff, rr, rtmp, thresh2;
    /* unused locals preserved from the Fortran declarations */
    double r;
    fcomplex z, ima4inv, ztmp;
    (void)r;
    (void)z;
    (void)ima4inv;
    (void)ztmp;

    thresh2 = (*thresh) * (*thresh);

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1, j, 2)] - sources[FA2(1, i, 2)];
            ydiff = targ[FA2(2, j, 2)] - sources[FA2(2, i, 2)];
            rr = xdiff * xdiff + ydiff * ydiff;
            if (rr <= thresh2) goto skip;

            rtmp = log(rr) / 2.0;
            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii, j, nd_v)] += rtmp * charge[FA2(ii, i, nd_v)];
            }
        skip:;
        }
    }
}

/*
 * c2d_directcg: charges -> potential + gradient
 *
 * Threshold: rr .lt. thresh2 (squared, strict)
 *
 * pot  += log(|z|) * charge
 * grad += (1/z)    * charge
 */
void FNAME(c2d_directcg)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *charge, const double *targ, const fint *nt,
                         fcomplex *pot, fcomplex *grad, const double *thresh)
{
    fint i, j, ii;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    double xdiff, ydiff, rr, rtmp, thresh2;
    double r;
    fcomplex z, zinv, ztmp;
    (void)r;
    (void)ztmp;

    thresh2 = (*thresh) * (*thresh);

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1, j, 2)] - sources[FA2(1, i, 2)];
            ydiff = targ[FA2(2, j, 2)] - sources[FA2(2, i, 2)];
            rr = xdiff * xdiff + ydiff * ydiff;
            if (rr < thresh2) goto skip;
            rtmp = log(rr) / 2.0;
            z = xdiff + ydiff * I;
            zinv = 1.0 / z;

            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii, j, nd_v)]  += rtmp * charge[FA2(ii, i, nd_v)];
                grad[FA2(ii, j, nd_v)] += zinv * charge[FA2(ii, i, nd_v)];
            }
        skip:;
        }
    }
}

/*
 * c2d_directch: charges -> potential + gradient + Hessian
 *
 * Threshold: rr .lt. thresh2 (squared, strict)
 *
 * pot  += log(|z|)   * charge
 * grad += (1/z)      * charge
 * hess += -(1/z^2)   * charge
 */
void FNAME(c2d_directch)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *charge, const double *targ, const fint *nt,
                         fcomplex *pot, fcomplex *grad, fcomplex *hess,
                         const double *thresh)
{
    fint i, j, ii;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    fint ifexpon;
    double xdiff, ydiff, rr, rtmp, thresh2;
    double r;
    fcomplex z, zinv, ztmp;
    (void)r;
    (void)ifexpon;

    thresh2 = (*thresh) * (*thresh);

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1, j, 2)] - sources[FA2(1, i, 2)];
            ydiff = targ[FA2(2, j, 2)] - sources[FA2(2, i, 2)];
            rr = xdiff * xdiff + ydiff * ydiff;
            if (rr < thresh2) goto skip;
            rtmp = log(rr) / 2.0;
            z = xdiff + ydiff * I;
            zinv = 1.0 / z;
            ztmp = -zinv * zinv;

            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii, j, nd_v)]  += rtmp * charge[FA2(ii, i, nd_v)];
                grad[FA2(ii, j, nd_v)] += zinv * charge[FA2(ii, i, nd_v)];
                hess[FA2(ii, j, nd_v)] += ztmp * charge[FA2(ii, i, nd_v)];
            }
        skip:;
        }
    }
}

/*
 * c2d_directdp: dipoles -> potential
 *
 * Threshold: cabs(z) .le. thresh (norm, not squared)
 *
 * pot += dipstr / z
 */
void FNAME(c2d_directdp)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *dipstr, const double *targ, const fint *nt,
                         fcomplex *pot, const double *thresh)
{
    fint i, j, ii;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    double xdiff, ydiff, rr, rtmp;
    double r;
    fcomplex z, zinv;
    (void)rr;
    (void)rtmp;
    (void)r;

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1, j, 2)] - sources[FA2(1, i, 2)];
            ydiff = targ[FA2(2, j, 2)] - sources[FA2(2, i, 2)];
            z = xdiff + ydiff * I;
            if (cabs(z) <= *thresh) goto skip;

            zinv = 1.0 / z;
            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii, j, nd_v)] += zinv * dipstr[FA2(ii, i, nd_v)];
            }
        skip:;
        }
    }
}

/*
 * c2d_directdg: dipoles -> potential + gradient
 *
 * Threshold: cabs(z) .lt. thresh (norm, strict)
 *
 * pot  += dipstr / z
 * grad += -dipstr / z^2
 */
void FNAME(c2d_directdg)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *dipstr, const double *targ, const fint *nt,
                         fcomplex *pot, fcomplex *grad, const double *thresh)
{
    fint i, j, ii;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    double xdiff, ydiff, rr, rtmp;
    double r;
    fcomplex z, zinv, ztmp2;
    (void)rr;
    (void)rtmp;
    (void)r;

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1, j, 2)] - sources[FA2(1, i, 2)];
            ydiff = targ[FA2(2, j, 2)] - sources[FA2(2, i, 2)];
            z = xdiff + ydiff * I;
            if (cabs(z) < *thresh) goto skip;
            zinv = 1.0 / z;
            ztmp2 = -zinv * zinv;

            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii, j, nd_v)]  += zinv  * dipstr[FA2(ii, i, nd_v)];
                grad[FA2(ii, j, nd_v)] += ztmp2 * dipstr[FA2(ii, i, nd_v)];
            }
        skip:;
        }
    }
}

/*
 * c2d_directdh: dipoles -> potential + gradient + Hessian
 *
 * Threshold: cabs(z) .lt. thresh (norm, strict)
 *
 * pot  += dipstr / z
 * grad += -dipstr / z^2
 * hess += 2 * dipstr / z^3
 */
void FNAME(c2d_directdh)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *dipstr, const double *targ, const fint *nt,
                         fcomplex *pot, fcomplex *grad, fcomplex *hess,
                         const double *thresh)
{
    fint i, j, ii;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    fint ifexpon;
    double xdiff, ydiff, rr, rtmp;
    double r;
    fcomplex z, zinv, ztmp2, ztmp3;
    (void)rr;
    (void)rtmp;
    (void)r;
    (void)ifexpon;

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1, j, 2)] - sources[FA2(1, i, 2)];
            ydiff = targ[FA2(2, j, 2)] - sources[FA2(2, i, 2)];
            z = xdiff + ydiff * I;
            if (cabs(z) < *thresh) goto skip;
            zinv = 1.0 / z;
            ztmp2 = -zinv * zinv;
            ztmp3 = -2.0 * ztmp2 * zinv;

            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii, j, nd_v)]  += zinv  * dipstr[FA2(ii, i, nd_v)];
                grad[FA2(ii, j, nd_v)] += ztmp2 * dipstr[FA2(ii, i, nd_v)];
                hess[FA2(ii, j, nd_v)] += ztmp3 * dipstr[FA2(ii, i, nd_v)];
            }
        skip:;
        }
    }
}

/*
 * c2d_directcdp: charges + dipoles -> potential
 *
 * Threshold: rr .le. thresh2 (squared)
 *
 * pot += log(|z|) * charge + dipstr / z
 */
void FNAME(c2d_directcdp)(const fint *nd, const double *sources, const fint *ns,
                          const fcomplex *charge, const fcomplex *dipstr,
                          const double *targ, const fint *nt,
                          fcomplex *pot, const double *thresh)
{
    fint i, j, ii;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    double xdiff, ydiff, rr, rtmp, thresh2;
    double r;
    fcomplex z, zinv, ztmp;
    (void)r;
    (void)ztmp;

    thresh2 = (*thresh) * (*thresh);

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1, j, 2)] - sources[FA2(1, i, 2)];
            ydiff = targ[FA2(2, j, 2)] - sources[FA2(2, i, 2)];
            rr = xdiff * xdiff + ydiff * ydiff;
            if (rr <= thresh2) goto skip;

            rtmp = log(rr) / 2.0;
            z = xdiff + ydiff * I;
            zinv = 1.0 / z;
            for (ii = 1; ii <= nd_v; ii++) {
                /* Fortran:
                 *   pot(ii,j) = pot(ii,j) + rtmp*charge(ii,i)
                 *                         + zinv*dipstr(ii,i)
                 * is evaluated left-to-right as
                 *   ((pot + rtmp*charge) + zinv*dipstr).
                 * Writing this as `pot += a + b` in C parses as
                 *   pot + (a + b),
                 * which has a different rounding order. To preserve
                 * bit-for-bit numerics, accumulate one term at a time.
                 */
                pot[FA2(ii, j, nd_v)] += rtmp * charge[FA2(ii, i, nd_v)];
                pot[FA2(ii, j, nd_v)] += zinv * dipstr[FA2(ii, i, nd_v)];
            }
        skip:;
        }
    }
}

/*
 * c2d_directcdg: charges + dipoles -> potential + gradient
 *
 * Threshold: rr .lt. thresh2 (squared, strict)
 *
 * pot  += log(|z|) * charge + dipstr / z
 * grad += (1/z)    * charge + (-1/z^2) * dipstr
 */
void FNAME(c2d_directcdg)(const fint *nd, const double *sources, const fint *ns,
                          const fcomplex *charge, const fcomplex *dipstr,
                          const double *targ, const fint *nt,
                          fcomplex *pot, fcomplex *grad, const double *thresh)
{
    fint i, j, ii;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    double xdiff, ydiff, rr, rtmp, thresh2;
    double r;
    fcomplex z, zinv, ztmp, ztmp2;
    (void)r;
    (void)ztmp;

    thresh2 = (*thresh) * (*thresh);

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1, j, 2)] - sources[FA2(1, i, 2)];
            ydiff = targ[FA2(2, j, 2)] - sources[FA2(2, i, 2)];
            rr = xdiff * xdiff + ydiff * ydiff;
            if (rr < thresh2) goto skip;
            rtmp = log(rr) / 2.0;
            z = xdiff + ydiff * I;
            zinv = 1.0 / z;
            ztmp2 = -zinv * zinv;

            for (ii = 1; ii <= nd_v; ii++) {
                /* Fortran left-to-right associativity for `pot + a + b`
                 * forces `(pot + a) + b`. Splitting the C accumulations
                 * into two += statements preserves that order; the
                 * `pot += a + b` shorthand groups as `pot + (a + b)`
                 * and gives a different last-bit result. */
                pot[FA2(ii, j, nd_v)]  += rtmp  * charge[FA2(ii, i, nd_v)];
                pot[FA2(ii, j, nd_v)]  += zinv  * dipstr[FA2(ii, i, nd_v)];
                grad[FA2(ii, j, nd_v)] += zinv  * charge[FA2(ii, i, nd_v)];
                grad[FA2(ii, j, nd_v)] += ztmp2 * dipstr[FA2(ii, i, nd_v)];
            }
        skip:;
        }
    }
}

/*
 * c2d_directcdh: charges + dipoles -> potential + gradient + Hessian
 *
 * Threshold: rr .lt. thresh2 (squared, strict)
 *
 * pot  += log(|z|)   * charge + dipstr / z
 * grad += (1/z)      * charge + (-1/z^2) * dipstr
 * hess += (-1/z^2)   * charge + (2/z^3)  * dipstr
 */
void FNAME(c2d_directcdh)(const fint *nd, const double *sources, const fint *ns,
                          const fcomplex *charge, const fcomplex *dipstr,
                          const double *targ, const fint *nt,
                          fcomplex *pot, fcomplex *grad, fcomplex *hess,
                          const double *thresh)
{
    fint i, j, ii;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    fint ifexpon;
    double xdiff, ydiff, rr, rtmp, thresh2;
    double r;
    fcomplex z, zinv, ztmp2, ztmp3;
    (void)r;
    (void)ifexpon;

    thresh2 = (*thresh) * (*thresh);

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1, j, 2)] - sources[FA2(1, i, 2)];
            ydiff = targ[FA2(2, j, 2)] - sources[FA2(2, i, 2)];
            rr = xdiff * xdiff + ydiff * ydiff;
            if (rr < thresh2) goto skip;
            rtmp = log(rr) / 2.0;
            z = xdiff + ydiff * I;
            zinv = 1.0 / z;
            ztmp2 = -zinv * zinv;
            ztmp3 = -2.0 * ztmp2 * zinv;

            for (ii = 1; ii <= nd_v; ii++) {
                /* See c2d_directcdg comment: split into separate
                 * += statements to mirror Fortran's left-to-right
                 * (pot + a) + b ordering. */
                pot[FA2(ii, j, nd_v)]  += rtmp  * charge[FA2(ii, i, nd_v)];
                pot[FA2(ii, j, nd_v)]  += zinv  * dipstr[FA2(ii, i, nd_v)];
                grad[FA2(ii, j, nd_v)] += zinv  * charge[FA2(ii, i, nd_v)];
                grad[FA2(ii, j, nd_v)] += ztmp2 * dipstr[FA2(ii, i, nd_v)];
                hess[FA2(ii, j, nd_v)] += ztmp2 * charge[FA2(ii, i, nd_v)];
                hess[FA2(ii, j, nd_v)] += ztmp3 * dipstr[FA2(ii, i, nd_v)];
            }
        skip:;
        }
    }
}
