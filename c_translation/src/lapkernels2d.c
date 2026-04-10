/*
 * lapkernels2d.c - C translation of src/laplace/lapkernels2d.f
 *
 * Direct evaluation kernels for the 2D Laplace FMM with complex-valued
 * charge, dipstr, pot, grad, hess and real-valued dipvec.
 * Each routine INCREMENTS its output arrays.
 * The unscaled log response (log|z|) is used.
 *
 * Threshold comparisons (.le. vs .lt.) are preserved exactly from
 * the Fortran source.
 */

#include "lapkernels2d.h"

/*
 * l2d_directcp: charges -> potential
 * Threshold: rr .le. thresh2
 */
void FNAME(l2d_directcp)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *charge, const double *targ, const fint *nt,
                         fcomplex *pot, const double *thresh)
{
    fint i, j, ii;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    double xdiff, ydiff, rr, rtmp, thresh2;

    thresh2 = (*thresh) * (*thresh);

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1, j, 2)] - sources[FA2(1, i, 2)];
            ydiff = targ[FA2(2, j, 2)] - sources[FA2(2, i, 2)];
            rr = xdiff * xdiff + ydiff * ydiff;
            if (rr <= thresh2) goto skip;

            rtmp = log(rr) / 2;
            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii, j, nd_v)] += rtmp * charge[FA2(ii, i, nd_v)];
            }
        skip:;
        }
    }
}

/*
 * l2d_directcg: charges -> potential + gradient
 * Threshold: rr .lt. thresh2
 */
void FNAME(l2d_directcg)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *charge, const double *targ, const fint *nt,
                         fcomplex *pot, fcomplex *grad, const double *thresh)
{
    fint i, j, ii;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    double xdiff, ydiff, rr, rtmp, thresh2, dx, dy;

    thresh2 = (*thresh) * (*thresh);

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1, j, 2)] - sources[FA2(1, i, 2)];
            ydiff = targ[FA2(2, j, 2)] - sources[FA2(2, i, 2)];
            rr = xdiff * xdiff + ydiff * ydiff;
            if (rr < thresh2) goto skip;
            rtmp = log(rr) / 2;
            dx = xdiff / rr;
            dy = ydiff / rr;

            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii, j, nd_v)] += rtmp * charge[FA2(ii, i, nd_v)];
                grad[FA3(ii, 1, j, nd_v, 2)] += dx * charge[FA2(ii, i, nd_v)];
                grad[FA3(ii, 2, j, nd_v, 2)] += dy * charge[FA2(ii, i, nd_v)];
            }
        skip:;
        }
    }
}

/*
 * l2d_directch: charges -> potential + gradient + hessian
 * Threshold: rr .lt. thresh2
 */
void FNAME(l2d_directch)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *charge, const double *targ, const fint *nt,
                         fcomplex *pot, fcomplex *grad, fcomplex *hess,
                         const double *thresh)
{
    fint i, j, ii;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    double xdiff, ydiff, rr, thresh2;
    double rtmp, xdiff2, ydiff2, rr2, dx, dy, dxx, dxy, dyy;

    thresh2 = (*thresh) * (*thresh);

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1, j, 2)] - sources[FA2(1, i, 2)];
            ydiff = targ[FA2(2, j, 2)] - sources[FA2(2, i, 2)];
            xdiff2 = xdiff * xdiff;
            ydiff2 = ydiff * ydiff;
            rr = xdiff2 + ydiff2;
            rr2 = rr * rr;
            if (rr < thresh2) goto skip;
            rtmp = log(rr) / 2;
            dx = xdiff / rr;
            dy = ydiff / rr;
            dxx = (rr - 2 * xdiff2) / rr2;
            dxy = -2 * xdiff * ydiff / rr2;
            dyy = (rr - 2 * ydiff2) / rr2;

            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii, j, nd_v)] += rtmp * charge[FA2(ii, i, nd_v)];
                grad[FA3(ii, 1, j, nd_v, 2)] += dx * charge[FA2(ii, i, nd_v)];
                grad[FA3(ii, 2, j, nd_v, 2)] += dy * charge[FA2(ii, i, nd_v)];
                hess[FA3(ii, 1, j, nd_v, 3)] += dxx * charge[FA2(ii, i, nd_v)];
                hess[FA3(ii, 2, j, nd_v, 3)] += dxy * charge[FA2(ii, i, nd_v)];
                hess[FA3(ii, 3, j, nd_v, 3)] += dyy * charge[FA2(ii, i, nd_v)];
            }
        skip:;
        }
    }
}

/*
 * l2d_directdp: dipoles -> potential
 * Threshold: rr .le. thresh2
 */
void FNAME(l2d_directdp)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *dipstr, const double *dipvec,
                         const double *targ, const fint *nt,
                         fcomplex *pot, const double *thresh)
{
    fint i, j, ii;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    double xdiff, ydiff, rr, thresh2, p1, p2;

    thresh2 = (*thresh) * (*thresh);

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1, j, 2)] - sources[FA2(1, i, 2)];
            ydiff = targ[FA2(2, j, 2)] - sources[FA2(2, i, 2)];
            rr = xdiff * xdiff + ydiff * ydiff;
            if (rr <= thresh2) goto skip;

            p1 = -xdiff / rr;
            p2 = -ydiff / rr;

            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii, j, nd_v)] += dipstr[FA2(ii, i, nd_v)] *
                    (dipvec[FA3(ii, 1, i, nd_v, 2)] * p1
                   + dipvec[FA3(ii, 2, i, nd_v, 2)] * p2);
            }
        skip:;
        }
    }
}

/*
 * l2d_directdg: dipoles -> potential + gradient
 * Threshold: rr .le. thresh2
 */
void FNAME(l2d_directdg)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *dipstr, const double *dipvec,
                         const double *targ, const fint *nt,
                         fcomplex *pot, fcomplex *grad, const double *thresh)
{
    fint i, j, ii;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    double xdiff, ydiff, rr, thresh2;
    double xdiff2, ydiff2, p1, p2, rr2;
    double dx1, dx2, dy1, dy2;
    fcomplex d1, d2;

    thresh2 = (*thresh) * (*thresh);

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1, j, 2)] - sources[FA2(1, i, 2)];
            ydiff = targ[FA2(2, j, 2)] - sources[FA2(2, i, 2)];
            xdiff2 = xdiff * xdiff;
            ydiff2 = ydiff * ydiff;
            rr = xdiff2 + ydiff2;
            if (rr <= thresh2) goto skip;

            rr2 = rr * rr;

            p1 = -xdiff / rr;
            p2 = -ydiff / rr;

            dx1 = -(rr - 2 * xdiff2) / rr2;
            dx2 = 2 * xdiff * ydiff / rr2;
            dy1 = dx2;
            dy2 = -(rr - 2 * ydiff2) / rr2;

            for (ii = 1; ii <= nd_v; ii++) {
                d1 = dipstr[FA2(ii, i, nd_v)] * dipvec[FA3(ii, 1, i, nd_v, 2)];
                d2 = dipstr[FA2(ii, i, nd_v)] * dipvec[FA3(ii, 2, i, nd_v, 2)];
                /* pot = pot + d1*p1 + d2*p2  (split for FP order) */
                pot[FA2(ii, j, nd_v)] += d1 * p1;
                pot[FA2(ii, j, nd_v)] += d2 * p2;
                /* grad(1) = grad(1) + d1*dx1 + d2*dx2 */
                grad[FA3(ii, 1, j, nd_v, 2)] += d1 * dx1;
                grad[FA3(ii, 1, j, nd_v, 2)] += d2 * dx2;
                /* grad(2) = grad(2) + d1*dy1 + d2*dy2 */
                grad[FA3(ii, 2, j, nd_v, 2)] += d1 * dy1;
                grad[FA3(ii, 2, j, nd_v, 2)] += d2 * dy2;
            }
        skip:;
        }
    }
}

/*
 * l2d_directdh: dipoles -> potential + gradient + hessian
 * Threshold: rr .le. thresh2
 */
void FNAME(l2d_directdh)(const fint *nd, const double *sources, const fint *ns,
                         const fcomplex *dipstr, const double *dipvec,
                         const double *targ, const fint *nt,
                         fcomplex *pot, fcomplex *grad, fcomplex *hess,
                         const double *thresh)
{
    fint i, j, ii;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    double xdiff, ydiff, rr, thresh2;
    double xdiff2, ydiff2, p1, p2;
    double rr2, rr3;
    double dx1, dx2, dy1, dy2;
    double dxx1, dxx2, dxy1, dxy2, dyy1, dyy2;
    fcomplex d1, d2;

    thresh2 = (*thresh) * (*thresh);

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1, j, 2)] - sources[FA2(1, i, 2)];
            ydiff = targ[FA2(2, j, 2)] - sources[FA2(2, i, 2)];
            xdiff2 = xdiff * xdiff;
            ydiff2 = ydiff * ydiff;
            rr = xdiff2 + ydiff2;
            if (rr <= thresh2) goto skip;

            rr2 = rr * rr;
            rr3 = rr2 * rr;

            p1 = -xdiff / rr;
            p2 = -ydiff / rr;

            dx1 = -(rr - 2 * xdiff2) / rr2;
            dx2 = 2 * xdiff * ydiff / rr2;
            dy1 = dx2;
            dy2 = -(rr - 2 * ydiff2) / rr2;

            dxx1 = -2 * xdiff * (xdiff2 - 3 * ydiff2) / rr3;
            dxx2 = 2 * ydiff * (ydiff2 - 3 * xdiff2) / rr3;
            dxy1 = dxx2;
            dxy2 = -dxx1;
            dyy1 = dxy2;
            dyy2 = -dxx2;

            for (ii = 1; ii <= nd_v; ii++) {
                d1 = dipstr[FA2(ii, i, nd_v)] * dipvec[FA3(ii, 1, i, nd_v, 2)];
                d2 = dipstr[FA2(ii, i, nd_v)] * dipvec[FA3(ii, 2, i, nd_v, 2)];
                /* pot = pot + d1*p1 + d2*p2 */
                pot[FA2(ii, j, nd_v)] += d1 * p1;
                pot[FA2(ii, j, nd_v)] += d2 * p2;
                /* grad = grad + d1*dk1 + d2*dk2 */
                grad[FA3(ii, 1, j, nd_v, 2)] += d1 * dx1;
                grad[FA3(ii, 1, j, nd_v, 2)] += d2 * dx2;
                grad[FA3(ii, 2, j, nd_v, 2)] += d1 * dy1;
                grad[FA3(ii, 2, j, nd_v, 2)] += d2 * dy2;
                /* hess = hess + d1*dkk1 + d2*dkk2 */
                hess[FA3(ii, 1, j, nd_v, 3)] += d1 * dxx1;
                hess[FA3(ii, 1, j, nd_v, 3)] += d2 * dxx2;
                hess[FA3(ii, 2, j, nd_v, 3)] += d1 * dxy1;
                hess[FA3(ii, 2, j, nd_v, 3)] += d2 * dxy2;
                hess[FA3(ii, 3, j, nd_v, 3)] += d1 * dyy1;
                hess[FA3(ii, 3, j, nd_v, 3)] += d2 * dyy2;
            }
        skip:;
        }
    }
}

/*
 * l2d_directcdp: charges + dipoles -> potential
 * Threshold: rr .le. thresh2
 */
void FNAME(l2d_directcdp)(const fint *nd, const double *sources, const fint *ns,
                          const fcomplex *charge, const fcomplex *dipstr,
                          const double *dipvec,
                          const double *targ, const fint *nt,
                          fcomplex *pot, const double *thresh)
{
    fint i, j, ii;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    double xdiff, ydiff, rr, rtmp, thresh2, p1, p2;

    thresh2 = (*thresh) * (*thresh);

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1, j, 2)] - sources[FA2(1, i, 2)];
            ydiff = targ[FA2(2, j, 2)] - sources[FA2(2, i, 2)];
            rr = xdiff * xdiff + ydiff * ydiff;
            if (rr <= thresh2) goto skip;

            rtmp = log(rr) / 2.0;

            p1 = -xdiff / rr;
            p2 = -ydiff / rr;

            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii, j, nd_v)] += rtmp * charge[FA2(ii, i, nd_v)];
                pot[FA2(ii, j, nd_v)] += dipstr[FA2(ii, i, nd_v)] *
                    (dipvec[FA3(ii, 1, i, nd_v, 2)] * p1
                   + dipvec[FA3(ii, 2, i, nd_v, 2)] * p2);
            }
        skip:;
        }
    }
}

/*
 * l2d_directcdg: charges + dipoles -> potential + gradient
 * Threshold: rr .lt. thresh2
 */
void FNAME(l2d_directcdg)(const fint *nd, const double *sources, const fint *ns,
                          const fcomplex *charge, const fcomplex *dipstr,
                          const double *dipvec,
                          const double *targ, const fint *nt,
                          fcomplex *pot, fcomplex *grad, const double *thresh)
{
    fint i, j, ii;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    double xdiff, ydiff, rr, rtmp, thresh2;
    double xdiff2, ydiff2, p1, p2;
    double rr2, dx, dy, dx1, dx2, dy1, dy2;
    fcomplex d1, d2;

    thresh2 = (*thresh) * (*thresh);

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1, j, 2)] - sources[FA2(1, i, 2)];
            ydiff = targ[FA2(2, j, 2)] - sources[FA2(2, i, 2)];
            xdiff2 = xdiff * xdiff;
            ydiff2 = ydiff * ydiff;
            rr = xdiff2 + ydiff2;
            if (rr < thresh2) goto skip;
            rtmp = log(rr) / 2;

            rr2 = rr * rr;

            dx = xdiff / rr;
            dy = ydiff / rr;

            p1 = -dx;
            p2 = -dy;

            dx1 = -(rr - 2 * xdiff2) / rr2;
            dx2 = 2 * xdiff * ydiff / rr2;
            dy1 = dx2;
            dy2 = -(rr - 2 * ydiff2) / rr2;

            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii, j, nd_v)] += rtmp * charge[FA2(ii, i, nd_v)];
                grad[FA3(ii, 1, j, nd_v, 2)] += dx * charge[FA2(ii, i, nd_v)];
                grad[FA3(ii, 2, j, nd_v, 2)] += dy * charge[FA2(ii, i, nd_v)];

                d1 = dipstr[FA2(ii, i, nd_v)] * dipvec[FA3(ii, 1, i, nd_v, 2)];
                d2 = dipstr[FA2(ii, i, nd_v)] * dipvec[FA3(ii, 2, i, nd_v, 2)];

                /* pot = pot + d1*p1 + d2*p2 */
                pot[FA2(ii, j, nd_v)] += d1 * p1;
                pot[FA2(ii, j, nd_v)] += d2 * p2;
                /* grad(1) = grad(1) + d1*dx1 + d2*dx2 */
                grad[FA3(ii, 1, j, nd_v, 2)] += d1 * dx1;
                grad[FA3(ii, 1, j, nd_v, 2)] += d2 * dx2;
                /* grad(2) = grad(2) + d1*dy1 + d2*dy2 */
                grad[FA3(ii, 2, j, nd_v, 2)] += d1 * dy1;
                grad[FA3(ii, 2, j, nd_v, 2)] += d2 * dy2;
            }
        skip:;
        }
    }
}

/*
 * l2d_directcdh: charges + dipoles -> potential + gradient + hessian
 * Threshold: rr .lt. thresh2
 */
void FNAME(l2d_directcdh)(const fint *nd, const double *sources, const fint *ns,
                          const fcomplex *charge, const fcomplex *dipstr,
                          const double *dipvec,
                          const double *targ, const fint *nt,
                          fcomplex *pot, fcomplex *grad, fcomplex *hess,
                          const double *thresh)
{
    fint i, j, ii;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    double xdiff, ydiff, rr, thresh2, rtmp;
    double xdiff2, ydiff2, p1, p2;
    double rr2, rr3, dx, dy, dxx, dxy, dyy;
    double dx1, dx2, dy1, dy2;
    double dxx1, dxx2, dxy1, dxy2, dyy1, dyy2;
    fcomplex d1, d2;

    thresh2 = (*thresh) * (*thresh);

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1, j, 2)] - sources[FA2(1, i, 2)];
            ydiff = targ[FA2(2, j, 2)] - sources[FA2(2, i, 2)];
            xdiff2 = xdiff * xdiff;
            ydiff2 = ydiff * ydiff;
            rr = xdiff2 + ydiff2;
            if (rr < thresh2) goto skip;
            rtmp = log(rr) / 2;

            rr2 = rr * rr;
            rr3 = rr2 * rr;

            dx = xdiff / rr;
            dy = ydiff / rr;
            dxx = (rr - 2 * xdiff2) / rr2;
            dxy = -2 * xdiff * ydiff / rr2;
            dyy = (rr - 2 * ydiff2) / rr2;

            p1 = -dx;
            p2 = -dy;

            dx1 = -dxx;
            dx2 = -dxy;
            dy1 = dx2;
            dy2 = -dyy;

            dxx1 = -2 * xdiff * (xdiff2 - 3 * ydiff2) / rr3;
            dxx2 = 2 * ydiff * (ydiff2 - 3 * xdiff2) / rr3;
            dxy1 = dxx2;
            dxy2 = -dxx1;
            dyy1 = dxy2;
            dyy2 = -dxx2;

            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii, j, nd_v)] += rtmp * charge[FA2(ii, i, nd_v)];
                grad[FA3(ii, 1, j, nd_v, 2)] += dx * charge[FA2(ii, i, nd_v)];
                grad[FA3(ii, 2, j, nd_v, 2)] += dy * charge[FA2(ii, i, nd_v)];
                hess[FA3(ii, 1, j, nd_v, 3)] += dxx * charge[FA2(ii, i, nd_v)];
                hess[FA3(ii, 2, j, nd_v, 3)] += dxy * charge[FA2(ii, i, nd_v)];
                hess[FA3(ii, 3, j, nd_v, 3)] += dyy * charge[FA2(ii, i, nd_v)];

                d1 = dipstr[FA2(ii, i, nd_v)] * dipvec[FA3(ii, 1, i, nd_v, 2)];
                d2 = dipstr[FA2(ii, i, nd_v)] * dipvec[FA3(ii, 2, i, nd_v, 2)];

                /* pot = pot + d1*p1 + d2*p2 */
                pot[FA2(ii, j, nd_v)] += d1 * p1;
                pot[FA2(ii, j, nd_v)] += d2 * p2;
                /* grad = grad + d1*dk1 + d2*dk2 */
                grad[FA3(ii, 1, j, nd_v, 2)] += d1 * dx1;
                grad[FA3(ii, 1, j, nd_v, 2)] += d2 * dx2;
                grad[FA3(ii, 2, j, nd_v, 2)] += d1 * dy1;
                grad[FA3(ii, 2, j, nd_v, 2)] += d2 * dy2;
                /* hess = hess + d1*dkk1 + d2*dkk2 */
                hess[FA3(ii, 1, j, nd_v, 3)] += d1 * dxx1;
                hess[FA3(ii, 1, j, nd_v, 3)] += d2 * dxx2;
                hess[FA3(ii, 2, j, nd_v, 3)] += d1 * dxy1;
                hess[FA3(ii, 2, j, nd_v, 3)] += d2 * dxy2;
                hess[FA3(ii, 3, j, nd_v, 3)] += d1 * dyy1;
                hess[FA3(ii, 3, j, nd_v, 3)] += d2 * dyy2;
            }
        skip:;
        }
    }
}
