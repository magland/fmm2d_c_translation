/*
 * helmkernels2d.c - C translation of src/helmholtz/helmkernels2d.f
 *
 * Direct evaluation kernels for the 2D Helmholtz FMM. Each routine
 * INCREMENTS its output arrays (pot, grad, hess). Uses H_0(k*r)*(i/4)
 * as the Green's function.
 *
 * Strict 1:1 translation. Multi-term Fortran additions are split into
 * separate += statements to preserve left-to-right evaluation order.
 */

#include "helmkernels2d.h"

/* Cross-file dependency */
extern void hank103_(const fcomplex *z, fcomplex *h0, fcomplex *h1,
                     const fint *ifexpon);

/* ================================================================
 * h2d_directcp: charges -> potential
 * ================================================================ */
void FNAME(h2d_directcp)(const fint *nd_p, const fcomplex *wavek_p,
    const double *sources, const fint *ns_p, const fcomplex *charge,
    const double *targ, const fint *nt_p, fcomplex *pot,
    const double *thresh_p)
{
    fint nd_v = *nd_p, ns_v = *ns_p, nt_v = *nt_p;
    fcomplex wavek_v = *wavek_p;
    double thresh_v = *thresh_p;
    fcomplex ima4inv = I / 4.0;
    fint i, j, ii, ifexpon;
    double xdiff, ydiff, rr, r;
    fcomplex z, h0, h1;

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1,j,2)] - sources[FA2(1,i,2)];
            ydiff = targ[FA2(2,j,2)] - sources[FA2(2,i,2)];
            rr = xdiff * xdiff + ydiff * ydiff;
            r = sqrt(rr);
            z = wavek_v * r;
            if (cabs(z) <= thresh_v) goto skip_cp;
            ifexpon = 1;
            hank103_(&z, &h0, &h1, &ifexpon);
            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii,j,nd_v)] += h0 * charge[FA2(ii,i,nd_v)] * ima4inv;
            }
        skip_cp:;
        }
    }
}

/* ================================================================
 * h2d_directcg: charges -> potential + gradient
 * ================================================================ */
void FNAME(h2d_directcg)(const fint *nd_p, const fcomplex *wavek_p,
    const double *sources, const fint *ns_p, const fcomplex *charge,
    const double *targ, const fint *nt_p, fcomplex *pot,
    fcomplex *grad, const double *thresh_p)
{
    fint nd_v = *nd_p, ns_v = *ns_p, nt_v = *nt_p;
    fcomplex wavek_v = *wavek_p;
    double thresh_v = *thresh_p;
    fcomplex ima4inv = I / 4.0;
    fint i, j, ii, ifexpon;
    double xdiff, ydiff, rr, r;
    fcomplex z, h0, h1, cd, cdd;

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1,j,2)] - sources[FA2(1,i,2)];
            ydiff = targ[FA2(2,j,2)] - sources[FA2(2,i,2)];
            rr = xdiff * xdiff + ydiff * ydiff;
            r = sqrt(rr);
            z = wavek_v * r;
            if (cabs(z) < thresh_v) goto skip_cg;
            ifexpon = 1;
            hank103_(&z, &h0, &h1, &ifexpon);
            cdd = -h1 * (wavek_v * ima4inv / r);
            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii,j,nd_v)] += h0 * charge[FA2(ii,i,nd_v)] * ima4inv;
                cd = cdd * charge[FA2(ii,i,nd_v)];
                grad[FA3(ii,1,j,nd_v,2)] += cd * xdiff;
                grad[FA3(ii,2,j,nd_v,2)] += cd * ydiff;
            }
        skip_cg:;
        }
    }
}

/* ================================================================
 * h2d_directch: charges -> potential + gradient + hessian
 * ================================================================ */
void FNAME(h2d_directch)(const fint *nd_p, const fcomplex *wavek_p,
    const double *sources, const fint *ns_p, const fcomplex *charge,
    const double *targ, const fint *nt_p, fcomplex *pot,
    fcomplex *grad, fcomplex *hess, const double *thresh_p)
{
    fint nd_v = *nd_p, ns_v = *ns_p, nt_v = *nt_p;
    fcomplex wavek_v = *wavek_p;
    double thresh_v = *thresh_p;
    fcomplex ima4inv = I / 4.0;
    fint i, j, ii, ifexpon;
    double xdiff, ydiff, rr, r;
    fcomplex z, h0, h1, cd, cdd, cdd2, h2z, hf1, hf2, hf3;

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1,j,2)] - sources[FA2(1,i,2)];
            ydiff = targ[FA2(2,j,2)] - sources[FA2(2,i,2)];
            rr = xdiff * xdiff + ydiff * ydiff;
            r = sqrt(rr);
            z = wavek_v * r;
            if (cabs(z) < thresh_v) goto skip_ch;
            ifexpon = 1;
            hank103_(&z, &h0, &h1, &ifexpon);
            cdd = -h1 * (wavek_v * ima4inv / r);
            cdd2 = (wavek_v * ima4inv / r) / rr;
            h2z = (-z * h0 + 2 * h1);
            hf1 = (h2z * xdiff * xdiff - rr * h1);
            hf2 = (h2z * xdiff * ydiff);
            hf3 = (h2z * ydiff * ydiff - rr * h1);
            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii,j,nd_v)] += h0 * charge[FA2(ii,i,nd_v)] * ima4inv;
                cd = cdd * charge[FA2(ii,i,nd_v)];
                grad[FA3(ii,1,j,nd_v,2)] += cd * xdiff;
                grad[FA3(ii,2,j,nd_v,2)] += cd * ydiff;
                cd = cdd2 * charge[FA2(ii,i,nd_v)];
                hess[FA3(ii,1,j,nd_v,3)] += cd * hf1;
                hess[FA3(ii,2,j,nd_v,3)] += cd * hf2;
                hess[FA3(ii,3,j,nd_v,3)] += cd * hf3;
            }
        skip_ch:;
        }
    }
}

/* ================================================================
 * h2d_directdp: dipoles -> potential
 * ================================================================ */
void FNAME(h2d_directdp)(const fint *nd_p, const fcomplex *wavek_p,
    const double *sources, const fint *ns_p, const fcomplex *dipstr,
    const double *dipvec, const double *targ, const fint *nt_p,
    fcomplex *pot, const double *thresh_p)
{
    fint nd_v = *nd_p, ns_v = *ns_p, nt_v = *nt_p;
    fcomplex wavek_v = *wavek_p;
    double thresh_v = *thresh_p;
    fcomplex ima4inv = I / 4.0;
    fint i, j, ii, ifexpon;
    double xdiff, ydiff, rr, r, ctheta, stheta, dotprod;
    fcomplex z, h0, h1, h2, h3, cd, cdd;

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1,j,2)] - sources[FA2(1,i,2)];
            ydiff = targ[FA2(2,j,2)] - sources[FA2(2,i,2)];
            rr = xdiff * xdiff + ydiff * ydiff;
            r = sqrt(rr);
            z = wavek_v * r;
            if (cabs(z) < thresh_v) goto skip_dp;
            ifexpon = 1;
            hank103_(&z, &h0, &h1, &ifexpon);
            cdd = h1 / r * wavek_v * ima4inv;
            ctheta = xdiff / r;
            stheta = ydiff / r;
            h2 = 2 * h1 / z - h0;
            h3 = 4 * h2 / z - h1;
            for (ii = 1; ii <= nd_v; ii++) {
                cd = cdd * dipstr[FA2(ii,i,nd_v)];
                dotprod = xdiff * dipvec[FA3(ii,1,i,nd_v,2)]
                        + ydiff * dipvec[FA3(ii,2,i,nd_v,2)];
                pot[FA2(ii,j,nd_v)] += cd * dotprod;
            }
        skip_dp:;
        }
    }
}

/* ================================================================
 * h2d_directdg: dipoles -> potential + gradient
 * ================================================================ */
void FNAME(h2d_directdg)(const fint *nd_p, const fcomplex *wavek_p,
    const double *sources, const fint *ns_p, const fcomplex *dipstr,
    const double *dipvec, const double *targ, const fint *nt_p,
    fcomplex *pot, fcomplex *grad, const double *thresh_p)
{
    fint nd_v = *nd_p, ns_v = *ns_p, nt_v = *nt_p;
    fcomplex wavek_v = *wavek_p;
    double thresh_v = *thresh_p;
    fcomplex ima4inv = I / 4.0;
    fint i, j, ii, ifexpon;
    double xdiff, ydiff, rr, r, ctheta, stheta, dotprod;
    fcomplex z, h0, h1, h2, h3, cd, cdd, cdd2;
    fcomplex hxx, hxy, hyy, hf1, hf2, hf3;

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1,j,2)] - sources[FA2(1,i,2)];
            ydiff = targ[FA2(2,j,2)] - sources[FA2(2,i,2)];
            rr = xdiff * xdiff + ydiff * ydiff;
            r = sqrt(rr);
            z = wavek_v * r;
            if (cabs(z) < thresh_v) goto skip_dg;
            ifexpon = 1;
            hank103_(&z, &h0, &h1, &ifexpon);
            ctheta = xdiff / r;
            stheta = ydiff / r;
            h2 = 2 * h1 / z - h0;
            h3 = 4 * h2 / z - h1;
            cdd = h1 / r * wavek_v * ima4inv;
            cdd2 = -wavek_v * wavek_v * ima4inv;
            hf1 = (h2 * (ctheta * ctheta - 0.5) - h0 / 2);
            hf2 = (h2 * ctheta * stheta);
            hf3 = (h2 * (stheta * stheta - 0.5) - h0 / 2);
            for (ii = 1; ii <= nd_v; ii++) {
                cd = cdd * dipstr[FA2(ii,i,nd_v)];
                dotprod = xdiff * dipvec[FA3(ii,1,i,nd_v,2)]
                        + ydiff * dipvec[FA3(ii,2,i,nd_v,2)];
                pot[FA2(ii,j,nd_v)] += cd * dotprod;
                cd = cdd2 * dipstr[FA2(ii,i,nd_v)];
                hxx = cd * hf1;
                hxy = cd * hf2;
                hyy = cd * hf3;
                /* Split: grad = grad + hxx*dv1 + hxy*dv2 */
                grad[FA3(ii,1,j,nd_v,2)] += hxx * dipvec[FA3(ii,1,i,nd_v,2)];
                grad[FA3(ii,1,j,nd_v,2)] += hxy * dipvec[FA3(ii,2,i,nd_v,2)];
                grad[FA3(ii,2,j,nd_v,2)] += hxy * dipvec[FA3(ii,1,i,nd_v,2)];
                grad[FA3(ii,2,j,nd_v,2)] += hyy * dipvec[FA3(ii,2,i,nd_v,2)];
            }
        skip_dg:;
        }
    }
}

/* ================================================================
 * h2d_directdh: dipoles -> potential + gradient + hessian
 * ================================================================ */
void FNAME(h2d_directdh)(const fint *nd_p, const fcomplex *wavek_p,
    const double *sources, const fint *ns_p, const fcomplex *dipstr,
    const double *dipvec, const double *targ, const fint *nt_p,
    fcomplex *pot, fcomplex *grad, fcomplex *hess,
    const double *thresh_p)
{
    fint nd_v = *nd_p, ns_v = *ns_p, nt_v = *nt_p;
    fcomplex wavek_v = *wavek_p;
    double thresh_v = *thresh_p;
    fcomplex ima4inv = I / 4.0;
    fint i, j, ii, ifexpon;
    double xdiff, ydiff, rr, r, ctheta, stheta, dotprod;
    fcomplex z, h0, h1, h2, h3, cd, cdd, cdd2;
    fcomplex hxx, hxy, hyy, hf1, hf2, hf3;
    fcomplex hxxx, hxxy, hxyy, hyyy;
    fcomplex hxxx1, hxxy1, hxyy1, hyyy1;

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1,j,2)] - sources[FA2(1,i,2)];
            ydiff = targ[FA2(2,j,2)] - sources[FA2(2,i,2)];
            rr = xdiff * xdiff + ydiff * ydiff;
            r = sqrt(rr);
            z = wavek_v * r;
            if (cabs(z) < thresh_v) goto skip_dh;
            ifexpon = 1;
            hank103_(&z, &h0, &h1, &ifexpon);
            ctheta = xdiff / r;
            stheta = ydiff / r;
            h2 = 2 * h1 / z - h0;
            h3 = 4 * h2 / z - h1;
            cdd = h1 / r * wavek_v * ima4inv;
            cdd2 = -wavek_v * wavek_v * ima4inv;
            hf1 = (h2 * (ctheta * ctheta - 0.5) - h0 / 2);
            hf2 = (h2 * ctheta * stheta);
            hf3 = (h2 * (stheta * stheta - 0.5) - h0 / 2);
            hxxx1 = (-h1 / 2 * (-1.5)
                - h3 / 2 * (ctheta * ctheta - 0.5 - stheta * stheta)
                ) * ctheta;
            hxxy1 = (-h1 / 2 * (-0.5)
                - h3 / 2 * (1.5 * ctheta * ctheta - 0.5 * stheta * stheta)
                ) * stheta;
            hxyy1 = (-h1 / 2 * (-0.5)
                - h3 / 2 * (1.5 * stheta * stheta - 0.5 * ctheta * ctheta)
                ) * ctheta;
            hyyy1 = (-h1 / 2 * (-1.5)
                - h3 / 2 * (stheta * stheta - 0.5 - ctheta * ctheta)
                ) * stheta;

            for (ii = 1; ii <= nd_v; ii++) {
                dotprod = xdiff * dipvec[FA3(ii,1,i,nd_v,2)]
                        + ydiff * dipvec[FA3(ii,2,i,nd_v,2)];
                cd = cdd * dipstr[FA2(ii,i,nd_v)];
                pot[FA2(ii,j,nd_v)] += cd * dotprod;
                cd = cdd2 * dipstr[FA2(ii,i,nd_v)];
                hxx = cd * hf1;
                hxy = cd * hf2;
                hyy = cd * hf3;
                /* Split multi-term additions */
                grad[FA3(ii,1,j,nd_v,2)] += hxx * dipvec[FA3(ii,1,i,nd_v,2)];
                grad[FA3(ii,1,j,nd_v,2)] += hxy * dipvec[FA3(ii,2,i,nd_v,2)];
                grad[FA3(ii,2,j,nd_v,2)] += hxy * dipvec[FA3(ii,1,i,nd_v,2)];
                grad[FA3(ii,2,j,nd_v,2)] += hyy * dipvec[FA3(ii,2,i,nd_v,2)];

                cd = -wavek_v * wavek_v * wavek_v * dipstr[FA2(ii,i,nd_v)] * ima4inv;
                hxxx = cd * hxxx1;
                hxxy = cd * hxxy1;
                hxyy = cd * hxyy1;
                hyyy = cd * hyyy1;
                /* Split multi-term additions */
                hess[FA3(ii,1,j,nd_v,3)] += hxxx * dipvec[FA3(ii,1,i,nd_v,2)];
                hess[FA3(ii,1,j,nd_v,3)] += hxxy * dipvec[FA3(ii,2,i,nd_v,2)];
                hess[FA3(ii,2,j,nd_v,3)] += hxxy * dipvec[FA3(ii,1,i,nd_v,2)];
                hess[FA3(ii,2,j,nd_v,3)] += hxyy * dipvec[FA3(ii,2,i,nd_v,2)];
                hess[FA3(ii,3,j,nd_v,3)] += hxyy * dipvec[FA3(ii,1,i,nd_v,2)];
                hess[FA3(ii,3,j,nd_v,3)] += hyyy * dipvec[FA3(ii,2,i,nd_v,2)];
            }
        skip_dh:;
        }
    }
}

/* ================================================================
 * h2d_directcdp: charges + dipoles -> potential
 * ================================================================ */
void FNAME(h2d_directcdp)(const fint *nd_p, const fcomplex *wavek_p,
    const double *sources, const fint *ns_p, const fcomplex *charge,
    const fcomplex *dipstr, const double *dipvec,
    const double *targ, const fint *nt_p, fcomplex *pot,
    const double *thresh_p)
{
    fint nd_v = *nd_p, ns_v = *ns_p, nt_v = *nt_p;
    fcomplex wavek_v = *wavek_p;
    double thresh_v = *thresh_p;
    fcomplex ima4inv = I / 4.0;
    fint i, j, ii, ifexpon;
    double xdiff, ydiff, rr, r, ctheta, stheta, dotprod;
    fcomplex z, h0, h1, h2, h3, cd, cdd;

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1,j,2)] - sources[FA2(1,i,2)];
            ydiff = targ[FA2(2,j,2)] - sources[FA2(2,i,2)];
            rr = xdiff * xdiff + ydiff * ydiff;
            r = sqrt(rr);
            z = wavek_v * r;
            if (cabs(z) < thresh_v) goto skip_cdp;
            ifexpon = 1;
            hank103_(&z, &h0, &h1, &ifexpon);
            ctheta = xdiff / r;
            stheta = ydiff / r;
            h2 = 2 * h1 / z - h0;
            h3 = 4 * h2 / z - h1;
            cdd = h1 / r * wavek_v * ima4inv;
            for (ii = 1; ii <= nd_v; ii++) {
                cd = cdd * dipstr[FA2(ii,i,nd_v)];
                pot[FA2(ii,j,nd_v)] += h0 * charge[FA2(ii,i,nd_v)] * ima4inv;
                dotprod = xdiff * dipvec[FA3(ii,1,i,nd_v,2)]
                        + ydiff * dipvec[FA3(ii,2,i,nd_v,2)];
                pot[FA2(ii,j,nd_v)] += cd * dotprod;
            }
        skip_cdp:;
        }
    }
}

/* ================================================================
 * h2d_directcdg: charges + dipoles -> potential + gradient
 * ================================================================ */
void FNAME(h2d_directcdg)(const fint *nd_p, const fcomplex *wavek_p,
    const double *sources, const fint *ns_p, const fcomplex *charge,
    const fcomplex *dipstr, const double *dipvec,
    const double *targ, const fint *nt_p, fcomplex *pot,
    fcomplex *grad, const double *thresh_p)
{
    fint nd_v = *nd_p, ns_v = *ns_p, nt_v = *nt_p;
    fcomplex wavek_v = *wavek_p;
    double thresh_v = *thresh_p;
    fcomplex ima4inv = I / 4.0;
    fint i, j, ii, ifexpon;
    double xdiff, ydiff, rr, r, ctheta, stheta, dotprod;
    fcomplex z, h0, h1, h2, h3, cd, cdd, cdd2;
    fcomplex hxx, hxy, hyy, hf1, hf2, hf3;

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1,j,2)] - sources[FA2(1,i,2)];
            ydiff = targ[FA2(2,j,2)] - sources[FA2(2,i,2)];
            rr = xdiff * xdiff + ydiff * ydiff;
            r = sqrt(rr);
            z = wavek_v * r;
            if (cabs(z) < thresh_v) goto skip_cdg;
            ifexpon = 1;
            hank103_(&z, &h0, &h1, &ifexpon);
            cdd = -h1 * (wavek_v * ima4inv / r);
            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii,j,nd_v)] += h0 * charge[FA2(ii,i,nd_v)] * ima4inv;
                cd = cdd * charge[FA2(ii,i,nd_v)];
                grad[FA3(ii,1,j,nd_v,2)] += cd * xdiff;
                grad[FA3(ii,2,j,nd_v,2)] += cd * ydiff;
            }
            ctheta = xdiff / r;
            stheta = ydiff / r;
            h2 = 2.0 * h1 / z - h0;
            h3 = 4.0 * h2 / z - h1;
            cdd = h1 / r * wavek_v * ima4inv;
            cdd2 = -wavek_v * wavek_v * ima4inv;
            hf1 = (h2 * (ctheta * ctheta - 0.5) - h0 / 2);
            hf2 = (h2 * ctheta * stheta);
            hf3 = (h2 * (stheta * stheta - 0.5) - h0 / 2);
            for (ii = 1; ii <= nd_v; ii++) {
                cd = cdd * dipstr[FA2(ii,i,nd_v)];
                dotprod = xdiff * dipvec[FA3(ii,1,i,nd_v,2)]
                        + ydiff * dipvec[FA3(ii,2,i,nd_v,2)];
                pot[FA2(ii,j,nd_v)] += cd * dotprod;
                cd = cdd2 * dipstr[FA2(ii,i,nd_v)];
                hxx = cd * hf1;
                hxy = cd * hf2;
                hyy = cd * hf3;
                /* Split multi-term additions */
                grad[FA3(ii,1,j,nd_v,2)] += hxx * dipvec[FA3(ii,1,i,nd_v,2)];
                grad[FA3(ii,1,j,nd_v,2)] += hxy * dipvec[FA3(ii,2,i,nd_v,2)];
                grad[FA3(ii,2,j,nd_v,2)] += hxy * dipvec[FA3(ii,1,i,nd_v,2)];
                grad[FA3(ii,2,j,nd_v,2)] += hyy * dipvec[FA3(ii,2,i,nd_v,2)];
            }
        skip_cdg:;
        }
    }
}

/* ================================================================
 * h2d_directcdh: charges + dipoles -> potential + gradient + hessian
 * ================================================================ */
void FNAME(h2d_directcdh)(const fint *nd_p, const fcomplex *wavek_p,
    const double *sources, const fint *ns_p, const fcomplex *charge,
    const fcomplex *dipstr, const double *dipvec,
    const double *targ, const fint *nt_p, fcomplex *pot,
    fcomplex *grad, fcomplex *hess, const double *thresh_p)
{
    fint nd_v = *nd_p, ns_v = *ns_p, nt_v = *nt_p;
    fcomplex wavek_v = *wavek_p;
    double thresh_v = *thresh_p;
    fcomplex ima4inv = I / 4.0;
    fint i, j, ii, ifexpon;
    double xdiff, ydiff, rr, r, ctheta, stheta, dotprod;
    fcomplex z, h0, h1, h2, h3, cd, cdd, cdd2, h2z;
    fcomplex hxx, hxy, hyy, hf1, hf2, hf3;
    fcomplex hxxx, hxxy, hxyy, hyyy;
    fcomplex hxxx1, hxxy1, hxyy1, hyyy1;

    for (j = 1; j <= nt_v; j++) {
        for (i = 1; i <= ns_v; i++) {
            xdiff = targ[FA2(1,j,2)] - sources[FA2(1,i,2)];
            ydiff = targ[FA2(2,j,2)] - sources[FA2(2,i,2)];
            rr = xdiff * xdiff + ydiff * ydiff;
            r = sqrt(rr);
            z = wavek_v * r;
            if (cabs(z) < thresh_v) goto skip_cdh;
            ifexpon = 1;
            hank103_(&z, &h0, &h1, &ifexpon);
            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii,j,nd_v)] += h0 * charge[FA2(ii,i,nd_v)] * ima4inv;
                cd = -h1 * (wavek_v * charge[FA2(ii,i,nd_v)] * ima4inv / r);
                grad[FA3(ii,1,j,nd_v,2)] += cd * xdiff;
                grad[FA3(ii,2,j,nd_v,2)] += cd * ydiff;
                cd = (wavek_v * charge[FA2(ii,i,nd_v)] * ima4inv / r) / rr;
                h2z = (-z * h0 + 2 * h1);
                hess[FA3(ii,1,j,nd_v,3)] += cd * (h2z * xdiff * xdiff - rr * h1);
                hess[FA3(ii,2,j,nd_v,3)] += cd * (h2z * xdiff * ydiff);
                hess[FA3(ii,3,j,nd_v,3)] += cd * (h2z * ydiff * ydiff - rr * h1);
            }
            ctheta = xdiff / r;
            stheta = ydiff / r;
            h2 = 2 * h1 / z - h0;
            h3 = 4 * h2 / z - h1;
            cdd = h1 / r * wavek_v * ima4inv;
            cdd2 = -wavek_v * wavek_v * ima4inv;
            hf1 = (h2 * (ctheta * ctheta - 0.5) - h0 / 2);
            hf2 = (h2 * ctheta * stheta);
            hf3 = (h2 * (stheta * stheta - 0.5) - h0 / 2);
            hxxx1 = (-h1 / 2 * (-1.5)
                - h3 / 2 * (ctheta * ctheta - 0.5 - stheta * stheta)
                ) * ctheta;
            hxxy1 = (-h1 / 2 * (-0.5)
                - h3 / 2 * (1.5 * ctheta * ctheta - 0.5 * stheta * stheta)
                ) * stheta;
            hxyy1 = (-h1 / 2 * (-0.5)
                - h3 / 2 * (1.5 * stheta * stheta - 0.5 * ctheta * ctheta)
                ) * ctheta;
            hyyy1 = (-h1 / 2 * (-1.5)
                - h3 / 2 * (stheta * stheta - 0.5 - ctheta * ctheta)
                ) * stheta;
            for (ii = 1; ii <= nd_v; ii++) {
                cd = cdd * dipstr[FA2(ii,i,nd_v)];
                dotprod = xdiff * dipvec[FA3(ii,1,i,nd_v,2)]
                        + ydiff * dipvec[FA3(ii,2,i,nd_v,2)];
                pot[FA2(ii,j,nd_v)] += cd * dotprod;
                cd = cdd2 * dipstr[FA2(ii,i,nd_v)];
                hxx = cd * hf1;
                hxy = cd * hf2;
                hyy = cd * hf3;
                /* Split multi-term additions */
                grad[FA3(ii,1,j,nd_v,2)] += hxx * dipvec[FA3(ii,1,i,nd_v,2)];
                grad[FA3(ii,1,j,nd_v,2)] += hxy * dipvec[FA3(ii,2,i,nd_v,2)];
                grad[FA3(ii,2,j,nd_v,2)] += hxy * dipvec[FA3(ii,1,i,nd_v,2)];
                grad[FA3(ii,2,j,nd_v,2)] += hyy * dipvec[FA3(ii,2,i,nd_v,2)];
                cd = -wavek_v * wavek_v * wavek_v * dipstr[FA2(ii,i,nd_v)] * ima4inv;
                hxxx = cd * hxxx1;
                hxxy = cd * hxxy1;
                hxyy = cd * hxyy1;
                hyyy = cd * hyyy1;
                /* Split multi-term additions */
                hess[FA3(ii,1,j,nd_v,3)] += hxxx * dipvec[FA3(ii,1,i,nd_v,2)];
                hess[FA3(ii,1,j,nd_v,3)] += hxxy * dipvec[FA3(ii,2,i,nd_v,2)];
                hess[FA3(ii,2,j,nd_v,3)] += hxxy * dipvec[FA3(ii,1,i,nd_v,2)];
                hess[FA3(ii,2,j,nd_v,3)] += hxyy * dipvec[FA3(ii,2,i,nd_v,2)];
                hess[FA3(ii,3,j,nd_v,3)] += hxyy * dipvec[FA3(ii,1,i,nd_v,2)];
                hess[FA3(ii,3,j,nd_v,3)] += hyyy * dipvec[FA3(ii,2,i,nd_v,2)];
            }
        skip_cdh:;
        }
    }
}
