/*
 * h2dcommon.c - C translation of src/helmholtz/h2dcommon.f
 *
 * Common routines for Helmholtz 2D FMM:
 *   h2cart2polar - Cartesian to polar conversion
 *   h2dall      - scaled Hankel functions H_0..H_n and optional derivatives
 *   h2dmpzero   - zero out a multipole expansion
 *   h2dsigzero  - zero out a signature expansion
 *
 * Strict 1:1 translation.
 */

#include "h2dcommon.h"

/* Cross-file dependency: hank103 from hank103.f */
extern void hank103_(const fcomplex *z, fcomplex *h0, fcomplex *h1,
                     const fint *ifexpon);

/* ================================================================
 * h2cart2polar - convert from Cartesian to polar coordinates
 * ================================================================ */
void FNAME(h2cart2polar)(const double *zat, double *r, double *theta)
{
    *r = sqrt(zat[0] * zat[0] + zat[1] * zat[1]);
    if (fabs(zat[0]) == 0 && fabs(zat[1]) == 0) {
        *theta = 0;
    } else {
        *theta = atan2(zat[1], zat[0]);
    }
}

/* ================================================================
 * h2dall - scaled Hankel functions H_n(z)*rscale^n
 * ================================================================ */
void FNAME(h2dall)(const fint *nterms_p, const fcomplex *z_p,
                   const double *rscale_p, fcomplex *hvec,
                   const fint *ifder_p, fcomplex *hder)
{
    fint nterms_v = *nterms_p;
    fcomplex z_v = *z_p;
    double rscale_v = *rscale_p;
    fint ifder_v = *ifder_p;

    double thresh = 1.0e-200;
    double done = 1.0;
    fint i, ifexpon;
    double scal2, dtmp;
    fcomplex zinv, ztmp, h0, h1;

    /* If |z| < thresh, return zeros */
    if (cabs(z_v) < thresh) {
        for (i = 0; i <= nterms_v; i++) {
            hvec[i] = 0;
            hder[i] = 0;
        }
        return;
    }

    /* Get H_0 and H_1 via hank103 */
    ifexpon = 1;
    hank103_(&z_v, &h0, &h1, &ifexpon);
    hvec[0] = h0;
    hvec[1] = h1 * rscale_v;

    /* Recursion: H_{n+1}(z) = rscale*(2*n/z*H_n(z) - H_{n-1}(z)*rscale) */
    scal2 = rscale_v * rscale_v;
    zinv = rscale_v / z_v;
    for (i = 1; i <= nterms_v - 1; i++) {
        dtmp = 2 * i;
        ztmp = zinv * dtmp;
        hvec[i + 1] = ztmp * hvec[i] - scal2 * hvec[i - 1];
    }

    /* Derivatives: H_n'(z) = scale*H_{n-1}(z) - n/z * H_n(z) */
    if (ifder_v == 1) {
        hder[0] = -hvec[1] / rscale_v;
        zinv = 1.0 / z_v;
        for (i = 1; i <= nterms_v; i++) {
            dtmp = (double)i;
            ztmp = zinv * dtmp;
            hder[i] = rscale_v * hvec[i - 1] - ztmp * hvec[i];
        }
    }

    return;
}

/* ================================================================
 * h2dmpzero - zero out a multipole expansion
 *   mpole(nd, -nterms:nterms)
 * ================================================================ */
void FNAME(h2dmpzero)(const fint *nd_p, fcomplex *mpole, const fint *nterms_p)
{
    fint nd_v = *nd_p;
    fint nterms_v = *nterms_p;
    fint n, idim;

    /* mpole is dimensioned (nd, -nterms:nterms) in Fortran,
       which is (nd, 2*nterms+1) contiguous complex values.
       Index: mpole[idim-1 + (n+nterms)*nd] for idim=1..nd, n=-nterms..nterms */
    for (n = -nterms_v; n <= nterms_v; n++) {
        for (idim = 1; idim <= nd_v; idim++) {
            mpole[(idim - 1) + (n + nterms_v) * nd_v] = 0.0;
        }
    }
}

/* ================================================================
 * h2dsigzero - zero out a signature expansion
 *   sig(nd, nsig)
 * ================================================================ */
void FNAME(h2dsigzero)(const fint *nd_p, fcomplex *sig, const fint *nsig_p)
{
    fint nd_v = *nd_p;
    fint nsig_v = *nsig_p;
    fint n, idim;

    for (n = 1; n <= nsig_v; n++) {
        for (idim = 1; idim <= nd_v; idim++) {
            sig[FA2(idim, n, nd_v)] = 0.0;
        }
    }
}
