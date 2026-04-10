/*
 * cdjseval2d.c - C translation of src/common/cdjseval2d.f
 *
 * Bessel function evaluation via downward recurrence.
 *
 * Strict 1:1 translation. The routine uses implicit real*8 (a-h,o-z).
 */

#include "cdjseval2d.h"

void FNAME(jbessel2d)(const fint *nterms_p, const fcomplex *z_p,
                      const double *rscale_p, fcomplex *fjs,
                      const fint *ifder_p, fcomplex *fjder)
{
    fint nterms_v = *nterms_p;
    fcomplex z_v = *z_p;
    double rscale_v = *rscale_p;
    fint ifder_v = *ifder_p;

    double upbound = 1.0e+32;
    double upbound2 = 1.0e+40;
    double upbound2inv = 1.0e-40;
    double tiny_val = 1.0e-200;
    double done = 1.0;
    double zero = 0.0;
    fcomplex ima = I;

    fint ier, i;
    int ntop, nextra, ncntr;
    double dcoef, dd, scalinv, sctot, dc1, dc2;
    fcomplex zinv, ztmp, ffn, ffnm1, ffnp1;
    fcomplex psi, zmul, zsn, zmulinv, zscale;

    int *iscale;
    fcomplex *fjtemp;

    ier = 0;

    /* set to asymptotic values if argument is sufficiently small */
    if (cabs(z_v) < tiny_val) {
        fjs[0] = done;
        for (i = 1; i <= nterms_v; i++) {
            fjs[i] = zero;
        }

        if (ifder_v == 1) {
            for (i = 0; i <= nterms_v; i++) {
                fjder[i] = zero;
            }
            fjder[1] = done / (2 * rscale_v);
        }

        return;
    }

    /* Step 1: carry out upward recursion starting at nterms until the
       magnitude has grown by upbound2 = 10^40 */
    ntop = nterms_v + 10000;  /* initial estimate */
    zinv = done / z_v;
    ffn = done;
    ffnm1 = zero;

    nextra = 10000;
    for (i = nterms_v; i <= nterms_v + nextra; i++) {
        dcoef = 2 * i;
        ztmp = dcoef * zinv * ffn - ffnm1;
        ffnp1 = ztmp;

        dd = creal(ztmp) * creal(ztmp) + cimag(ztmp) * cimag(ztmp);
        if (dd > upbound2) {
            ntop = i + 1;
            goto label_1300;
        }
        ffnm1 = ffn;
        ffn = ffnp1;
    }
label_1300:

    iscale = (int *)malloc((ntop + 1) * sizeof(int));
    fjtemp = (fcomplex *)malloc((ntop + 1) * sizeof(fcomplex));

    /* Step 2: Recursion back down to generate the unscaled jfuns */
    for (i = 0; i <= ntop; i++) {
        iscale[i] = 0;
    }

    fjtemp[ntop] = zero;
    fjtemp[ntop - 1] = done;
    for (i = ntop - 1; i >= 1; i--) {
        dcoef = 2 * i;
        ztmp = dcoef * zinv * fjtemp[i] - fjtemp[i + 1];
        fjtemp[i - 1] = ztmp;

        dd = creal(ztmp) * creal(ztmp) + cimag(ztmp) * cimag(ztmp);
        if (dd > upbound2) {
            fjtemp[i] = fjtemp[i] * upbound2inv;
            fjtemp[i - 1] = fjtemp[i - 1] * upbound2inv;
            iscale[i] = 1;
        }
    }

    /* Step 3: go back up and scale uniformly */
    ncntr = 0;
    scalinv = done / rscale_v;
    sctot = 1.0;
    for (i = 1; i <= ntop; i++) {
        sctot = sctot * scalinv;
        if (iscale[i - 1] == 1) sctot = sctot * upbound2inv;
        fjtemp[i] = fjtemp[i] * sctot;
    }

    /* Determine the normalization parameter */
    psi = 0;

    if (cimag(z_v) < 0) zmul = +ima;
    if (cimag(z_v) >= 0) zmul = -ima;

    /* zmul**(mod(ntop,4)) */
    {
        int r = ntop % 4;
        zsn = 1.0;
        for (i = 0; i < r; i++) {
            zsn = zsn * zmul;
        }
    }

    zmulinv = 1.0 / zmul;
    for (i = ntop; i >= 1; i--) {
        psi = rscale_v * psi + fjtemp[i] * zsn;
        zsn = zsn * zmulinv;
    }
    psi = 2 * psi * rscale_v + fjtemp[0];

    if (cimag(z_v) < 0) zscale = cexp(+ima * z_v) / psi;
    if (cimag(z_v) >= 0) zscale = cexp(-ima * z_v) / psi;

    /* Scale the jfuns by zscale */
    ztmp = zscale;
    for (i = 0; i <= nterms_v; i++) {
        fjs[i] = fjtemp[i] * ztmp;
    }

    /* Calculate the derivatives if desired */
    if (ifder_v == 1) {
        fjs[nterms_v + 1] = fjtemp[nterms_v + 1] * ztmp;

        fjder[0] = -fjs[1] * rscale_v;
        dc1 = 0.5;
        dc2 = done - dc1;
        dc1 = dc1 * scalinv;
        dc2 = dc2 * rscale_v;
        for (i = 1; i <= nterms_v; i++) {
            fjder[i] = dc1 * fjs[i - 1] - dc2 * fjs[i + 1];
        }
    }

    free(iscale);
    free(fjtemp);
    return;
}
