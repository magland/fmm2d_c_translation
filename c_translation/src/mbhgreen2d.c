/*
 * mbhgreen2d.c
 *
 * C translation of src/modified-biharmonic/mbhgreen2d.f
 * Green's function and difference kernel evaluations for the
 * modified biharmonic equation.
 *
 * 9 routines:
 *   modbhgreen_all  - Green's function with derivatives up to 5th order
 *   modbhgreen      - convenience wrapper (pot/grad/hess only)
 *   modbhgreend1    - one source directional derivative
 *   modbhgreend2    - two source directional derivatives
 *   difflogbk       - K_0 difference kernel (power series)
 *   diffslogbk      - scaled difference kernel modes
 *   diffslogbk_fast - fast version with precomputed coefficients
 *   diffszkik       - I_0 difference kernel modes
 *   diffszkik_fast  - fast version
 */

#include "fmm2d_c.h"
#include "mbhgreen2d.h"

/* Cross-file calls (already translated in cdjseval2d.c) */
extern void h2dall_(const fint *nterms, const fcomplex *z,
                    const double *rscale, fcomplex *hvec,
                    const fint *ifder, fcomplex *hder);
extern void jbessel2d_(const fint *nterms, const fcomplex *z,
                       const double *rscale, fcomplex *fjs,
                       const fint *ifder, double *fjder);

/* ================================================================
 * modbhgreen_all
 * ================================================================ */
void FNAME(modbhgreen_all)(
    const double *beta_p, const double *zx, const double *zy,
    const fint *ifpot_p, double *pot,
    const fint *ifgrad_p, double *grad,
    const fint *ifhess_p, double *hess,
    const fint *ifder3_p, double *der3,
    const fint *ifder4_p, double *der4,
    const fint *ifder5_p, double *der5)
{
    double beta = *beta_p;
    fint ifpot = *ifpot_p, ifgrad = *ifgrad_p, ifhess = *ifhess_p;
    fint ifder3 = *ifder3_p, ifder4 = *ifder4_p, ifder5 = *ifder5_p;

    double dx = zx[0] - zy[0];
    double dy = zx[1] - zy[1];

    double dx2 = dx * dx;
    double dx3 = dx2 * dx;
    double dx4 = dx3 * dx;
    double dx5 = dx4 * dx;
    double dy2 = dy * dy;
    double dy3 = dy2 * dy;
    double dy4 = dy3 * dy;
    double dy5 = dy4 * dy;

    double r2 = dx2 + dy2;
    double r = sqrt(r2);
    double r3 = r2 * r;
    double r4 = r3 * r;
    double r5 = r4 * r;

    fint nmax = 5;

    /* get values of difference kernel and higher modes */
    double diffs[6]; /* 0:5 */
    double ders_loc[6]; /* 0:5, unused */
    double kvec[11]; /* 0:10 */
    double rscale = 1.0;
    fint ifders = 0;

    FNAME(diffslogbk_fast)(&r, &beta, &rscale, diffs, &ifders,
                           ders_loc, kvec, &nmax);

    /* combine to get values of derivatives of G with respect to r */
    double pih = 2.0 * atan(1.0);
    double dtemp = 1.0 / (4.0 * pih * beta * beta);
    double betah = 0.5 * beta;

    double g0 = diffs[0] * dtemp;
    dtemp = -dtemp * beta;
    double g1 = diffs[1] * dtemp;
    dtemp = -dtemp * betah;
    double g2 = (diffs[2] + kvec[0]) * dtemp;
    dtemp = -dtemp * betah;
    double g3 = (diffs[3] + 3.0 * kvec[1]) * dtemp;
    dtemp = -dtemp * betah;
    double g4 = (diffs[4] + 3.0 * kvec[0] + 4.0 * kvec[2]) * dtemp;
    dtemp = -dtemp * betah;
    double g5 = (diffs[5] + 10.0 * kvec[1] + 5.0 * kvec[3]) * dtemp;

    /* evaluate potential and derivatives */
    if (ifpot == 1) {
        *pot = g0;
    }

    if (ifgrad == 1) {
        grad[0] = dx * g1 / r;
        grad[1] = dy * g1 / r;
    }

    if (ifhess == 1) {
        hess[0] = dx2 * g2 / r2 + g1 * (1.0 / r - dx2 / r3);
        hess[1] = dx * dy * (g2 / r2 - g1 / r3);
        hess[2] = dy2 * g2 / r2 + g1 * (1.0 / r - dy2 / r3);
    }

    if (ifder3 == 1) {
        der3[0] = (dx3 * g3 + 3.0 * dy2 * dx * (g2 / r - g1 / r2)) / r3;
        der3[1] = dx2 * dy * (g3 / r3 - 3.0 * (g2 / r4 - g1 / r5)) +
                  dy * (g2 / r2 - g1 / r3);
        der3[2] = dx * dy2 * (g3 / r3 - 3.0 * (g2 / r4 - g1 / r5)) +
                  dx * (g2 / r2 - g1 / r3);
        der3[3] = (dy3 * g3 + 3.0 * dx2 * dy * (g2 / r - g1 / r2)) / r3;
    }

    if (ifder4 == 1) {
        der4[0] = (dx4 * (g4 - 6.0 * g3 / r + 15.0 * (g2 / r2 - g1 / r3))) / r4 +
                  (dx2 * 6.0 * (g3 - 3.0 * (g2 / r - g1 / r2))) / r3 +
                  (3.0 * (g2 - g1 / r)) / r2;
        der4[1] = (dx3 * dy * (g4 - 6.0 * g3 / r + 15.0 * (g2 / r2 - g1 / r3))) / r4 +
                  (3.0 * dx * dy * (g3 - 3.0 * (g2 / r - g1 / r2))) / r3;
        der4[2] = dx2 * dy2 * (g4 - 6.0 * g3 / r + 15.0 * g2 / r2 - 15.0 * g1 / r3) / r4
                  + g3 / r - 2.0 * g2 / r2 + 2.0 * g1 / r3;
        der4[3] = dx * dy3 * (g4 - 6.0 * g3 / r + 15.0 * (g2 / r2 - g1 / r3)) / r4 +
                  3.0 * dx * dy * (g3 - 3.0 * (g2 / r - g1 / r2)) / r3;
        der4[4] = dy4 * (g4 - 6.0 * g3 / r + 15.0 * (g2 / r2 - g1 / r3)) / r4 +
                  dy2 * 6.0 * (g3 - 3.0 * (g2 / r - g1 / r2)) / r3 +
                  3.0 * (g2 - g1 / r) / r2;
    }

    if (ifder5 == 1) {
        der5[0] = (dx5 * g5 + 10.0 * dy2 * dx3 * g4 / r +
                   (15.0 * dy4 * dx - 30.0 * dy2 * dx3) * g3 / r2 +
                   (60.0 * dy2 * dx3 - 45.0 * dy4 * dx) * g2 / r3 +
                   (45.0 * dy4 * dx - 60.0 * dy2 * dx3) * g1 / r4) / r5;
        der5[1] = (dy * dx4 * g5 + (6.0 * dy3 * dx2 - 4.0 * dy * dx4) * g4 / r +
                   (3.0 * dy5 + 12.0 * dy * dx4 - 30.0 * dy3 * dx2) * g3 / r2 +
                   (72.0 * dy3 * dx2 - 9.0 * dy5 - 24.0 * dy * dx4) * g2 / r3 +
                   (9.0 * dy5 - 72.0 * dy3 * dx2 + 24.0 * dy * dx4) * g1 / r4) / r5;
        der5[2] = (dy2 * dx3 * g5 + (3.0 * dy4 * dx - 6.0 * dy2 * dx3 + dx5) * g4 / r +
                   (27.0 * dy2 * dx3 - 15.0 * dy4 * dx - 3.0 * dx5) * g3 / r2 +
                   (36.0 * dy4 * dx - 63.0 * dy2 * dx3 + 6.0 * dx5) * g2 / r3 +
                   (63.0 * dy2 * dx3 - 36.0 * dy4 * dx - 6.0 * dx5) * g1 / r4) / r5;
        der5[3] = (dx2 * dy3 * g5 + (3.0 * dx4 * dy - 6.0 * dx2 * dy3 + dy5) * g4 / r +
                   (27.0 * dx2 * dy3 - 15.0 * dx4 * dy - 3.0 * dy5) * g3 / r2 +
                   (36.0 * dx4 * dy - 63.0 * dx2 * dy3 + 6.0 * dy5) * g2 / r3 +
                   (63.0 * dx2 * dy3 - 36.0 * dx4 * dy - 6.0 * dy5) * g1 / r4) / r5;
        der5[4] = (dx * dy4 * g5 + (6.0 * dx3 * dy2 - 4.0 * dx * dy4) * g4 / r +
                   (3.0 * dx5 + 12.0 * dx * dy4 - 30.0 * dx3 * dy2) * g3 / r2 +
                   (72.0 * dx3 * dy2 - 9.0 * dx5 - 24.0 * dx * dy4) * g2 / r3 +
                   (9.0 * dx5 - 72.0 * dx3 * dy2 + 24.0 * dx * dy4) * g1 / r4) / r5;
        der5[5] = (dy5 * g5 + 10.0 * dx2 * dy3 * g4 / r +
                   (15.0 * dx4 * dy - 30.0 * dx2 * dy3) * g3 / r2 +
                   (60.0 * dx2 * dy3 - 45.0 * dx4 * dy) * g2 / r3 +
                   (45.0 * dx4 * dy - 60.0 * dx2 * dy3) * g1 / r4) / r5;
    }
}

/* ================================================================
 * modbhgreen
 * ================================================================ */
void FNAME(modbhgreen)(
    const double *beta_p, const double *zx, const double *zy,
    const fint *ifpot_p, double *pot,
    const fint *ifgrad_p, double *grad,
    const fint *ifhess_p, double *hess)
{
    double der3[4], der4[5], der5[6];
    fint ifder3 = 0, ifder4 = 0, ifder5 = 0;

    FNAME(modbhgreen_all)(beta_p, zx, zy, ifpot_p, pot, ifgrad_p, grad,
                          ifhess_p, hess, &ifder3, der3, &ifder4, der4,
                          &ifder5, der5);
}

/* ================================================================
 * modbhgreend1
 * ================================================================ */
void FNAME(modbhgreend1)(
    const double *beta_p, const double *zx, const double *zy,
    const fint *ifpot_p, double *pot,
    const fint *ifgrad_p, double *grad,
    const fint *ifhess_p, double *hess,
    const double *dir1)
{
    fint ifpot = *ifpot_p, ifgrad = *ifgrad_p, ifhess = *ifhess_p;

    double potloc, gradloc[2], hessloc[3], der3[4], der4[5], der5[6];
    fint ifpotloc = 0, ifgradloc = 0, ifhessloc = 0;
    fint ifder3 = 0, ifder4 = 0, ifder5 = 0;

    if (ifpot == 1) ifgradloc = 1;
    if (ifgrad == 1) ifhessloc = 1;
    if (ifhess == 1) ifder3 = 1;

    FNAME(modbhgreen_all)(beta_p, zx, zy, &ifpotloc, &potloc,
                          &ifgradloc, gradloc, &ifhessloc, hessloc,
                          &ifder3, der3, &ifder4, der4, &ifder5, der5);

    if (ifpot == 1) {
        *pot = -dir1[0] * gradloc[0] - dir1[1] * gradloc[1];
    }

    if (ifgrad == 1) {
        grad[0] = -dir1[0] * hessloc[0] - dir1[1] * hessloc[1];
        grad[1] = -dir1[0] * hessloc[1] - dir1[1] * hessloc[2];
    }

    if (ifhess == 1) {
        hess[0] = -dir1[0] * der3[0] - dir1[1] * der3[1];
        hess[1] = -dir1[0] * der3[1] - dir1[1] * der3[2];
        hess[2] = -dir1[0] * der3[2] - dir1[1] * der3[3];
    }
}

/* ================================================================
 * modbhgreend2
 * ================================================================ */
void FNAME(modbhgreend2)(
    const double *beta_p, const double *zx, const double *zy,
    const fint *ifpot_p, double *pot,
    const fint *ifgrad_p, double *grad,
    const fint *ifhess_p, double *hess,
    const double *dir1, const double *dir2)
{
    fint ifpot = *ifpot_p, ifgrad = *ifgrad_p, ifhess = *ifhess_p;

    /* convert two directional derivatives to xx, xy, yy */
    double quadvec[3];
    quadvec[0] = dir1[0] * dir2[0];
    quadvec[1] = dir1[0] * dir2[1] + dir1[1] * dir2[0];
    quadvec[2] = dir1[1] * dir2[1];

    double potloc, gradloc[2], hessloc[3], der3[4], der4[5], der5[6];
    fint ifpot1 = 0, ifgrad1 = 0, ifhess1 = 1;
    fint ifder3 = 1, ifder4 = 1, ifder5 = 0;

    FNAME(modbhgreen_all)(beta_p, zx, zy, &ifpot1, &potloc,
                          &ifgrad1, gradloc, &ifhess1, hessloc,
                          &ifder3, der3, &ifder4, der4, &ifder5, der5);

    if (ifpot == 1) {
        *pot = quadvec[0] * hessloc[0] + quadvec[1] * hessloc[1]
             + quadvec[2] * hessloc[2];
    }

    if (ifgrad == 1) {
        grad[0] = quadvec[0] * der3[0] + quadvec[1] * der3[1]
                + quadvec[2] * der3[2];
        grad[1] = quadvec[0] * der3[1] + quadvec[1] * der3[2]
                + quadvec[2] * der3[3];
    }

    if (ifhess == 1) {
        hess[0] = quadvec[0] * der4[0] + quadvec[1] * der4[1]
                + quadvec[2] * der4[2];
        hess[1] = quadvec[0] * der4[1] + quadvec[1] * der4[2]
                + quadvec[2] * der4[3];
        hess[2] = quadvec[0] * der4[2] + quadvec[1] * der4[3]
                + quadvec[2] * der4[4];
    }
}

/* ================================================================
 * difflogbk
 * ================================================================ */
void FNAME(difflogbk)(
    const double *x_p, const double *beta_p,
    const fint *if0_p, double *g0_p,
    const fint *if1_p, double *g1_p,
    const fint *if2_p, double *g2_p,
    const fint *if3_p, double *g3_p)
{
    double x = *x_p, beta = *beta_p;
    (void)if0_p; (void)if1_p; (void)if2_p; (void)if3_p;

    static const double gamma_em = 0.5772156649015328606;
    static const int pmax = 12;

    double y = x * beta;
    double yh = 0.5 * y;
    double yh2 = yh * yh;
    double yh3 = yh2 * yh;
    double lnyh = log(yh);

    double dfac2 = 1;
    double i0 = 1;
    double i1 = 1;
    double i2 = 1 / 2.0;
    double i3 = 1 / 6.0;
    double psi0 = -gamma_em;
    double psi1 = psi0 + 1 / 2.0;
    double psi2 = psi1 + 1 / 4.0;
    double psi3 = psi2 + 1 / 6.0;
    double k0ps = i0 * psi0;
    double k1ps = i1 * psi1;
    double k2ps = i2 * psi2;
    double k3ps = i3 * psi3;

    double yh2p = yh2;
    for (int p = 1; p <= pmax; p++) {
        double dtemp = yh2p / dfac2;
        psi0 = psi0 + 1.0 / p;
        i0 = i0 + dtemp;
        k0ps = k0ps + dtemp * psi0;

        dtemp = dtemp / (p + 1);
        psi1 = psi0 + 0.5 / (p + 1);
        i1 = i1 + dtemp;
        k1ps = k1ps + dtemp * psi1;

        dtemp = dtemp / (p + 2);
        psi2 = psi1 + 0.5 / (p + 2);
        i2 = i2 + dtemp;
        k2ps = k2ps + dtemp * psi2;

        dtemp = dtemp / (p + 3);
        psi3 = psi2 + 0.5 / (p + 3);
        i3 = i3 + dtemp;
        k3ps = k3ps + dtemp * psi3;

        yh2p = yh2p * yh2;
        dfac2 = dfac2 * (p + 1) * (p + 1);
    }

    i1 = yh * i1;
    i2 = yh2 * i2;
    i3 = yh3 * i3;

    double k0 = k0ps - lnyh * i0;
    double k1 = -yh * k1ps + lnyh * i1 + 1 / y;
    double k2 = yh2 * k2ps - lnyh * i2 + 0.5 / yh2 - 0.5;
    double k3 = -yh3 * k3ps + lnyh * i3 + 1 / yh3 - 0.5 / yh + 0.25 * yh;
    (void)k0; (void)k1; (void)k2; (void)k3;

    double betah = beta * 0.5;
    double lnbetah = log(betah);
    (void)lnbetah;

    /* g0 = k0+log(x) */
    *g0_p = k0ps - lnyh * (i0 - 1) - lnbetah;

    /* g1 = -beta*k1+1/x */
    *g1_p = beta * (yh * k1ps - lnyh * i1);

    /* g2 = beta^2/2*(k0+k2)-1/x^2 */
    *g2_p = beta * beta * 0.5 * (k0ps + yh2 * k2ps - lnyh * (i0 + i2) - 0.5);

    /* g3 = -beta^3/4*(3*k1+k3)+2/x^3 */
    *g3_p = beta * beta * beta * (0.75 * yh * k1ps + 0.25 * yh3 * k3ps
                                  - lnyh * (0.75 * i1 + 0.25 * i3)
                                  - 0.25 / yh - 0.25 * 0.25 * yh);
}

/* ================================================================
 * diffslogbk
 * ================================================================ */
void FNAME(diffslogbk)(
    const double *x_p, const double *beta_p,
    const double *rscale_p, double *diffs,
    const fint *n_p)
{
    double x = *x_p, beta = *beta_p, rscale = *rscale_p;
    fint n = *n_p;

    static const double gamma_em = 0.5772156649015328606;
    static const int pmax = 12;

    double y = x * beta;
    double yh = 0.5 * y;
    double lnyh = log(yh);

    if (yh <= 1.0) {
        double yh2 = yh * yh;

        /* psi(1..pmax+n+1) */
        int psilen = pmax + n + 2;
        double *psi = (double *)malloc(psilen * sizeof(double));
        psi[0] = -gamma_em;
        for (int j = 1; j < psilen; j++) {
            psi[j] = psi[j - 1] + 1.0 / j;
        }

        double *kps = (double *)calloc(n + 1, sizeof(double));
        double *ips = (double *)calloc(n + 1, sizeof(double));
        for (int j = 0; j <= n; j++) {
            diffs[j] = 0.0;
        }

        ips[0] = -1.0;

        double yh2p = 1.0;
        for (int p = 0; p <= pmax; p++) {
            double dtemp = yh2p;
            double psitemp1 = psi[p];
            for (int j = 0; j <= n; j++) {
                double psitemp2 = psi[p + j];
                kps[j] = kps[j] + (psitemp1 + psitemp2) * dtemp;
                ips[j] = ips[j] + dtemp;
                dtemp = dtemp / (j + p + 1.0);
            }
            yh2p = yh2p * yh2 / ((p + 1.0) * (p + 1.0));
        }

        double yhp = yh;
        int isign = -1;
        double rscalep = rscale;
        for (int j = 1; j <= n; j++) {
            diffs[j] = isign * (0.5 * yhp * kps[j] - lnyh * ips[j] * yhp) * rscalep;
            yhp = yh * yhp;
            rscalep = rscalep * rscale;
            isign = -isign;
        }

        diffs[0] = 0.5 * kps[0] - lnyh * ips[0] - log(0.5 * beta);

        double *fac = (double *)malloc((n + 1) * sizeof(double));
        fac[0] = 1.0;
        if (n >= 1) fac[1] = 1.0;
        for (int j = 2; j <= n; j++) {
            fac[j] = j * fac[j - 1];
        }

        double rscale2 = rscale * rscale;
        for (int j = 1; j <= n; j++) {
            double dtemp = -0.5 * pow(yh / rscale, 2 - j) * rscale2;
            for (int p = 1; p <= j - 1; p++) {
                diffs[j] = diffs[j] + fac[j - p - 1] * dtemp / fac[p];
                dtemp = -dtemp * yh2;
            }
        }

        free(psi);
        free(kps);
        free(ips);
        free(fac);

    } else {
        /* yh > 1 branch: use h2dall */
        fint nt = n + 5;
        fcomplex *hvec = (fcomplex *)malloc((nt + 1) * sizeof(fcomplex));
        fcomplex *hder = (fcomplex *)malloc((nt + 1) * sizeof(fcomplex));
        double *kvec = (double *)malloc((nt + 1) * sizeof(double));

        fcomplex eye = I;
        fcomplex z = eye * y;
        fint ifder = 0;
        h2dall_(&nt, &z, &rscale, hvec, &ifder, hder);

        double pih = 2.0 * atan(1.0);

        for (int j = 0; j <= n + 1; j += 4) {
            kvec[j]     = -cimag(hvec[j])     * pih;
            kvec[j + 1] = -creal(hvec[j + 1]) * pih;
            kvec[j + 2] =  cimag(hvec[j + 2]) * pih;
            kvec[j + 3] =  creal(hvec[j + 3]) * pih;
        }

        diffs[0] = kvec[0] + log(x);

        double dfac2 = 1.0 / yh;
        double dtemp = 1.0 / y;

        double rscalep = rscale;
        for (int j = 1; j <= n; j++) {
            diffs[j] = (kvec[j] - dtemp) * rscalep;
            dtemp = dtemp * dfac2 * j;
            rscalep = rscalep * rscale;
        }

        free(hvec);
        free(hder);
        free(kvec);
    }
}

/* ================================================================
 * diffslogbk_fast
 * ================================================================ */
void FNAME(diffslogbk_fast)(
    const double *x_p, const double *beta_p,
    const double *rscale_p, double *diffs,
    const fint *ifders_p, double *ders, double *kvec,
    const fint *n_p)
{
    double x = *x_p, beta = *beta_p, rscale = *rscale_p;
    fint n = *n_p, ifders = *ifders_p;

    static const int pmax = 12;

    static const double coeffs0[13] = {
        -1.1544313298030657,
         0.84556867019693427,
         1.8455686701969343,
         2.5122353368636010,
         3.0122353368636010,
         3.4122353368636009,
         3.7455686701969344,
         4.0312829559112204,
         4.2812829559112204,
         4.5035051781334428,
         4.7035051781334429,
         4.8853233599516246,
         5.0519900266182916
    };

    static const double coeffs1[13] = {
        -0.15443132980306573,
         1.3455686701969343,
         2.1789020035302675,
         2.7622353368636010,
         3.2122353368636007,
         3.5789020035302679,
         3.8884258130540772,
         4.1562829559112204,
         4.3923940670223320,
         4.6035051781334424,
         4.7944142690425338,
         4.9686566932849576,
         5.1289131035413682
    };

    double y = x * beta;
    double yh = 0.5 * y;
    double yh2 = yh * yh;
    double lnyh = log(yh);
    double betah_loc = beta / 2.0;

    fcomplex eye = I;
    fcomplex z = eye * y;
    fint ifder = 0;

    if (yh > 1.0 || n > 1 || ifders == 1) {
        fint nt = n + 4;
        fcomplex *hvec = (fcomplex *)malloc((nt + 1) * sizeof(fcomplex));
        fcomplex *hder_arr = (fcomplex *)malloc((nt + 1) * sizeof(fcomplex));

        h2dall_(&nt, &z, &rscale, hvec, &ifder, hder_arr);

        double pih = 2.0 * atan(1.0);
        for (int j = 0; j <= n + 1; j += 4) {
            kvec[j]     = -cimag(hvec[j])     * pih;
            kvec[j + 1] = -creal(hvec[j + 1]) * pih;
            kvec[j + 2] =  cimag(hvec[j + 2]) * pih;
            kvec[j + 3] =  creal(hvec[j + 3]) * pih;
        }

        free(hvec);
        free(hder_arr);
    }

    if (yh <= 1.0) {
        double kps0 = 0.0;
        double ips0 = -1.0;
        double kps1 = 0.0;
        double ips1 = 0.0;

        double yh2p = 1.0;
        for (int p = 0; p <= pmax; p++) {
            kps0 = kps0 + coeffs0[p] * yh2p;
            ips0 = ips0 + yh2p;
            yh2p = yh2p / (p + 1.0);
            kps1 = kps1 + coeffs1[p] * yh2p;
            ips1 = ips1 + yh2p;
            yh2p = yh2p * yh2 / (p + 1.0);
        }

        diffs[0] = 0.5 * kps0 - lnyh * ips0 - log(betah_loc);
        diffs[1] = (lnyh * ips1 * yh - 0.5 * yh * kps1) * rscale;

    } else {
        diffs[0] = kvec[0] + log(x);

        double dtemp = 1.0 / y;
        diffs[1] = kvec[1] - dtemp * rscale;
    }

    double rscale2 = rscale * rscale;

    double dtemp = diffs[1];
    for (int j = 2; j <= n; j++) {
        if (fabs(dtemp) < 1.0e-200) dtemp = 0.0;
        dtemp = 2.0 * (j - 1.0) * dtemp * (rscale / y)
              + kvec[j - 2] * rscale2;
        diffs[j] = dtemp;
    }

    if (ifders == 1) {
        double diffsnp1 = 2.0 * n * diffs[n] * (rscale / y) + kvec[n - 1] * rscale2;
        double betah = 0.5 * beta;

        ders[0] = -beta * diffs[1] / rscale;

        for (int j = 1; j <= n - 1; j++) {
            ders[j] = -betah * (diffs[j + 1] / rscale + kvec[j - 1] * rscale);
        }

        ders[n] = -betah * (diffsnp1 / rscale + kvec[n - 1] * rscale);
    }
}

/* ================================================================
 * diffszkik
 * ================================================================ */
void FNAME(diffszkik)(
    const double *x_p, const double *beta_p,
    const double *rscale_p, double *diffs,
    const fint *n_p)
{
    double x = *x_p, beta = *beta_p, rscale = *rscale_p;
    fint n = *n_p;

    static const int pmax = 12;

    double y = x * beta;
    double yh = 0.5 * y;

    if (yh <= 1.0) {
        double yh2 = yh * yh;

        for (int j = 0; j <= n; j++) {
            diffs[j] = 0.0;
        }

        double yh2p = yh2;
        for (int p = 1; p <= pmax; p++) {
            double dtemp = yh2p;
            for (int j = 0; j <= n; j++) {
                diffs[j] = diffs[j] + dtemp;
                dtemp = dtemp / (j + p + 1.0);
            }
            yh2p = yh2p * yh2 / ((p + 1.0) * (p + 1.0));
        }

        double yh2p2 = 1.0;
        for (int j = 0; j <= n; j++) {
            diffs[j] = yh2p2 * diffs[j];
            yh2p2 = yh2p2 * yh / rscale;
        }

    } else {
        fint lwfjs = n + 5 + 4 * n + 100;

        fcomplex *fjs = (fcomplex *)malloc((lwfjs + 1) * sizeof(fcomplex));
        double *ivals = (double *)malloc((n + 6) * sizeof(double));

        fcomplex eye = I;
        fcomplex z = eye * y;
        fint ifder = 0;
        double fjder[10];

        jbessel2d_(&n, &z, &rscale, fjs, &ifder, fjder);

        for (int j = 0; j <= n + 1; j += 4) {
            ivals[j]     =  creal(fjs[j]);
            ivals[j + 1] =  cimag(fjs[j + 1]);
            ivals[j + 2] = -creal(fjs[j + 2]);
            ivals[j + 3] = -cimag(fjs[j + 3]);
        }

        diffs[0] = ivals[0] - 1.0;

        double dtemp = yh / rscale;
        for (int j = 1; j <= n; j++) {
            diffs[j] = ivals[j] - dtemp;
            dtemp = dtemp * yh / (rscale * (j + 1));
        }

        free(fjs);
        free(ivals);
    }
}

/* ================================================================
 * diffszkik_fast
 * ================================================================ */
void FNAME(diffszkik_fast)(
    const double *x_p, const double *beta_p,
    const double *rscale_p, double *diffs,
    const fint *ifders_p, double *ders, double *ivec,
    const fint *n_p)
{
    double x = *x_p, beta = *beta_p, rscale = *rscale_p;
    fint n = *n_p, ifders = *ifders_p;

    static const int pmax = 12;

    double y = x * beta;
    double yh = 0.5 * y;

    fint lwfjs = n + 5 + 4 * n + 500;

    fcomplex *fjs = (fcomplex *)malloc((lwfjs + 1) * sizeof(fcomplex));
    double *ivals = (double *)malloc((n + 6) * sizeof(double));

    fcomplex eye = I;
    fcomplex z = eye * y;
    fint ifder = 0;
    double fjder[10];

    if (cabs(z) > 300) {
        for (int j = 0; j <= n + 5; j++) {
            fjs[j] = 0.0;
        }
    } else {
        fint nt = n + 5;
        jbessel2d_(&nt, &z, &rscale, fjs, &ifder, fjder);
    }

    for (int j = 0; j <= n + 1; j += 4) {
        ivals[j]     =  creal(fjs[j]);
        ivals[j + 1] =  cimag(fjs[j + 1]);
        ivals[j + 2] = -creal(fjs[j + 2]);
        ivals[j + 3] = -cimag(fjs[j + 3]);
    }

    for (int j = 0; j <= n; j++) {
        ivec[j] = ivals[j];
    }

    if (yh <= 1.0) {
        double yh2 = yh * yh;

        for (int j = 0; j <= n; j++) {
            diffs[j] = 0.0;
        }

        double yh2p = yh2;
        for (int p = 1; p <= pmax; p++) {
            double dtemp = yh2p;
            for (int j = 0; j <= n; j++) {
                diffs[j] = diffs[j] + dtemp;
                if (dtemp < 1.0e-200) dtemp = 0.0;
                dtemp = dtemp / (j + p + 1.0);
            }
            if (yh2p < 1.0e-200) yh2p = 0.0;
            yh2p = yh2p * yh2 / ((p + 1.0) * (p + 1.0));
        }

        double yh2p2 = 1.0;
        for (int j = 0; j <= n; j++) {
            diffs[j] = yh2p2 * diffs[j];
            if (yh2p2 < 1.0e-200) yh2p2 = 0.0;
            yh2p2 = yh2p2 * yh / rscale;
        }

        if (ifders == 1) {
            double betah = 0.5 * beta;
            ders[0] = beta * ivals[1] * rscale;

            for (int j = 1; j <= n; j++) {
                ders[j] = betah * (diffs[j - 1] / rscale + rscale * ivals[j + 1]);
            }
        }

    } else {
        diffs[0] = ivals[0] - 1.0;

        double dtemp = yh / rscale;
        for (int j = 1; j <= n; j++) {
            diffs[j] = ivals[j] - dtemp;
            dtemp = dtemp * yh / (rscale * (j + 1));
            if (dtemp > 1.0e250) dtemp = 0.0;
        }

        if (ifders == 1) {
            double betah = 0.5 * beta;
            ders[0] = beta * ivals[1] * rscale;

            for (int j = 1; j <= n; j++) {
                ders[j] = betah * (diffs[j - 1] / rscale + rscale * ivals[j + 1]);
            }
        }
    }

    free(fjs);
    free(ivals);
}
