/*
 * wideband2d.c - C translation of src/helmholtz/wideband2d.f
 *
 * High-frequency acceleration routines for the Helmholtz FMM.
 * Uses FFT-based diagonal translations for efficiency at large
 * box-to-wavelength ratios.
 *
 * FFT calls (zfftf_, zfftb_) are cross-file dependencies resolved
 * against the Fortran library.
 */

#include "wideband2d.h"

/* Cross-file FFT dependencies */
extern void zfftf_(const fint *n, fcomplex *c, fcomplex *wsave);
extern void zfftb_(const fint *n, fcomplex *c, fcomplex *wsave);

/* Cross-file dependencies */
extern void h2cart2polar_(const double *zat, double *r, double *theta);
extern void jbessel2d_(const fint *nterms, const fcomplex *z,
                       const double *rscale, fcomplex *fjs,
                       const fint *ifder, fcomplex *fjder);
extern void h2dall_(const fint *nterms, const fcomplex *z,
                    const double *rscale, fcomplex *hvec,
                    const fint *ifder, fcomplex *hder);

/* Index into mpole(nd, -nterms:nterms): ii 1-based, n in [-nt,nt] */
#define MIDX(ii, n, nd, nt) ((ii) - 1 + ((n) + (nt)) * (nd))

/* ================================================================
 * h2d_diagtrans: diagonal translation (sig2 incremented)
 * ================================================================ */
void FNAME(h2d_diagtrans)(const fint *nd_p, const fint *nsig_p,
    const fcomplex *sig, const fcomplex *transvec, fcomplex *sig2)
{
    fint nd_v = *nd_p, nsig_v = *nsig_p;
    fint i, ii;

    for (i = 1; i <= nsig_v; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            sig2[FA2(ii,i,nd_v)] += transvec[i-1] * sig[FA2(ii,i,nd_v)];
        }
    }
}

/* ================================================================
 * h2d_mptosig: convert expansion to signature via FFT
 * ================================================================ */
void FNAME(h2d_mptosig)(const fint *nd_p, const fint *nterms1_p,
    const fint *nsig_p, const fcomplex *hexp, fcomplex *sig,
    fcomplex *wsave)
{
    fint nd_v = *nd_p, nt1 = *nterms1_p, nsig_v = *nsig_p;
    fint j, ii;
    fcomplex *sigtmp = (fcomplex *)malloc(nsig_v * sizeof(fcomplex));

    for (ii = 1; ii <= nd_v; ii++) {
        for (j = 0; j < nsig_v; j++) {
            sigtmp[j] = 0;
        }
        for (j = 0; j <= nt1; j++) {
            sigtmp[j] = hexp[MIDX(ii,j,nd_v,nt1)];
        }
        for (j = 1; j <= nt1; j++) {
            sigtmp[nsig_v - j] = hexp[MIDX(ii,-j,nd_v,nt1)];
        }

        zfftf_(&nsig_v, sigtmp, wsave);
        for (j = 1; j <= nsig_v; j++) {
            sig[FA2(ii,j,nd_v)] = sigtmp[j-1];
        }
    }

    free(sigtmp);
}

/* ================================================================
 * h2d_sig2exp: convert signature to expansion via inverse FFT
 * ================================================================ */
void FNAME(h2d_sig2exp)(const fint *nd_p, const fint *nsig_p,
    fcomplex *sig, fcomplex *wsave, const fint *nterms_p,
    fcomplex *expans)
{
    fint nd_v = *nd_p, nsig_v = *nsig_p, nt_v = *nterms_p;
    fint j, ii;
    fcomplex *sig2 = (fcomplex *)malloc(nsig_v * sizeof(fcomplex));

    for (ii = 1; ii <= nd_v; ii++) {
        for (j = 0; j < nsig_v; j++) {
            sig2[j] = sig[FA2(ii,j+1,nd_v)];
        }
        zfftb_(&nsig_v, sig2, wsave);

        for (j = 0; j <= nt_v; j++) {
            expans[MIDX(ii,j,nd_v,nt_v)] += sig2[j];
        }
        for (j = 1; j <= nt_v; j++) {
            expans[MIDX(ii,-j,nd_v,nt_v)] += sig2[nsig_v - j];
        }
    }

    free(sig2);
}

/* ================================================================
 * h2d_mkmpshift: create shift translation operator via FFT
 * ================================================================ */
void FNAME(h2d_mkmpshift)(const fcomplex *zk_p, const double *center1,
    const fint *nterms1_p, const double *center2, const fint *nterms2_p,
    const fint *nsig_p, fcomplex *wsave, fcomplex *transvec)
{
    fint nt1 = *nterms1_p, nt2 = *nterms2_p, nsig_v = *nsig_p;
    fcomplex zk_v = *zk_p;
    fcomplex ima = I;
    fint nterms = nt1 + nt2;
    double done = 1, pi = 4 * atan(done);
    double zdiff[2], r, theta, rscale1 = 1;
    fcomplex z, zmul, zinv, ztemp1, ztemp2;
    fint j, ifder;

    fcomplex *jval = (fcomplex *)malloc((nterms + 4) * sizeof(fcomplex));
    fcomplex *jder = (fcomplex *)malloc((nterms + 4) * sizeof(fcomplex));
    fcomplex *jtemp = (fcomplex *)malloc((2 * nterms + 7) * sizeof(fcomplex));
    /* jtemp offset: index n maps to jtemp[n + nterms + 3] */
    #define JTL(n) ((n) + nterms + 3)

    zdiff[0] = center2[0] - center1[0];
    zdiff[1] = center2[1] - center1[1];
    h2cart2polar_(zdiff, &r, &theta);
    theta = theta - pi;
    z = zk_v * r;

    ifder = 0;
    jbessel2d_(&nterms, &z, &rscale1, jval, &ifder, jder);

    jtemp[JTL(0)] = jval[0];
    zmul = cexp(-ima * theta);
    zinv = conj(zmul);
    ztemp1 = zmul;
    ztemp2 = -zinv;
    for (j = 1; j <= nterms; j++) {
        jtemp[JTL(j)] = ztemp1 * jval[j];
        jtemp[JTL(-j)] = ztemp2 * jval[j];
        ztemp1 = ztemp1 * zmul;
        ztemp2 = -ztemp2 * zinv;
    }

    /* compute translation operator by FFT */
    for (j = 0; j < nsig_v; j++) {
        transvec[j] = 0;
    }
    for (j = 0; j <= nterms; j++) {
        transvec[j] = jtemp[JTL(j)] / nsig_v;
    }
    for (j = 1; j <= nterms; j++) {
        transvec[nsig_v - j] = jtemp[JTL(-j)] / nsig_v;
    }
    zfftf_(&nsig_v, transvec, wsave);

    free(jval); free(jder); free(jtemp);
    #undef JTL
}

/* ================================================================
 * h2d_mkm2ltrans: create M2L translation operator via FFT
 * ================================================================ */
void FNAME(h2d_mkm2ltrans)(const fcomplex *zk_p, const double *center1,
    const fint *nterms1_p, const double *center2, const fint *nterms2_p,
    const fint *nsig_p, fcomplex *wsave, fcomplex *transvec)
{
    fint nt1 = *nterms1_p, nt2 = *nterms2_p, nsig_v = *nsig_p;
    fcomplex zk_v = *zk_p;
    fcomplex ima = I;
    fint nterms = nt1 + nt2;
    double done = 1, pi = 4 * atan(done);
    double zdiff[2], r, theta, rscale1 = 1;
    fcomplex z, zmul, zinv, ztemp1, ztemp2;
    fint j, ifder;
    fint ntj = nterms + 1;

    fcomplex *hval = (fcomplex *)malloc((nterms + 6) * sizeof(fcomplex));
    fcomplex *hder = (fcomplex *)malloc((nterms + 6) * sizeof(fcomplex));
    fcomplex *htemp = (fcomplex *)malloc((2 * nterms + 11) * sizeof(fcomplex));
    #define HTL(n) ((n) + nterms + 5)

    zdiff[0] = center2[0] - center1[0];
    zdiff[1] = center2[1] - center1[1];
    h2cart2polar_(zdiff, &r, &theta);
    theta = theta - pi;
    z = zk_v * r;

    ifder = 0;
    h2dall_(&ntj, &z, &rscale1, hval, &ifder, hder);

    htemp[HTL(0)] = hval[0];
    zmul = cexp(-ima * theta);
    zinv = conj(zmul);
    ztemp1 = zmul;
    ztemp2 = -zinv;
    for (j = 1; j <= nterms; j++) {
        htemp[HTL(j)] = ztemp1 * hval[j];
        htemp[HTL(-j)] = ztemp2 * hval[j];
        ztemp1 = ztemp1 * zmul;
        ztemp2 = -ztemp2 * zinv;
    }

    /* compute transfer function by FFT */
    for (j = 0; j < nsig_v; j++) {
        transvec[j] = 0;
    }
    for (j = 0; j <= nterms; j++) {
        transvec[j] = htemp[HTL(j)] / nsig_v;
    }
    for (j = 1; j <= nterms; j++) {
        transvec[nsig_v - j] = htemp[HTL(-j)] / nsig_v;
    }

    zfftf_(&nsig_v, transvec, wsave);

    free(hval); free(hder); free(htemp);
    #undef HTL
}

/* ================================================================
 * h2dmpmphf: multipole to multipole (HF, via FFT)
 * ================================================================ */
void FNAME(h2dmpmphf)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale1, const double *center1, const fcomplex *hexp1,
    const fint *nterms1_p, const double *rscale2, const double *center2,
    fcomplex *sig2, const fint *nterms2_p, const fint *nsig_p,
    fcomplex *wsave, const fcomplex *transvec)
{
    fint nd_v = *nd_p, nt1 = *nterms1_p, nsig_v = *nsig_p;
    fint j, ii;

    fcomplex *sig = (fcomplex *)malloc(nd_v * nsig_v * sizeof(fcomplex));

    FNAME(h2d_mptosig)(&nd_v, &nt1, &nsig_v, hexp1, sig, wsave);
    FNAME(h2d_diagtrans)(&nd_v, &nsig_v, sig, transvec, sig2);

    free(sig);
}

/* ================================================================
 * h2dloclochf: local to local (HF, via FFT)
 * ================================================================ */
void FNAME(h2dloclochf)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale1, const double *center1, const fcomplex *sig,
    const fint *nterms1_p, const fint *nsig_p, const double *rscale2,
    const double *center2, fcomplex *hexp2, const fint *nterms2_p,
    const fcomplex *transvec, fcomplex *wsave)
{
    fint nd_v = *nd_p, nt2 = *nterms2_p, nsig_v = *nsig_p;
    fint j, ii;

    fcomplex *sig2 = (fcomplex *)malloc(nd_v * nsig_v * sizeof(fcomplex));

    for (ii = 1; ii <= nd_v; ii++) {
        for (j = 1; j <= nsig_v; j++) {
            sig2[FA2(ii,j,nd_v)] = 0.0;
        }
    }

    FNAME(h2d_diagtrans)(&nd_v, &nsig_v, sig, transvec, sig2);
    FNAME(h2d_sig2exp)(&nd_v, &nsig_v, sig2, wsave, &nt2, hexp2);

    free(sig2);
}

/* ================================================================
 * h2dmplochf: multipole to local (HF, via FFT)
 * ================================================================ */
void FNAME(h2dmplochf)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale1, const double *center1, const fcomplex *sig,
    const fint *nterms1_p, const double *rscale2, const double *center2,
    fcomplex *sig2, const fint *nterms2_p, const fint *nsig_p,
    fcomplex *wsave, const fcomplex *tvec)
{
    fint nd_v = *nd_p, nsig_v = *nsig_p;

    FNAME(h2d_diagtrans)(&nd_v, &nsig_v, sig, tvec, sig2);
}
