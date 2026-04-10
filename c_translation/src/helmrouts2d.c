/*
 * helmrouts2d.c - C translation of src/helmholtz/helmrouts2d.f
 *
 * Core FMM operations for 2D Helmholtz: multipole/local expansion
 * formation, evaluation, and translation operators.
 *
 * Strict 1:1 translation. All expansions use INCREMENTAL updates.
 *
 * Array indexing conventions:
 *   mpole(nd, -nterms:nterms) -> MIDX(ii, n, nd, nt)
 *   mpolex(nd, -nterms-1:nterms+1) -> MX(ii, n, nd, nt)
 *   mpolexx(nd, -nterms-2:nterms+2) -> MXX(ii, n, nd, nt)
 *   mptemp(-nterms-2:nterms+2) -> MT(n, nt)
 *   jtemp/htemp(-nterms-5:nterms+5) -> JT(n, nt)
 */

#include "helmrouts2d.h"

/* Index into mpole(nd, -nterms:nterms): ii 1-based, n in [-nt,nt] */
#define MIDX(ii, n, nd, nt) ((ii) - 1 + ((n) + (nt)) * (nd))
/* Index into mpolex(nd, -nterms-1:nterms+1) */
#define MX(ii, n, nd, nt)   ((ii) - 1 + ((n) + (nt) + 1) * (nd))
/* Index into mpolexx(nd, -nterms-2:nterms+2) */
#define MXX(ii, n, nd, nt)  ((ii) - 1 + ((n) + (nt) + 2) * (nd))
/* Index into mptemp(-nterms-2:nterms+2), 1D */
#define MT(n, nt)            ((n) + (nt) + 2)
/* Index into jtemp/htemp(-nterms-5:nterms+5), 1D */
#define JT(n, nt)            ((n) + (nt) + 5)

/* Cross-file dependencies */
extern void h2cart2polar_(const double *zat, double *r, double *theta);
extern void h2dall_(const fint *nterms, const fcomplex *z,
                    const double *rscale, fcomplex *hvec,
                    const fint *ifder, fcomplex *hder);
extern void jbessel2d_(const fint *nterms, const fcomplex *z,
                       const double *rscale, fcomplex *fjs,
                       const fint *ifder, fcomplex *fjder);

/* ================================================================
 * Internal helpers (same-file, use FNAME)
 * ================================================================ */

static void FNAME(ctompole)(const fint *nd_p, const fcomplex *zmul_p,
    const fcomplex *zinv_p, fcomplex *mpole, const fcomplex *jval,
    const fcomplex *charge, const fint *nterms_p)
{
    fint nd_v = *nd_p, nt_v = *nterms_p;
    fcomplex zmul_v = *zmul_p, zinv_v = *zinv_p;
    fint n, ii;
    fcomplex ztemp1, ztemp2;

    for (ii = 1; ii <= nd_v; ii++) {
        mpole[MIDX(ii,0,nd_v,nt_v)] += charge[ii-1] * jval[0];
    }
    ztemp1 = zmul_v;
    ztemp2 = -zinv_v;
    for (n = 1; n <= nt_v; n++) {
        for (ii = 1; ii <= nd_v; ii++) {
            mpole[MIDX(ii,n,nd_v,nt_v)] += jval[n] * ztemp1 * charge[ii-1];
            mpole[MIDX(ii,-n,nd_v,nt_v)] += jval[n] * ztemp2 * charge[ii-1];
        }
        ztemp1 = ztemp1 * zmul_v;
        ztemp2 = -ztemp2 * zinv_v;
    }
}

static void FNAME(dtompole)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale_p, const fcomplex *zmul_p, const fcomplex *zinv_p,
    fcomplex *mpole, const fcomplex *jval, const fcomplex *dipstr,
    const double *dipvec, const fint *nterms_p)
{
    fint nd_v = *nd_p, nt_v = *nterms_p;
    fcomplex zk_v = *zk_p, zmul_v = *zmul_p, zinv_v = *zinv_p;
    double rscale_v = *rscale_p;
    fcomplex ima = I;
    fint n, ii;
    fcomplex ztemp1, ztemp2, ztemp3, ztemp4, ztemp5, ztemp6;
    fcomplex ztemp7, ztemp8, ztemp9, ztemp10;

    ztemp3 = zinv_v / rscale_v;
    ztemp4 = zmul_v * rscale_v;
    ztemp5 = zmul_v / rscale_v;
    ztemp6 = zinv_v * rscale_v;
    for (ii = 1; ii <= nd_v; ii++) {
        mpole[MIDX(ii,0,nd_v,nt_v)] += -dipstr[ii-1] * jval[1] * zk_v / 2 * rscale_v *
            ((zmul_v + zinv_v) * dipvec[FA3(ii,1,1,nd_v,2) /* but dipvec is nd,2 for this source */]
            + (zmul_v - zinv_v) * ima * dipvec[FA3(ii,2,1,nd_v,2)]);
    }
    ztemp1 = -zmul_v * zk_v / 2;
    ztemp2 = +zinv_v * zk_v / 2;
    for (n = 1; n <= nt_v; n++) {
        for (ii = 1; ii <= nd_v; ii++) {
            ztemp7 = ztemp3 * (-dipvec[FA3(ii,1,1,nd_v,2)] + ima * dipvec[FA3(ii,2,1,nd_v,2)]);
            ztemp8 = ztemp4 * (+dipvec[FA3(ii,1,1,nd_v,2)] + ima * dipvec[FA3(ii,2,1,nd_v,2)]);
            ztemp9 = ztemp5 * (-dipvec[FA3(ii,1,1,nd_v,2)] - ima * dipvec[FA3(ii,2,1,nd_v,2)]);
            ztemp10 = ztemp6 * (+dipvec[FA3(ii,1,1,nd_v,2)] - ima * dipvec[FA3(ii,2,1,nd_v,2)]);
            mpole[MIDX(ii,n,nd_v,nt_v)] +=
                (jval[n-1] * ztemp7 + jval[n+1] * ztemp8) * ztemp1 * dipstr[ii-1];
            mpole[MIDX(ii,-n,nd_v,nt_v)] +=
                (jval[n-1] * ztemp9 + jval[n+1] * ztemp10) * ztemp2 * dipstr[ii-1];
        }
        ztemp1 = ztemp1 * zmul_v;
        ztemp2 = -ztemp2 * zinv_v;
    }
}

static void FNAME(mpole_evalp)(const fint *nd_p, const fcomplex *zmul_p,
    const fcomplex *zinv_p, const fcomplex *mpole, fcomplex *mptemp,
    const fcomplex *hval, const fint *nterms_p, fcomplex *pot1)
{
    fint nd_v = *nd_p, nt_v = *nterms_p;
    fcomplex zmul_v = *zmul_p, zinv_v = *zinv_p;
    fcomplex ima = I, ima4inv = I / 4.0;
    fint j, n, ii;
    fcomplex ztemp1, ztemp2;

    mptemp[MT(0,nt_v)] = hval[0];
    ztemp1 = zmul_v * ima4inv;
    ztemp2 = -zinv_v * ima4inv;
    for (j = 1; j <= nt_v + 2; j++) {
        mptemp[MT(j,nt_v)] = ztemp1 * hval[j];
        mptemp[MT(-j,nt_v)] = ztemp2 * hval[j];
        ztemp1 = ztemp1 * zmul_v;
        ztemp2 = -ztemp2 * zinv_v;
    }
    for (ii = 1; ii <= nd_v; ii++) {
        pot1[ii-1] += hval[0] * mpole[MIDX(ii,0,nd_v,nt_v)] * ima4inv;
    }
    for (n = 1; n <= nt_v; n++) {
        for (ii = 1; ii <= nd_v; ii++) {
            /* Split: pot1 = pot1 + mpole(n)*mptemp(n) + mpole(-n)*mptemp(-n) */
            pot1[ii-1] += mpole[MIDX(ii,n,nd_v,nt_v)] * mptemp[MT(n,nt_v)];
            pot1[ii-1] += mpole[MIDX(ii,-n,nd_v,nt_v)] * mptemp[MT(-n,nt_v)];
        }
    }
}

static void FNAME(mpole_evalg)(const fint *nd_p, const fcomplex *mpolex,
    const fcomplex *mpoley, const fcomplex *mptemp, const fint *nterms_p,
    fcomplex *grad1)
{
    fint nd_v = *nd_p, nt_v = *nterms_p;
    fcomplex ima4inv = I / 4.0;
    fint n, ii;

    for (ii = 1; ii <= nd_v; ii++) {
        grad1[FA2(ii,1,nd_v)] += mptemp[MT(0,nt_v)] * mpolex[MX(ii,0,nd_v,nt_v)] * ima4inv;
        grad1[FA2(ii,2,nd_v)] += mptemp[MT(0,nt_v)] * mpoley[MX(ii,0,nd_v,nt_v)] * ima4inv;
    }
    for (n = 1; n <= nt_v + 1; n++) {
        for (ii = 1; ii <= nd_v; ii++) {
            /* Split multi-term additions */
            grad1[FA2(ii,1,nd_v)] += mpolex[MX(ii,n,nd_v,nt_v)] * mptemp[MT(n,nt_v)];
            grad1[FA2(ii,1,nd_v)] += mpolex[MX(ii,-n,nd_v,nt_v)] * mptemp[MT(-n,nt_v)];
            grad1[FA2(ii,2,nd_v)] += mpoley[MX(ii,n,nd_v,nt_v)] * mptemp[MT(n,nt_v)];
            grad1[FA2(ii,2,nd_v)] += mpoley[MX(ii,-n,nd_v,nt_v)] * mptemp[MT(-n,nt_v)];
        }
    }
}

static void FNAME(mpole_evalh)(const fint *nd_p, const fcomplex *mpolexx,
    const fcomplex *mpolexy, const fcomplex *mpoleyy,
    const fcomplex *mptemp, const fint *nterms_p, fcomplex *hess1)
{
    fint nd_v = *nd_p, nt_v = *nterms_p;
    fcomplex ima4inv = I / 4.0;
    fint n, ii;

    for (ii = 1; ii <= nd_v; ii++) {
        hess1[FA2(ii,1,nd_v)] += mptemp[MT(0,nt_v)] * mpolexx[MXX(ii,0,nd_v,nt_v)] * ima4inv;
        hess1[FA2(ii,2,nd_v)] += mptemp[MT(0,nt_v)] * mpolexy[MXX(ii,0,nd_v,nt_v)] * ima4inv;
        hess1[FA2(ii,3,nd_v)] += mptemp[MT(0,nt_v)] * mpoleyy[MXX(ii,0,nd_v,nt_v)] * ima4inv;
    }
    for (n = 1; n <= nt_v + 2; n++) {
        for (ii = 1; ii <= nd_v; ii++) {
            /* Split multi-term additions */
            hess1[FA2(ii,1,nd_v)] += mpolexx[MXX(ii,n,nd_v,nt_v)] * mptemp[MT(n,nt_v)];
            hess1[FA2(ii,1,nd_v)] += mpolexx[MXX(ii,-n,nd_v,nt_v)] * mptemp[MT(-n,nt_v)];
            hess1[FA2(ii,2,nd_v)] += mpolexy[MXX(ii,n,nd_v,nt_v)] * mptemp[MT(n,nt_v)];
            hess1[FA2(ii,2,nd_v)] += mpolexy[MXX(ii,-n,nd_v,nt_v)] * mptemp[MT(-n,nt_v)];
            hess1[FA2(ii,3,nd_v)] += mpoleyy[MXX(ii,n,nd_v,nt_v)] * mptemp[MT(n,nt_v)];
            hess1[FA2(ii,3,nd_v)] += mpoleyy[MXX(ii,-n,nd_v,nt_v)] * mptemp[MT(-n,nt_v)];
        }
    }
}

static void FNAME(mk_mpoleg)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale_p, const fcomplex *mpole, fcomplex *mpolex,
    fcomplex *mpoley, const fint *nterms_p)
{
    fint nd_v = *nd_p, nt_v = *nterms_p;
    fcomplex zk_v = *zk_p;
    double rscale_v = *rscale_p;
    fcomplex ima = I;
    fint i, ii;
    fcomplex z1scale, z2scale, z3scale, z4scale;

    for (i = -nt_v - 1; i <= nt_v + 1; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            mpolex[MX(ii,i,nd_v,nt_v)] = 0;
            mpoley[MX(ii,i,nd_v,nt_v)] = 0;
        }
    }
    z1scale = zk_v / 2 / rscale_v;
    z2scale = zk_v / 2 * rscale_v;
    z3scale = zk_v / 2 / rscale_v * ima;
    z4scale = zk_v / 2 * rscale_v * ima;

    for (ii = 1; ii <= nd_v; ii++) {
        mpolex[MX(ii,-1,nd_v,nt_v)] += z1scale * mpole[MIDX(ii,0,nd_v,nt_v)];
        mpoley[MX(ii,-1,nd_v,nt_v)] += z3scale * mpole[MIDX(ii,0,nd_v,nt_v)];
        mpolex[MX(ii,1,nd_v,nt_v)] += -z1scale * mpole[MIDX(ii,0,nd_v,nt_v)];
        mpoley[MX(ii,1,nd_v,nt_v)] += z3scale * mpole[MIDX(ii,0,nd_v,nt_v)];
    }

    for (i = 1; i <= nt_v; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            mpolex[MX(ii,i-1,nd_v,nt_v)] += z2scale * mpole[MIDX(ii,i,nd_v,nt_v)];
            mpoley[MX(ii,i-1,nd_v,nt_v)] += z4scale * mpole[MIDX(ii,i,nd_v,nt_v)];
            mpolex[MX(ii,i+1,nd_v,nt_v)] += -z1scale * mpole[MIDX(ii,i,nd_v,nt_v)];
            mpoley[MX(ii,i+1,nd_v,nt_v)] += z3scale * mpole[MIDX(ii,i,nd_v,nt_v)];
        }
    }

    for (i = -nt_v; i <= -1; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            mpolex[MX(ii,i-1,nd_v,nt_v)] += z1scale * mpole[MIDX(ii,i,nd_v,nt_v)];
            mpoley[MX(ii,i-1,nd_v,nt_v)] += z3scale * mpole[MIDX(ii,i,nd_v,nt_v)];
            mpolex[MX(ii,i+1,nd_v,nt_v)] += -z2scale * mpole[MIDX(ii,i,nd_v,nt_v)];
            mpoley[MX(ii,i+1,nd_v,nt_v)] += z4scale * mpole[MIDX(ii,i,nd_v,nt_v)];
        }
    }
}

static void FNAME(mk_mpoleh)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale_p, const fcomplex *mpolex, const fcomplex *mpoley,
    fcomplex *mpolexx, fcomplex *mpolexy, fcomplex *mpoleyy,
    const fint *nterms_p)
{
    fint nd_v = *nd_p, nt_v = *nterms_p;
    fcomplex zk_v = *zk_p;
    double rscale_v = *rscale_p;
    fcomplex ima = I;
    fint i, ii;
    fcomplex z1scale, z2scale, z3scale, z4scale;

    for (i = -nt_v - 2; i <= nt_v + 2; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            mpolexx[MXX(ii,i,nd_v,nt_v)] = 0;
            mpolexy[MXX(ii,i,nd_v,nt_v)] = 0;
            mpoleyy[MXX(ii,i,nd_v,nt_v)] = 0;
        }
    }
    z1scale = zk_v / 2 / rscale_v;
    z2scale = zk_v / 2 * rscale_v;
    z3scale = zk_v / 2 / rscale_v * ima;
    z4scale = zk_v / 2 * rscale_v * ima;

    for (ii = 1; ii <= nd_v; ii++) {
        mpolexx[MXX(ii,-1,nd_v,nt_v)] += z1scale * mpolex[MX(ii,0,nd_v,nt_v)];
        mpolexy[MXX(ii,-1,nd_v,nt_v)] += z3scale * mpolex[MX(ii,0,nd_v,nt_v)];
        mpoleyy[MXX(ii,-1,nd_v,nt_v)] += z3scale * mpoley[MX(ii,0,nd_v,nt_v)];
        mpolexx[MXX(ii,1,nd_v,nt_v)] += -z1scale * mpolex[MX(ii,0,nd_v,nt_v)];
        mpolexy[MXX(ii,1,nd_v,nt_v)] += z3scale * mpolex[MX(ii,0,nd_v,nt_v)];
        mpoleyy[MXX(ii,1,nd_v,nt_v)] += z3scale * mpoley[MX(ii,0,nd_v,nt_v)];
    }
    for (i = 1; i <= nt_v + 1; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            mpolexx[MXX(ii,i-1,nd_v,nt_v)] += z2scale * mpolex[MX(ii,i,nd_v,nt_v)];
            mpolexy[MXX(ii,i-1,nd_v,nt_v)] += z4scale * mpolex[MX(ii,i,nd_v,nt_v)];
            mpoleyy[MXX(ii,i-1,nd_v,nt_v)] += z4scale * mpoley[MX(ii,i,nd_v,nt_v)];
            mpolexx[MXX(ii,i+1,nd_v,nt_v)] += -z1scale * mpolex[MX(ii,i,nd_v,nt_v)];
            mpolexy[MXX(ii,i+1,nd_v,nt_v)] += z3scale * mpolex[MX(ii,i,nd_v,nt_v)];
            mpoleyy[MXX(ii,i+1,nd_v,nt_v)] += z3scale * mpoley[MX(ii,i,nd_v,nt_v)];
        }
    }
    for (i = -nt_v - 1; i <= -1; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            mpolexx[MXX(ii,i-1,nd_v,nt_v)] += z1scale * mpolex[MX(ii,i,nd_v,nt_v)];
            mpolexy[MXX(ii,i-1,nd_v,nt_v)] += z3scale * mpolex[MX(ii,i,nd_v,nt_v)];
            mpoleyy[MXX(ii,i-1,nd_v,nt_v)] += z3scale * mpoley[MX(ii,i,nd_v,nt_v)];
            mpolexx[MXX(ii,i+1,nd_v,nt_v)] += -z2scale * mpolex[MX(ii,i,nd_v,nt_v)];
            mpolexy[MXX(ii,i+1,nd_v,nt_v)] += z4scale * mpolex[MX(ii,i,nd_v,nt_v)];
            mpoleyy[MXX(ii,i+1,nd_v,nt_v)] += z4scale * mpoley[MX(ii,i,nd_v,nt_v)];
        }
    }
}

/* ================================================================
 * Main routines
 * ================================================================ */

/* h2dformmpc: multipole from charges */
void FNAME(h2dformmpc)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale_p, const double *source, const fint *ns_p,
    const fcomplex *charge, const double *center, const fint *nterms_p,
    fcomplex *mpole)
{
    fint nd_v = *nd_p, ns_v = *ns_p, nt_v = *nterms_p;
    fcomplex zk_v = *zk_p;
    double rscale_v = *rscale_p;
    fcomplex ima = I;
    fint j, ifder;
    double zdiff[2], r, theta;
    fcomplex z, zmul, zinv;
    fcomplex *jval = (fcomplex *)malloc((nt_v + 6) * sizeof(fcomplex));
    fcomplex *jder = (fcomplex *)malloc((nt_v + 6) * sizeof(fcomplex));
    fint nt1 = nt_v + 1;

    for (j = 1; j <= ns_v; j++) {
        zdiff[0] = source[FA2(1,j,2)] - center[0];
        zdiff[1] = source[FA2(2,j,2)] - center[1];
        h2cart2polar_(zdiff, &r, &theta);
        z = zk_v * r;
        ifder = 0;
        jbessel2d_(&nt1, &z, &rscale_v, jval, &ifder, jder);
        zmul = cexp(-ima * theta);
        zinv = conj(zmul);
        FNAME(ctompole)(&nd_v, &zmul, &zinv, mpole, jval,
                        &charge[FA2(1,j,nd_v)], &nt_v);
    }
    free(jval); free(jder);
}

/* h2dformmpd: multipole from dipoles */
void FNAME(h2dformmpd)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale_p, const double *source, const fint *ns_p,
    const fcomplex *dipstr, const double *dipvec,
    const double *center, const fint *nterms_p, fcomplex *mpole)
{
    fint nd_v = *nd_p, ns_v = *ns_p, nt_v = *nterms_p;
    fcomplex zk_v = *zk_p;
    double rscale_v = *rscale_p;
    fcomplex ima = I;
    fint j, ifder;
    double zdiff[2], r, theta;
    fcomplex z, zmul, zinv;
    fcomplex *jval = (fcomplex *)malloc((nt_v + 6) * sizeof(fcomplex));
    fcomplex *jder = (fcomplex *)malloc((nt_v + 6) * sizeof(fcomplex));
    fint nt1 = nt_v + 1;

    for (j = 1; j <= ns_v; j++) {
        zdiff[0] = source[FA2(1,j,2)] - center[0];
        zdiff[1] = source[FA2(2,j,2)] - center[1];
        h2cart2polar_(zdiff, &r, &theta);
        z = zk_v * r;
        ifder = 0;
        jbessel2d_(&nt1, &z, &rscale_v, jval, &ifder, jder);
        zmul = cexp(-ima * theta);
        zinv = conj(zmul);
        FNAME(dtompole)(&nd_v, &zk_v, &rscale_v, &zmul, &zinv,
                        mpole, jval, &dipstr[FA2(1,j,nd_v)],
                        &dipvec[FA3(1,1,j,nd_v,2)], &nt_v);
    }
    free(jval); free(jder);
}

/* h2dformmpcd: multipole from charges + dipoles */
void FNAME(h2dformmpcd)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale_p, const double *source, const fint *ns_p,
    const fcomplex *charge, const fcomplex *dipstr, const double *dipvec,
    const double *center, const fint *nterms_p, fcomplex *mpole)
{
    fint nd_v = *nd_p, ns_v = *ns_p, nt_v = *nterms_p;
    fcomplex zk_v = *zk_p;
    double rscale_v = *rscale_p;
    fcomplex ima = I;
    fint j, ifder;
    double zdiff[2], r, theta;
    fcomplex z, zmul, zinv;
    fcomplex *jval = (fcomplex *)malloc((nt_v + 6) * sizeof(fcomplex));
    fcomplex *jder = (fcomplex *)malloc((nt_v + 6) * sizeof(fcomplex));
    fint nt1 = nt_v + 1;

    for (j = 1; j <= ns_v; j++) {
        zdiff[0] = source[FA2(1,j,2)] - center[0];
        zdiff[1] = source[FA2(2,j,2)] - center[1];
        h2cart2polar_(zdiff, &r, &theta);
        z = zk_v * r;
        ifder = 0;
        jbessel2d_(&nt1, &z, &rscale_v, jval, &ifder, jder);
        zmul = cexp(-ima * theta);
        zinv = conj(zmul);
        FNAME(ctompole)(&nd_v, &zmul, &zinv, mpole, jval,
                        &charge[FA2(1,j,nd_v)], &nt_v);
        FNAME(dtompole)(&nd_v, &zk_v, &rscale_v, &zmul, &zinv,
                        mpole, jval, &dipstr[FA2(1,j,nd_v)],
                        &dipvec[FA3(1,1,j,nd_v,2)], &nt_v);
    }
    free(jval); free(jder);
}

/* h2dmpevalp: multipole -> potential */
void FNAME(h2dmpevalp)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale_p, const double *center, const fcomplex *mpole,
    const fint *nterms_p, const double *ztarg, const fint *ntarg_p,
    fcomplex *pot1)
{
    fint nd_v = *nd_p, nt_v = *nterms_p, ntarg_v = *ntarg_p;
    fcomplex zk_v = *zk_p;
    double rscale_v = *rscale_p;
    fcomplex ima = I;
    fint k, ifder;
    double zdiff[2], r, theta;
    fcomplex z, zmul, zinv;
    fint nt3 = nt_v + 3;
    fcomplex *hval = (fcomplex *)malloc((nt_v + 6) * sizeof(fcomplex));
    fcomplex *hder = (fcomplex *)malloc((nt_v + 6) * sizeof(fcomplex));
    fcomplex *mptemp = (fcomplex *)malloc((2 * nt_v + 5) * sizeof(fcomplex));

    for (k = 1; k <= ntarg_v; k++) {
        zdiff[0] = ztarg[FA2(1,k,2)] - center[0];
        zdiff[1] = ztarg[FA2(2,k,2)] - center[1];
        h2cart2polar_(zdiff, &r, &theta);
        z = zk_v * r;
        ifder = 0;
        h2dall_(&nt3, &z, &rscale_v, hval, &ifder, hder);
        zmul = cexp(ima * theta);
        zinv = conj(zmul);
        FNAME(mpole_evalp)(&nd_v, &zmul, &zinv, mpole, mptemp, hval,
                           &nt_v, &pot1[FA2(1,k,nd_v)]);
    }
    free(hval); free(hder); free(mptemp);
}

/* h2dmpevalg: multipole -> potential + gradient */
void FNAME(h2dmpevalg)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale_p, const double *center, const fcomplex *mpole,
    const fint *nterms_p, const double *ztarg, const fint *ntarg_p,
    fcomplex *pot1, fcomplex *grad1)
{
    fint nd_v = *nd_p, nt_v = *nterms_p, ntarg_v = *ntarg_p;
    fcomplex zk_v = *zk_p;
    double rscale_v = *rscale_p;
    fcomplex ima = I;
    fint k, ifder;
    double zdiff[2], r, theta;
    fcomplex z, zmul, zinv;
    fint nt3 = nt_v + 3;
    fcomplex *hval = (fcomplex *)malloc((nt_v + 6) * sizeof(fcomplex));
    fcomplex *hder = (fcomplex *)malloc((nt_v + 6) * sizeof(fcomplex));
    fcomplex *mpolex = (fcomplex *)malloc(nd_v * (2 * nt_v + 3) * sizeof(fcomplex));
    fcomplex *mpoley = (fcomplex *)malloc(nd_v * (2 * nt_v + 3) * sizeof(fcomplex));
    fcomplex *mptemp = (fcomplex *)malloc((2 * nt_v + 5) * sizeof(fcomplex));

    FNAME(mk_mpoleg)(&nd_v, &zk_v, &rscale_v, mpole, mpolex, mpoley, &nt_v);

    for (k = 1; k <= ntarg_v; k++) {
        zdiff[0] = ztarg[FA2(1,k,2)] - center[0];
        zdiff[1] = ztarg[FA2(2,k,2)] - center[1];
        h2cart2polar_(zdiff, &r, &theta);
        z = zk_v * r;
        ifder = 0;
        h2dall_(&nt3, &z, &rscale_v, hval, &ifder, hder);
        zmul = cexp(ima * theta);
        zinv = conj(zmul);
        FNAME(mpole_evalp)(&nd_v, &zmul, &zinv, mpole, mptemp, hval,
                           &nt_v, &pot1[FA2(1,k,nd_v)]);
        FNAME(mpole_evalg)(&nd_v, mpolex, mpoley, mptemp, &nt_v,
                           &grad1[FA3(1,1,k,nd_v,2)]);
    }
    free(hval); free(hder); free(mpolex); free(mpoley); free(mptemp);
}

/* h2dmpevalh: multipole -> potential + gradient + hessian */
void FNAME(h2dmpevalh)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale_p, const double *center, const fcomplex *mpole,
    const fint *nterms_p, const double *ztarg, const fint *ntarg_p,
    fcomplex *pot1, fcomplex *grad1, fcomplex *hess1)
{
    fint nd_v = *nd_p, nt_v = *nterms_p, ntarg_v = *ntarg_p;
    fcomplex zk_v = *zk_p;
    double rscale_v = *rscale_p;
    fcomplex ima = I;
    fint k, ifder;
    double zdiff[2], r, theta;
    fcomplex z, zmul, zinv;
    fint nt3 = nt_v + 3;
    fcomplex *hval = (fcomplex *)malloc((nt_v + 6) * sizeof(fcomplex));
    fcomplex *hder = (fcomplex *)malloc((nt_v + 6) * sizeof(fcomplex));
    fcomplex *mpolex = (fcomplex *)malloc(nd_v * (2 * nt_v + 3) * sizeof(fcomplex));
    fcomplex *mpoley = (fcomplex *)malloc(nd_v * (2 * nt_v + 3) * sizeof(fcomplex));
    fcomplex *mpolexx = (fcomplex *)malloc(nd_v * (2 * nt_v + 5) * sizeof(fcomplex));
    fcomplex *mpolexy = (fcomplex *)malloc(nd_v * (2 * nt_v + 5) * sizeof(fcomplex));
    fcomplex *mpoleyy = (fcomplex *)malloc(nd_v * (2 * nt_v + 5) * sizeof(fcomplex));
    fcomplex *mptemp = (fcomplex *)malloc((2 * nt_v + 5) * sizeof(fcomplex));

    FNAME(mk_mpoleg)(&nd_v, &zk_v, &rscale_v, mpole, mpolex, mpoley, &nt_v);
    FNAME(mk_mpoleh)(&nd_v, &zk_v, &rscale_v, mpolex, mpoley,
                     mpolexx, mpolexy, mpoleyy, &nt_v);

    for (k = 1; k <= ntarg_v; k++) {
        zdiff[0] = ztarg[FA2(1,k,2)] - center[0];
        zdiff[1] = ztarg[FA2(2,k,2)] - center[1];
        h2cart2polar_(zdiff, &r, &theta);
        z = zk_v * r;
        ifder = 0;
        h2dall_(&nt3, &z, &rscale_v, hval, &ifder, hder);
        zmul = cexp(ima * theta);
        zinv = conj(zmul);
        FNAME(mpole_evalp)(&nd_v, &zmul, &zinv, mpole, mptemp, hval,
                           &nt_v, &pot1[FA2(1,k,nd_v)]);
        FNAME(mpole_evalg)(&nd_v, mpolex, mpoley, mptemp, &nt_v,
                           &grad1[FA3(1,1,k,nd_v,2)]);
        FNAME(mpole_evalh)(&nd_v, mpolexx, mpolexy, mpoleyy, mptemp, &nt_v,
                           &hess1[FA3(1,1,k,nd_v,3)]);
    }
    free(hval); free(hder); free(mpolex); free(mpoley);
    free(mpolexx); free(mpolexy); free(mpoleyy); free(mptemp);
}

/* h2dformtac: local from charges */
void FNAME(h2dformtac)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale_p, const double *source, const fint *ns_p,
    const fcomplex *charge, const double *center, const fint *nterms_p,
    fcomplex *local)
{
    fint nd_v = *nd_p, ns_v = *ns_p, nt_v = *nterms_p;
    fcomplex zk_v = *zk_p;
    double rscale_v = *rscale_p;
    fcomplex ima = I;
    fint j, ifder;
    double zdiff[2], r, theta;
    fcomplex z, zmul, zinv;
    fint nt2 = nt_v + 2;
    fcomplex *hval = (fcomplex *)malloc((nt_v + 6) * sizeof(fcomplex));
    fcomplex *hder = (fcomplex *)malloc((nt_v + 6) * sizeof(fcomplex));

    for (j = 1; j <= ns_v; j++) {
        zdiff[0] = source[FA2(1,j,2)] - center[0];
        zdiff[1] = source[FA2(2,j,2)] - center[1];
        h2cart2polar_(zdiff, &r, &theta);
        z = zk_v * r;
        ifder = 0;
        h2dall_(&nt2, &z, &rscale_v, hval, &ifder, hder);
        zmul = cexp(-ima * theta);
        zinv = conj(zmul);
        FNAME(ctompole)(&nd_v, &zmul, &zinv, local, hval,
                        &charge[FA2(1,j,nd_v)], &nt_v);
    }
    free(hval); free(hder);
}

/* h2dformtad: local from dipoles */
void FNAME(h2dformtad)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale_p, const double *source, const fint *ns_p,
    const fcomplex *dipstr, const double *dipvec,
    const double *center, const fint *nterms_p, fcomplex *local)
{
    fint nd_v = *nd_p, ns_v = *ns_p, nt_v = *nterms_p;
    fcomplex zk_v = *zk_p;
    double rscale_v = *rscale_p;
    fcomplex ima = I;
    fint j, ifder;
    double zdiff[2], r, theta;
    fcomplex z, zmul, zinv;
    double rsinv = 1.0 / rscale_v;
    fint nt2 = nt_v + 2;
    fcomplex *hval = (fcomplex *)malloc((nt_v + 6) * sizeof(fcomplex));
    fcomplex *hder = (fcomplex *)malloc((nt_v + 6) * sizeof(fcomplex));

    for (j = 1; j <= ns_v; j++) {
        zdiff[0] = source[FA2(1,j,2)] - center[0];
        zdiff[1] = source[FA2(2,j,2)] - center[1];
        h2cart2polar_(zdiff, &r, &theta);
        z = zk_v * r;
        ifder = 0;
        h2dall_(&nt2, &z, &rscale_v, hval, &ifder, hder);
        zmul = cexp(-ima * theta);
        zinv = conj(zmul);
        FNAME(dtompole)(&nd_v, &zk_v, &rsinv, &zmul, &zinv,
                        local, hval, &dipstr[FA2(1,j,nd_v)],
                        &dipvec[FA3(1,1,j,nd_v,2)], &nt_v);
    }
    free(hval); free(hder);
}

/* h2dformtacd: local from charges + dipoles */
void FNAME(h2dformtacd)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale_p, const double *source, const fint *ns_p,
    const fcomplex *charge, const fcomplex *dipstr, const double *dipvec,
    const double *center, const fint *nterms_p, fcomplex *local)
{
    fint nd_v = *nd_p, ns_v = *ns_p, nt_v = *nterms_p;
    fcomplex zk_v = *zk_p;
    double rscale_v = *rscale_p;
    fcomplex ima = I;
    fint j, ifder;
    double zdiff[2], r, theta;
    fcomplex z, zmul, zinv;
    double rsinv = 1.0 / rscale_v;
    fint nt2 = nt_v + 2;
    fcomplex *hval = (fcomplex *)malloc((nt_v + 6) * sizeof(fcomplex));
    fcomplex *hder = (fcomplex *)malloc((nt_v + 6) * sizeof(fcomplex));

    for (j = 1; j <= ns_v; j++) {
        zdiff[0] = source[FA2(1,j,2)] - center[0];
        zdiff[1] = source[FA2(2,j,2)] - center[1];
        h2cart2polar_(zdiff, &r, &theta);
        z = zk_v * r;
        ifder = 0;
        h2dall_(&nt2, &z, &rscale_v, hval, &ifder, hder);
        zmul = cexp(-ima * theta);
        zinv = conj(zmul);
        FNAME(ctompole)(&nd_v, &zmul, &zinv, local, hval,
                        &charge[FA2(1,j,nd_v)], &nt_v);
        FNAME(dtompole)(&nd_v, &zk_v, &rsinv, &zmul, &zinv,
                        local, hval, &dipstr[FA2(1,j,nd_v)],
                        &dipvec[FA3(1,1,j,nd_v,2)], &nt_v);
    }
    free(hval); free(hder);
}

/* h2dtaevalp: local -> potential */
void FNAME(h2dtaevalp)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale_p, const double *center, const fcomplex *local,
    const fint *nterms_p, const double *ztarg, const fint *ntarg_p,
    fcomplex *pot1)
{
    fint nd_v = *nd_p, nt_v = *nterms_p, ntarg_v = *ntarg_p;
    fcomplex zk_v = *zk_p;
    double rscale_v = *rscale_p;
    fcomplex ima = I;
    fint k, ifder;
    double zdiff[2], r, theta;
    fcomplex z, zmul, zinv;
    fint nt3 = nt_v + 3;
    fcomplex *jval = (fcomplex *)malloc((nt_v + 11) * sizeof(fcomplex));
    fcomplex *jder = (fcomplex *)malloc((nt_v + 11) * sizeof(fcomplex));
    fcomplex *mptemp = (fcomplex *)malloc((2 * nt_v + 5) * sizeof(fcomplex));

    for (k = 1; k <= ntarg_v; k++) {
        zdiff[0] = ztarg[FA2(1,k,2)] - center[0];
        zdiff[1] = ztarg[FA2(2,k,2)] - center[1];
        h2cart2polar_(zdiff, &r, &theta);
        z = zk_v * r;
        ifder = 0;
        jbessel2d_(&nt3, &z, &rscale_v, jval, &ifder, jder);
        zmul = cexp(ima * theta);
        zinv = conj(zmul);
        FNAME(mpole_evalp)(&nd_v, &zmul, &zinv, local, mptemp, jval,
                           &nt_v, &pot1[FA2(1,k,nd_v)]);
    }
    free(jval); free(jder); free(mptemp);
}

/* h2dtaevalg: local -> potential + gradient */
void FNAME(h2dtaevalg)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale_p, const double *center, const fcomplex *local,
    const fint *nterms_p, const double *ztarg, const fint *ntarg_p,
    fcomplex *pot1, fcomplex *grad1)
{
    fint nd_v = *nd_p, nt_v = *nterms_p, ntarg_v = *ntarg_p;
    fcomplex zk_v = *zk_p;
    double rscale_v = *rscale_p;
    fcomplex ima = I;
    fint k, ifder;
    double zdiff[2], r, theta;
    fcomplex z, zmul, zinv;
    double rsinv = 1.0 / rscale_v;
    fint nt3 = nt_v + 3;
    fcomplex *jval = (fcomplex *)malloc((nt_v + 11) * sizeof(fcomplex));
    fcomplex *jder = (fcomplex *)malloc((nt_v + 11) * sizeof(fcomplex));
    fcomplex *localx = (fcomplex *)malloc(nd_v * (2 * nt_v + 3) * sizeof(fcomplex));
    fcomplex *localy = (fcomplex *)malloc(nd_v * (2 * nt_v + 3) * sizeof(fcomplex));
    fcomplex *mptemp = (fcomplex *)malloc((2 * nt_v + 5) * sizeof(fcomplex));

    FNAME(mk_mpoleg)(&nd_v, &zk_v, &rsinv, local, localx, localy, &nt_v);

    for (k = 1; k <= ntarg_v; k++) {
        zdiff[0] = ztarg[FA2(1,k,2)] - center[0];
        zdiff[1] = ztarg[FA2(2,k,2)] - center[1];
        h2cart2polar_(zdiff, &r, &theta);
        z = zk_v * r;
        ifder = 0;
        jbessel2d_(&nt3, &z, &rscale_v, jval, &ifder, jder);
        zmul = cexp(ima * theta);
        zinv = conj(zmul);
        FNAME(mpole_evalp)(&nd_v, &zmul, &zinv, local, mptemp, jval,
                           &nt_v, &pot1[FA2(1,k,nd_v)]);
        FNAME(mpole_evalg)(&nd_v, localx, localy, mptemp, &nt_v,
                           &grad1[FA3(1,1,k,nd_v,2)]);
    }
    free(jval); free(jder); free(localx); free(localy); free(mptemp);
}

/* h2dtaevalh: local -> potential + gradient + hessian */
void FNAME(h2dtaevalh)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale_p, const double *center, const fcomplex *local,
    const fint *nterms_p, const double *ztarg, const fint *ntarg_p,
    fcomplex *pot1, fcomplex *grad1, fcomplex *hess1)
{
    fint nd_v = *nd_p, nt_v = *nterms_p, ntarg_v = *ntarg_p;
    fcomplex zk_v = *zk_p;
    double rscale_v = *rscale_p;
    fcomplex ima = I;
    fint k, ifder;
    double zdiff[2], r, theta;
    fcomplex z, zmul, zinv;
    double rsinv = 1.0 / rscale_v;
    fint nt3 = nt_v + 3;
    fcomplex *jval = (fcomplex *)malloc((nt_v + 11) * sizeof(fcomplex));
    fcomplex *jder = (fcomplex *)malloc((nt_v + 11) * sizeof(fcomplex));
    fcomplex *localx = (fcomplex *)malloc(nd_v * (2 * nt_v + 3) * sizeof(fcomplex));
    fcomplex *localy = (fcomplex *)malloc(nd_v * (2 * nt_v + 3) * sizeof(fcomplex));
    fcomplex *localxx = (fcomplex *)malloc(nd_v * (2 * nt_v + 5) * sizeof(fcomplex));
    fcomplex *localxy = (fcomplex *)malloc(nd_v * (2 * nt_v + 5) * sizeof(fcomplex));
    fcomplex *localyy = (fcomplex *)malloc(nd_v * (2 * nt_v + 5) * sizeof(fcomplex));
    fcomplex *mptemp = (fcomplex *)malloc((2 * nt_v + 5) * sizeof(fcomplex));

    FNAME(mk_mpoleg)(&nd_v, &zk_v, &rsinv, local, localx, localy, &nt_v);
    FNAME(mk_mpoleh)(&nd_v, &zk_v, &rsinv, localx, localy,
                     localxx, localxy, localyy, &nt_v);

    for (k = 1; k <= ntarg_v; k++) {
        zdiff[0] = ztarg[FA2(1,k,2)] - center[0];
        zdiff[1] = ztarg[FA2(2,k,2)] - center[1];
        h2cart2polar_(zdiff, &r, &theta);
        z = zk_v * r;
        ifder = 0;
        jbessel2d_(&nt3, &z, &rscale_v, jval, &ifder, jder);
        zmul = cexp(ima * theta);
        zinv = conj(zmul);
        FNAME(mpole_evalp)(&nd_v, &zmul, &zinv, local, mptemp, jval,
                           &nt_v, &pot1[FA2(1,k,nd_v)]);
        FNAME(mpole_evalg)(&nd_v, localx, localy, mptemp, &nt_v,
                           &grad1[FA3(1,1,k,nd_v,2)]);
        FNAME(mpole_evalh)(&nd_v, localxx, localxy, localyy, mptemp, &nt_v,
                           &hess1[FA3(1,1,k,nd_v,3)]);
    }
    free(jval); free(jder); free(localx); free(localy);
    free(localxx); free(localxy); free(localyy); free(mptemp);
}

/* h2dmpmp: multipole to multipole translation */
void FNAME(h2dmpmp)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale1_p, const double *center1, const fcomplex *hexp1,
    const fint *nterms1_p, const double *rscale2_p, const double *center2,
    fcomplex *hexp2, const fint *nterms2_p)
{
    fint nd_v = *nd_p, nt1 = *nterms1_p, nt2 = *nterms2_p;
    fcomplex zk_v = *zk_p;
    double rscale1 = *rscale1_p, rscale2 = *rscale2_p;
    fint nterms = nt1 + nt2;
    fcomplex ima = I;
    fint i, j, ii, ifder;
    double zdiff[2], r, theta, pi, done;
    double rsj, rsi5, fs2, rsj2, fs3;
    fcomplex z, zmul, zinv, ztemp1, ztemp2;
    fint ntj = nterms + 3;

    fcomplex *jval = (fcomplex *)malloc((nterms + 11) * sizeof(fcomplex));
    fcomplex *jder = (fcomplex *)malloc((nterms + 11) * sizeof(fcomplex));
    fcomplex *jtemp = (fcomplex *)malloc((2 * nterms + 11) * sizeof(fcomplex));

    done = 1;
    pi = 4 * atan(done);

    zdiff[0] = center2[0] - center1[0];
    zdiff[1] = center2[1] - center1[1];
    h2cart2polar_(zdiff, &r, &theta);
    theta = theta - pi;
    z = zk_v * r;
    ifder = 0;
    jbessel2d_(&ntj, &z, &rscale1, jval, &ifder, jder);

    jtemp[JT(0,nterms)] = jval[0];
    zmul = cexp(-ima * theta);
    zinv = conj(zmul);
    ztemp1 = zmul;
    ztemp2 = -zinv;
    for (j = 1; j <= nterms; j++) {
        jtemp[JT(j,nterms)] = ztemp1 * jval[j];
        jtemp[JT(-j,nterms)] = ztemp2 * jval[j];
        ztemp1 = ztemp1 * zmul;
        ztemp2 = -ztemp2 * zinv;
    }

    for (ii = 1; ii <= nd_v; ii++) {
        hexp2[MIDX(ii,0,nd_v,nt2)] += hexp1[MIDX(ii,0,nd_v,nt1)] * jtemp[JT(0,nterms)];
    }
    rsj = rscale1 * rscale1;
    rsj2 = rsj;
    for (j = 1; j <= nt1; j++) {
        for (ii = 1; ii <= nd_v; ii++) {
            hexp2[MIDX(ii,0,nd_v,nt2)] += (hexp1[MIDX(ii,j,nd_v,nt1)] * jtemp[JT(-j,nterms)]) * rsj;
            hexp2[MIDX(ii,0,nd_v,nt2)] += (hexp1[MIDX(ii,-j,nd_v,nt1)] * jtemp[JT(j,nterms)]) * rsj;
        }
        rsj = rsj * rsj2;
    }

    rsi5 = rscale1 / rscale2;
    for (i = 1; i <= nt2; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            hexp2[MIDX(ii,i,nd_v,nt2)] += hexp1[MIDX(ii,0,nd_v,nt1)] * jtemp[JT(i,nterms)] * rsi5;
            hexp2[MIDX(ii,-i,nd_v,nt2)] += hexp1[MIDX(ii,0,nd_v,nt1)] * jtemp[JT(-i,nterms)] * rsi5;
        }
        rsi5 = rsi5 * rscale1 / rscale2;
    }
    rsi5 = rscale1 / rscale2;
    for (i = 1; i <= nt2; i++) {
        rsj = rsj2;
        for (j = 1; j <= (nt1 < i ? nt1 : i); j++) {
            fs3 = rsi5 * rsj;
            for (ii = 1; ii <= nd_v; ii++) {
                hexp2[MIDX(ii,i,nd_v,nt2)] += (hexp1[MIDX(ii,j,nd_v,nt1)] * jtemp[JT(i-j,nterms)]) * rsi5;
                hexp2[MIDX(ii,i,nd_v,nt2)] += (hexp1[MIDX(ii,-j,nd_v,nt1)] * jtemp[JT(i+j,nterms)]) * fs3;
                hexp2[MIDX(ii,-i,nd_v,nt2)] += (hexp1[MIDX(ii,j,nd_v,nt1)] * jtemp[JT(-i-j,nterms)]) * fs3;
                hexp2[MIDX(ii,-i,nd_v,nt2)] += (hexp1[MIDX(ii,-j,nd_v,nt1)] * jtemp[JT(-i+j,nterms)]) * rsi5;
            }
            rsj = rsj * rsj2;
        }
        {
            double rsj_pow = 1.0;
            for (int p = 0; p < i + 1; p++) rsj_pow *= rsj2;
            fs2 = rsi5 * rsj2;
            fs3 = rsi5 * rsj_pow;
        }
        for (j = i + 1; j <= nt1; j++) {
            for (ii = 1; ii <= nd_v; ii++) {
                hexp2[MIDX(ii,i,nd_v,nt2)] += (hexp1[MIDX(ii,j,nd_v,nt1)] * jtemp[JT(i-j,nterms)]) * fs2;
                hexp2[MIDX(ii,i,nd_v,nt2)] += (hexp1[MIDX(ii,-j,nd_v,nt1)] * jtemp[JT(i+j,nterms)]) * fs3;
                hexp2[MIDX(ii,-i,nd_v,nt2)] += (hexp1[MIDX(ii,j,nd_v,nt1)] * jtemp[JT(-i-j,nterms)]) * fs3;
                hexp2[MIDX(ii,-i,nd_v,nt2)] += (hexp1[MIDX(ii,-j,nd_v,nt1)] * jtemp[JT(-i+j,nterms)]) * fs2;
            }
            fs3 = fs3 * rsj2;
            fs2 = fs2 * rsj2;
        }
        rsi5 = rsi5 * rscale1 / rscale2;
    }
    free(jval); free(jder); free(jtemp);
}

/* h2dlocloc: local to local translation */
void FNAME(h2dlocloc)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale1_p, const double *center1, const fcomplex *jexp1,
    const fint *nterms1_p, const double *rscale2_p, const double *center2,
    fcomplex *jexp2, const fint *nterms2_p)
{
    fint nd_v = *nd_p, nt1 = *nterms1_p, nt2 = *nterms2_p;
    fcomplex zk_v = *zk_p;
    double rscale1 = *rscale1_p, rscale2 = *rscale2_p;
    fint nterms = nt1 + nt2;
    fcomplex ima = I;
    fint i, j, ii, ifder;
    double zdiff[2], r, theta, pi, done;
    double rsi, rsj, rsi7, rsi5, fs2;
    fcomplex z, zmul, zinv, ztemp1, ztemp2;
    fint ntj = nterms + 3;

    fcomplex *jval = (fcomplex *)malloc((nterms + 11) * sizeof(fcomplex));
    fcomplex *jder = (fcomplex *)malloc((nterms + 11) * sizeof(fcomplex));
    fcomplex *jtemp = (fcomplex *)malloc((2 * nterms + 11) * sizeof(fcomplex));

    done = 1;
    pi = 4 * atan(done);

    zdiff[0] = center2[0] - center1[0];
    zdiff[1] = center2[1] - center1[1];
    h2cart2polar_(zdiff, &r, &theta);
    theta = theta - pi;
    z = zk_v * r;
    ifder = 0;
    jbessel2d_(&ntj, &z, &rscale1, jval, &ifder, jder);

    jtemp[JT(0,nterms)] = jval[0];
    zmul = cexp(-ima * theta);
    zinv = conj(zmul);
    ztemp1 = zmul;
    ztemp2 = -zinv;
    for (j = 1; j <= nterms; j++) {
        jtemp[JT(j,nterms)] = ztemp1 * jval[j];
        jtemp[JT(-j,nterms)] = ztemp2 * jval[j];
        ztemp1 = ztemp1 * zmul;
        ztemp2 = -ztemp2 * zinv;
    }

    for (ii = 1; ii <= nd_v; ii++) {
        jexp2[MIDX(ii,0,nd_v,nt2)] += jexp1[MIDX(ii,0,nd_v,nt1)] * jtemp[JT(0,nterms)];
    }
    for (j = 1; j <= nt1; j++) {
        for (ii = 1; ii <= nd_v; ii++) {
            jexp2[MIDX(ii,0,nd_v,nt2)] += (jexp1[MIDX(ii,j,nd_v,nt1)] * jtemp[JT(-j,nterms)]);
            jexp2[MIDX(ii,0,nd_v,nt2)] += (jexp1[MIDX(ii,-j,nd_v,nt1)] * jtemp[JT(j,nterms)]);
        }
    }

    rsi = rscale1;
    rsi7 = rscale2;
    rsi5 = rscale2 / rscale1;
    for (i = 1; i <= nt2; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            jexp2[MIDX(ii,i,nd_v,nt2)] += jexp1[MIDX(ii,0,nd_v,nt1)] * jtemp[JT(i,nterms)] * rsi7 * rsi;
            jexp2[MIDX(ii,-i,nd_v,nt2)] += jexp1[MIDX(ii,0,nd_v,nt1)] * jtemp[JT(-i,nterms)] * rsi7 * rsi;
        }
        rsi = rsi * rscale1;
        rsi7 = rsi7 * rscale2;
    }
    rsi = rscale1;
    rsi7 = rscale2;
    rsi5 = rscale2 / rscale1;
    for (i = 1; i <= nt2; i++) {
        fs2 = rsi5;
        if (nt1 <= i) {
            double pw = 1.0;
            for (int p = 0; p < 2 * (i - nt1); p++) pw *= rscale1;
            fs2 = fs2 * pw;
        }
        for (j = (nt1 < i ? nt1 : i); j >= 1; j--) {
            for (ii = 1; ii <= nd_v; ii++) {
                jexp2[MIDX(ii,i,nd_v,nt2)] += (jexp1[MIDX(ii,j,nd_v,nt1)] * jtemp[JT(i-j,nterms)]) * fs2;
                jexp2[MIDX(ii,i,nd_v,nt2)] += (jexp1[MIDX(ii,-j,nd_v,nt1)] * jtemp[JT(i+j,nterms)]) * rsi7 * rsi;
                jexp2[MIDX(ii,-i,nd_v,nt2)] += (jexp1[MIDX(ii,j,nd_v,nt1)] * jtemp[JT(-i-j,nterms)]) * rsi7 * rsi;
                jexp2[MIDX(ii,-i,nd_v,nt2)] += (jexp1[MIDX(ii,-j,nd_v,nt1)] * jtemp[JT(-i+j,nterms)]) * fs2;
            }
            fs2 = fs2 * rscale1 * rscale1;
        }
        rsj = 1.0;
        for (int p = 0; p < i + 1; p++) rsj *= rscale1;
        for (j = i + 1; j <= nt1; j++) {
            for (ii = 1; ii <= nd_v; ii++) {
                jexp2[MIDX(ii,i,nd_v,nt2)] += (jexp1[MIDX(ii,j,nd_v,nt1)] * jtemp[JT(i-j,nterms)]) * rsi5;
                jexp2[MIDX(ii,i,nd_v,nt2)] += (jexp1[MIDX(ii,-j,nd_v,nt1)] * jtemp[JT(i+j,nterms)]) * rsi7 * rsi;
                jexp2[MIDX(ii,-i,nd_v,nt2)] += (jexp1[MIDX(ii,j,nd_v,nt1)] * jtemp[JT(-i-j,nterms)]) * rsi7 * rsi;
                jexp2[MIDX(ii,-i,nd_v,nt2)] += (jexp1[MIDX(ii,-j,nd_v,nt1)] * jtemp[JT(-i+j,nterms)]) * rsi5;
            }
            rsj = rsj * rscale1;
        }
        rsi = rsi * rscale1;
        rsi7 = rsi7 * rscale2;
        rsi5 = rsi5 * rscale2 / rscale1;
    }
    free(jval); free(jder); free(jtemp);
}

/* h2dmploc: multipole to local translation */
void FNAME(h2dmploc)(const fint *nd_p, const fcomplex *zk_p,
    const double *rscale1_p, const double *center1, const fcomplex *hexp,
    const fint *nterms1_p, const double *rscale2_p, const double *center2,
    fcomplex *jexp, const fint *nterms2_p)
{
    fint nd_v = *nd_p, nt1 = *nterms1_p, nt2 = *nterms2_p;
    fcomplex zk_v = *zk_p;
    double rscale1 = *rscale1_p, rscale2 = *rscale2_p;
    fint nterms = nt1 + nt2;
    fcomplex ima = I;
    fint i, j, ii, ifder;
    double zdiff[2], r, theta, pi, done, rfac, rs12;
    fcomplex z, zmul, zinv, ztemp1, ztemp2;
    fint ntj = nterms + 1;

    fcomplex *hval = (fcomplex *)malloc((nterms + 6) * sizeof(fcomplex));
    fcomplex *hder = (fcomplex *)malloc((nterms + 6) * sizeof(fcomplex));
    fcomplex *htemp = (fcomplex *)malloc((2 * nterms + 11) * sizeof(fcomplex));
    fcomplex *jexptemp = (fcomplex *)malloc(nd_v * (2 * nt2 + 1) * sizeof(fcomplex));
    double *rs2pow = (double *)malloc((nt2 + 1) * sizeof(double));

    done = 1;
    pi = 4 * atan(done);

    zdiff[0] = center2[0] - center1[0];
    zdiff[1] = center2[1] - center1[1];
    h2cart2polar_(zdiff, &r, &theta);
    theta = theta - pi;
    z = zk_v * r;
    ifder = 0;
    h2dall_(&ntj, &z, &rscale1, hval, &ifder, hder);

    htemp[JT(0,nterms)] = hval[0];
    zmul = cexp(-ima * theta);
    zinv = conj(zmul);
    ztemp1 = zmul;
    ztemp2 = -zinv;
    for (j = 1; j <= nterms; j++) {
        htemp[JT(j,nterms)] = ztemp1 * hval[j];
        htemp[JT(-j,nterms)] = ztemp2 * hval[j];
        ztemp1 = ztemp1 * zmul;
        ztemp2 = -ztemp2 * zinv;
    }

    for (i = -nt2; i <= nt2; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            jexptemp[MIDX(ii,i,nd_v,nt2)] = 0;
        }
    }

    rfac = rscale1 * rscale1;
    rs2pow[0] = rfac; /* rs2pow(1) in Fortran = rs2pow[0] in C */
    for (i = 1; i < nt2; i++) {
        rs2pow[i] = rs2pow[i - 1] * rfac;
    }

    for (ii = 1; ii <= nd_v; ii++) {
        jexptemp[MIDX(ii,0,nd_v,nt2)] += hexp[MIDX(ii,0,nd_v,nt1)] * htemp[JT(0,nterms)];
    }
    for (j = 1; j <= nt1; j++) {
        for (ii = 1; ii <= nd_v; ii++) {
            jexptemp[MIDX(ii,0,nd_v,nt2)] += (hexp[MIDX(ii,j,nd_v,nt1)] * htemp[JT(-j,nterms)]);
            jexptemp[MIDX(ii,0,nd_v,nt2)] += (hexp[MIDX(ii,-j,nd_v,nt1)] * htemp[JT(j,nterms)]);
        }
    }

    for (i = 1; i <= nt2; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            jexptemp[MIDX(ii,i,nd_v,nt2)] += hexp[MIDX(ii,0,nd_v,nt1)] * htemp[JT(i,nterms)];
            jexptemp[MIDX(ii,-i,nd_v,nt2)] += hexp[MIDX(ii,0,nd_v,nt1)] * htemp[JT(-i,nterms)];
        }
        for (j = 1; j <= (nt1 < i ? nt1 : i); j++) {
            for (ii = 1; ii <= nd_v; ii++) {
                jexptemp[MIDX(ii,i,nd_v,nt2)] += (hexp[MIDX(ii,j,nd_v,nt1)] * htemp[JT(i-j,nterms)]) * rs2pow[j-1];
                jexptemp[MIDX(ii,i,nd_v,nt2)] += (hexp[MIDX(ii,-j,nd_v,nt1)] * htemp[JT(i+j,nterms)]);
                jexptemp[MIDX(ii,-i,nd_v,nt2)] += (hexp[MIDX(ii,j,nd_v,nt1)] * htemp[JT(-i-j,nterms)]);
                jexptemp[MIDX(ii,-i,nd_v,nt2)] += (hexp[MIDX(ii,-j,nd_v,nt1)] * htemp[JT(-i+j,nterms)]) * rs2pow[j-1];
            }
        }
        for (j = i + 1; j <= nt1; j++) {
            for (ii = 1; ii <= nd_v; ii++) {
                jexptemp[MIDX(ii,i,nd_v,nt2)] += (hexp[MIDX(ii,j,nd_v,nt1)] * htemp[JT(i-j,nterms)]) * rs2pow[i-1];
                jexptemp[MIDX(ii,i,nd_v,nt2)] += (hexp[MIDX(ii,-j,nd_v,nt1)] * htemp[JT(i+j,nterms)]);
                jexptemp[MIDX(ii,-i,nd_v,nt2)] += (hexp[MIDX(ii,j,nd_v,nt1)] * htemp[JT(-i-j,nterms)]);
                jexptemp[MIDX(ii,-i,nd_v,nt2)] += (hexp[MIDX(ii,-j,nd_v,nt1)] * htemp[JT(-i+j,nterms)]) * rs2pow[i-1];
            }
        }
    }

    rfac = rscale2 / rscale1;
    rs12 = rfac;

    for (ii = 1; ii <= nd_v; ii++) {
        jexp[MIDX(ii,0,nd_v,nt2)] += jexptemp[MIDX(ii,0,nd_v,nt2)];
    }

    for (i = 1; i <= nt2; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            jexp[MIDX(ii,i,nd_v,nt2)] += jexptemp[MIDX(ii,i,nd_v,nt2)] * rs12;
            jexp[MIDX(ii,-i,nd_v,nt2)] += jexptemp[MIDX(ii,-i,nd_v,nt2)] * rs12;
        }
        rs12 = rs12 * rfac;
    }

    free(hval); free(hder); free(htemp); free(jexptemp); free(rs2pow);
}
