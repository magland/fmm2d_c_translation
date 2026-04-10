/*
 * dfft_threadsafe.c - C translation of src/common/dfft_threadsafe.f
 *
 * Complex FFT routines from FFTPACK version 4 (April 1985), by
 * Paul N. Swarztrauber, NCAR. Modified for thread safety by
 * Leslie Greengard (1/28/2023).
 *
 * Only the complex-FFT subset (zffti/zfftf/zfftb) and the butterfly
 * helpers they call (dpassb/f 2/3/4/5 and general) are translated.
 * Real FFT, sine/cosine transforms are not needed by fmm2d.
 *
 * All routines operate on flat double arrays with 1-based Fortran
 * indexing convention. Complex data c(1:2n) stores n complex values
 * as interleaved (re,im) pairs.
 *
 * Strict 1:1 translation — no refactoring.
 */

#include "dfft_threadsafe.h"

/* 3D col-major index: arr(i,j,k) with dimensions (ld1,ld2,*) */
#define IDX3(i, j, k, ld1, ld2) \
    ((((k) - 1) * (ld2) + ((j) - 1)) * (ld1) + ((i) - 1))

/* 2D col-major index: arr(i,j) with dimension (ld1,*) */
#define IDX2(i, j, ld1) (((j) - 1) * (ld1) + ((i) - 1))

/* ================================================================
 * dpassb2 - backward butterfly, radix 2
 * ================================================================ */
void FNAME(dpassb2)(const fint *ido_p, const fint *l1_p,
                    const double *cc, double *ch, const double *wa1)
{
    fint ido_v = *ido_p;
    fint l1_v = *l1_p;
    fint k, i;
    double tr2, ti2;

    if (ido_v > 2) goto L102;
    for (k = 1; k <= l1_v; k++) {
        ch[IDX3(1,k,1, ido_v,l1_v)] = cc[IDX3(1,1,k, ido_v,2)] + cc[IDX3(1,2,k, ido_v,2)];
        ch[IDX3(1,k,2, ido_v,l1_v)] = cc[IDX3(1,1,k, ido_v,2)] - cc[IDX3(1,2,k, ido_v,2)];
        ch[IDX3(2,k,1, ido_v,l1_v)] = cc[IDX3(2,1,k, ido_v,2)] + cc[IDX3(2,2,k, ido_v,2)];
        ch[IDX3(2,k,2, ido_v,l1_v)] = cc[IDX3(2,1,k, ido_v,2)] - cc[IDX3(2,2,k, ido_v,2)];
    }
    return;
L102:
    for (k = 1; k <= l1_v; k++) {
        for (i = 2; i <= ido_v; i += 2) {
            ch[IDX3(i-1,k,1, ido_v,l1_v)] = cc[IDX3(i-1,1,k, ido_v,2)] + cc[IDX3(i-1,2,k, ido_v,2)];
            tr2 = cc[IDX3(i-1,1,k, ido_v,2)] - cc[IDX3(i-1,2,k, ido_v,2)];
            ch[IDX3(i,k,1, ido_v,l1_v)] = cc[IDX3(i,1,k, ido_v,2)] + cc[IDX3(i,2,k, ido_v,2)];
            ti2 = cc[IDX3(i,1,k, ido_v,2)] - cc[IDX3(i,2,k, ido_v,2)];
            ch[IDX3(i,k,2, ido_v,l1_v)] = wa1[i-1-1]*ti2 + wa1[i-1]*tr2;
            ch[IDX3(i-1,k,2, ido_v,l1_v)] = wa1[i-1-1]*tr2 - wa1[i-1]*ti2;
        }
    }
}

/* ================================================================
 * dpassb3 - backward butterfly, radix 3
 * ================================================================ */
void FNAME(dpassb3)(const fint *ido_p, const fint *l1_p,
                    const double *cc, double *ch,
                    const double *wa1, const double *wa2)
{
    fint ido_v = *ido_p;
    fint l1_v = *l1_p;
    fint k, i;
    double tr2, cr2, ti2, ci2, cr3, ci3, dr2, dr3, di2, di3;

    static const double taur = -0.5;
    static const double taui = 0.86602540378443864676372317075293618;

    if (ido_v != 2) goto L102;
    for (k = 1; k <= l1_v; k++) {
        tr2 = cc[IDX3(1,2,k, ido_v,3)] + cc[IDX3(1,3,k, ido_v,3)];
        cr2 = cc[IDX3(1,1,k, ido_v,3)] + taur*tr2;
        ch[IDX3(1,k,1, ido_v,l1_v)] = cc[IDX3(1,1,k, ido_v,3)] + tr2;
        ti2 = cc[IDX3(2,2,k, ido_v,3)] + cc[IDX3(2,3,k, ido_v,3)];
        ci2 = cc[IDX3(2,1,k, ido_v,3)] + taur*ti2;
        ch[IDX3(2,k,1, ido_v,l1_v)] = cc[IDX3(2,1,k, ido_v,3)] + ti2;
        cr3 = taui*(cc[IDX3(1,2,k, ido_v,3)] - cc[IDX3(1,3,k, ido_v,3)]);
        ci3 = taui*(cc[IDX3(2,2,k, ido_v,3)] - cc[IDX3(2,3,k, ido_v,3)]);
        ch[IDX3(1,k,2, ido_v,l1_v)] = cr2 - ci3;
        ch[IDX3(1,k,3, ido_v,l1_v)] = cr2 + ci3;
        ch[IDX3(2,k,2, ido_v,l1_v)] = ci2 + cr3;
        ch[IDX3(2,k,3, ido_v,l1_v)] = ci2 - cr3;
    }
    return;
L102:
    for (k = 1; k <= l1_v; k++) {
        for (i = 2; i <= ido_v; i += 2) {
            tr2 = cc[IDX3(i-1,2,k, ido_v,3)] + cc[IDX3(i-1,3,k, ido_v,3)];
            cr2 = cc[IDX3(i-1,1,k, ido_v,3)] + taur*tr2;
            ch[IDX3(i-1,k,1, ido_v,l1_v)] = cc[IDX3(i-1,1,k, ido_v,3)] + tr2;
            ti2 = cc[IDX3(i,2,k, ido_v,3)] + cc[IDX3(i,3,k, ido_v,3)];
            ci2 = cc[IDX3(i,1,k, ido_v,3)] + taur*ti2;
            ch[IDX3(i,k,1, ido_v,l1_v)] = cc[IDX3(i,1,k, ido_v,3)] + ti2;
            cr3 = taui*(cc[IDX3(i-1,2,k, ido_v,3)] - cc[IDX3(i-1,3,k, ido_v,3)]);
            ci3 = taui*(cc[IDX3(i,2,k, ido_v,3)] - cc[IDX3(i,3,k, ido_v,3)]);
            dr2 = cr2 - ci3;
            dr3 = cr2 + ci3;
            di2 = ci2 + cr3;
            di3 = ci2 - cr3;
            ch[IDX3(i,k,2, ido_v,l1_v)] = wa1[i-1-1]*di2 + wa1[i-1]*dr2;
            ch[IDX3(i-1,k,2, ido_v,l1_v)] = wa1[i-1-1]*dr2 - wa1[i-1]*di2;
            ch[IDX3(i,k,3, ido_v,l1_v)] = wa2[i-1-1]*di3 + wa2[i-1]*dr3;
            ch[IDX3(i-1,k,3, ido_v,l1_v)] = wa2[i-1-1]*dr3 - wa2[i-1]*di3;
        }
    }
}

/* ================================================================
 * dpassb4 - backward butterfly, radix 4
 * ================================================================ */
void FNAME(dpassb4)(const fint *ido_p, const fint *l1_p,
                    const double *cc, double *ch,
                    const double *wa1, const double *wa2,
                    const double *wa3)
{
    fint ido_v = *ido_p;
    fint l1_v = *l1_p;
    fint k, i;
    double ti1, ti2, tr4, ti3, tr1, tr2, ti4, tr3;
    double cr2, cr3, cr4, ci2, ci3, ci4;

    if (ido_v != 2) goto L102;
    for (k = 1; k <= l1_v; k++) {
        ti1 = cc[IDX3(2,1,k, ido_v,4)] - cc[IDX3(2,3,k, ido_v,4)];
        ti2 = cc[IDX3(2,1,k, ido_v,4)] + cc[IDX3(2,3,k, ido_v,4)];
        tr4 = cc[IDX3(2,4,k, ido_v,4)] - cc[IDX3(2,2,k, ido_v,4)];
        ti3 = cc[IDX3(2,2,k, ido_v,4)] + cc[IDX3(2,4,k, ido_v,4)];
        tr1 = cc[IDX3(1,1,k, ido_v,4)] - cc[IDX3(1,3,k, ido_v,4)];
        tr2 = cc[IDX3(1,1,k, ido_v,4)] + cc[IDX3(1,3,k, ido_v,4)];
        ti4 = cc[IDX3(1,2,k, ido_v,4)] - cc[IDX3(1,4,k, ido_v,4)];
        tr3 = cc[IDX3(1,2,k, ido_v,4)] + cc[IDX3(1,4,k, ido_v,4)];
        ch[IDX3(1,k,1, ido_v,l1_v)] = tr2 + tr3;
        ch[IDX3(1,k,3, ido_v,l1_v)] = tr2 - tr3;
        ch[IDX3(2,k,1, ido_v,l1_v)] = ti2 + ti3;
        ch[IDX3(2,k,3, ido_v,l1_v)] = ti2 - ti3;
        ch[IDX3(1,k,2, ido_v,l1_v)] = tr1 + tr4;
        ch[IDX3(1,k,4, ido_v,l1_v)] = tr1 - tr4;
        ch[IDX3(2,k,2, ido_v,l1_v)] = ti1 + ti4;
        ch[IDX3(2,k,4, ido_v,l1_v)] = ti1 - ti4;
    }
    return;
L102:
    for (k = 1; k <= l1_v; k++) {
        for (i = 2; i <= ido_v; i += 2) {
            ti1 = cc[IDX3(i,1,k, ido_v,4)] - cc[IDX3(i,3,k, ido_v,4)];
            ti2 = cc[IDX3(i,1,k, ido_v,4)] + cc[IDX3(i,3,k, ido_v,4)];
            ti3 = cc[IDX3(i,2,k, ido_v,4)] + cc[IDX3(i,4,k, ido_v,4)];
            tr4 = cc[IDX3(i,4,k, ido_v,4)] - cc[IDX3(i,2,k, ido_v,4)];
            tr1 = cc[IDX3(i-1,1,k, ido_v,4)] - cc[IDX3(i-1,3,k, ido_v,4)];
            tr2 = cc[IDX3(i-1,1,k, ido_v,4)] + cc[IDX3(i-1,3,k, ido_v,4)];
            ti4 = cc[IDX3(i-1,2,k, ido_v,4)] - cc[IDX3(i-1,4,k, ido_v,4)];
            tr3 = cc[IDX3(i-1,2,k, ido_v,4)] + cc[IDX3(i-1,4,k, ido_v,4)];
            ch[IDX3(i-1,k,1, ido_v,l1_v)] = tr2 + tr3;
            cr3 = tr2 - tr3;
            ch[IDX3(i,k,1, ido_v,l1_v)] = ti2 + ti3;
            ci3 = ti2 - ti3;
            cr2 = tr1 + tr4;
            cr4 = tr1 - tr4;
            ci2 = ti1 + ti4;
            ci4 = ti1 - ti4;
            ch[IDX3(i-1,k,2, ido_v,l1_v)] = wa1[i-1-1]*cr2 - wa1[i-1]*ci2;
            ch[IDX3(i,k,2, ido_v,l1_v)] = wa1[i-1-1]*ci2 + wa1[i-1]*cr2;
            ch[IDX3(i-1,k,3, ido_v,l1_v)] = wa2[i-1-1]*cr3 - wa2[i-1]*ci3;
            ch[IDX3(i,k,3, ido_v,l1_v)] = wa2[i-1-1]*ci3 + wa2[i-1]*cr3;
            ch[IDX3(i-1,k,4, ido_v,l1_v)] = wa3[i-1-1]*cr4 - wa3[i-1]*ci4;
            ch[IDX3(i,k,4, ido_v,l1_v)] = wa3[i-1-1]*ci4 + wa3[i-1]*cr4;
        }
    }
}

/* ================================================================
 * dpassb5 - backward butterfly, radix 5
 * ================================================================ */
void FNAME(dpassb5)(const fint *ido_p, const fint *l1_p,
                    const double *cc, double *ch,
                    const double *wa1, const double *wa2,
                    const double *wa3, const double *wa4)
{
    fint ido_v = *ido_p;
    fint l1_v = *l1_p;
    fint k, i;
    double ti5, ti2, ti4, ti3, tr5, tr2, tr4, tr3;
    double cr2, ci2, cr3, ci3, cr5, ci5, cr4, ci4;
    double dr3, dr4, di3, di4, dr5, dr2, di5, di2;

    static const double tr11 =  0.30901699437494742410229341718281905;
    static const double ti11 =  0.95105651629515357211643933337938214;
    static const double tr12 = -0.80901699437494742410229341718281906;
    static const double ti12 =  0.58778525229247312916870595463907276;

    if (ido_v != 2) goto L102;
    for (k = 1; k <= l1_v; k++) {
        ti5 = cc[IDX3(2,2,k, ido_v,5)] - cc[IDX3(2,5,k, ido_v,5)];
        ti2 = cc[IDX3(2,2,k, ido_v,5)] + cc[IDX3(2,5,k, ido_v,5)];
        ti4 = cc[IDX3(2,3,k, ido_v,5)] - cc[IDX3(2,4,k, ido_v,5)];
        ti3 = cc[IDX3(2,3,k, ido_v,5)] + cc[IDX3(2,4,k, ido_v,5)];
        tr5 = cc[IDX3(1,2,k, ido_v,5)] - cc[IDX3(1,5,k, ido_v,5)];
        tr2 = cc[IDX3(1,2,k, ido_v,5)] + cc[IDX3(1,5,k, ido_v,5)];
        tr4 = cc[IDX3(1,3,k, ido_v,5)] - cc[IDX3(1,4,k, ido_v,5)];
        tr3 = cc[IDX3(1,3,k, ido_v,5)] + cc[IDX3(1,4,k, ido_v,5)];
        ch[IDX3(1,k,1, ido_v,l1_v)] = cc[IDX3(1,1,k, ido_v,5)] + tr2 + tr3;
        ch[IDX3(2,k,1, ido_v,l1_v)] = cc[IDX3(2,1,k, ido_v,5)] + ti2 + ti3;
        cr2 = cc[IDX3(1,1,k, ido_v,5)] + tr11*tr2 + tr12*tr3;
        ci2 = cc[IDX3(2,1,k, ido_v,5)] + tr11*ti2 + tr12*ti3;
        cr3 = cc[IDX3(1,1,k, ido_v,5)] + tr12*tr2 + tr11*tr3;
        ci3 = cc[IDX3(2,1,k, ido_v,5)] + tr12*ti2 + tr11*ti3;
        cr5 = ti11*tr5 + ti12*tr4;
        ci5 = ti11*ti5 + ti12*ti4;
        cr4 = ti12*tr5 - ti11*tr4;
        ci4 = ti12*ti5 - ti11*ti4;
        ch[IDX3(1,k,2, ido_v,l1_v)] = cr2 - ci5;
        ch[IDX3(1,k,5, ido_v,l1_v)] = cr2 + ci5;
        ch[IDX3(2,k,2, ido_v,l1_v)] = ci2 + cr5;
        ch[IDX3(2,k,3, ido_v,l1_v)] = ci3 + cr4;
        ch[IDX3(1,k,3, ido_v,l1_v)] = cr3 - ci4;
        ch[IDX3(1,k,4, ido_v,l1_v)] = cr3 + ci4;
        ch[IDX3(2,k,4, ido_v,l1_v)] = ci3 - cr4;
        ch[IDX3(2,k,5, ido_v,l1_v)] = ci2 - cr5;
    }
    return;
L102:
    for (k = 1; k <= l1_v; k++) {
        for (i = 2; i <= ido_v; i += 2) {
            ti5 = cc[IDX3(i,2,k, ido_v,5)] - cc[IDX3(i,5,k, ido_v,5)];
            ti2 = cc[IDX3(i,2,k, ido_v,5)] + cc[IDX3(i,5,k, ido_v,5)];
            ti4 = cc[IDX3(i,3,k, ido_v,5)] - cc[IDX3(i,4,k, ido_v,5)];
            ti3 = cc[IDX3(i,3,k, ido_v,5)] + cc[IDX3(i,4,k, ido_v,5)];
            tr5 = cc[IDX3(i-1,2,k, ido_v,5)] - cc[IDX3(i-1,5,k, ido_v,5)];
            tr2 = cc[IDX3(i-1,2,k, ido_v,5)] + cc[IDX3(i-1,5,k, ido_v,5)];
            tr4 = cc[IDX3(i-1,3,k, ido_v,5)] - cc[IDX3(i-1,4,k, ido_v,5)];
            tr3 = cc[IDX3(i-1,3,k, ido_v,5)] + cc[IDX3(i-1,4,k, ido_v,5)];
            ch[IDX3(i-1,k,1, ido_v,l1_v)] = cc[IDX3(i-1,1,k, ido_v,5)] + tr2 + tr3;
            ch[IDX3(i,k,1, ido_v,l1_v)] = cc[IDX3(i,1,k, ido_v,5)] + ti2 + ti3;
            cr2 = cc[IDX3(i-1,1,k, ido_v,5)] + tr11*tr2 + tr12*tr3;
            ci2 = cc[IDX3(i,1,k, ido_v,5)] + tr11*ti2 + tr12*ti3;
            cr3 = cc[IDX3(i-1,1,k, ido_v,5)] + tr12*tr2 + tr11*tr3;
            ci3 = cc[IDX3(i,1,k, ido_v,5)] + tr12*ti2 + tr11*ti3;
            cr5 = ti11*tr5 + ti12*tr4;
            ci5 = ti11*ti5 + ti12*ti4;
            cr4 = ti12*tr5 - ti11*tr4;
            ci4 = ti12*ti5 - ti11*ti4;
            dr3 = cr3 - ci4;
            dr4 = cr3 + ci4;
            di3 = ci3 + cr4;
            di4 = ci3 - cr4;
            dr5 = cr2 + ci5;
            dr2 = cr2 - ci5;
            di5 = ci2 - cr5;
            di2 = ci2 + cr5;
            ch[IDX3(i-1,k,2, ido_v,l1_v)] = wa1[i-1-1]*dr2 - wa1[i-1]*di2;
            ch[IDX3(i,k,2, ido_v,l1_v)] = wa1[i-1-1]*di2 + wa1[i-1]*dr2;
            ch[IDX3(i-1,k,3, ido_v,l1_v)] = wa2[i-1-1]*dr3 - wa2[i-1]*di3;
            ch[IDX3(i,k,3, ido_v,l1_v)] = wa2[i-1-1]*di3 + wa2[i-1]*dr3;
            ch[IDX3(i-1,k,4, ido_v,l1_v)] = wa3[i-1-1]*dr4 - wa3[i-1]*di4;
            ch[IDX3(i,k,4, ido_v,l1_v)] = wa3[i-1-1]*di4 + wa3[i-1]*dr4;
            ch[IDX3(i-1,k,5, ido_v,l1_v)] = wa4[i-1-1]*dr5 - wa4[i-1]*di5;
            ch[IDX3(i,k,5, ido_v,l1_v)] = wa4[i-1-1]*di5 + wa4[i-1]*dr5;
        }
    }
}

/* ================================================================
 * dpassb - backward butterfly, general radix
 *
 * Note: cc, c1, c2 may alias the same buffer (3 views).
 *       ch, ch2 may alias the same buffer (2 views).
 * ================================================================ */
void FNAME(dpassb)(fint *nac_p, const fint *ido_p, const fint *ip_p,
                   const fint *l1_p, const fint *idl1_p,
                   double *cc, double *c1, double *c2,
                   double *ch, double *ch2, const double *wa)
{
    fint ido_v = *ido_p;
    fint ip_v = *ip_p;
    fint l1_v = *l1_p;
    fint idl1_v = *idl1_p;
    fint idot, ipp2, ipph, idp;
    fint j, jc, k, i, l, lc, ik, idl, inc, idlj, idij, idj;
    double war, wai;

    idot = ido_v / 2;
    ipp2 = ip_v + 2;
    ipph = (ip_v + 1) / 2;
    idp = ip_v * ido_v;

    if (ido_v < l1_v) goto L106;
    for (j = 2; j <= ipph; j++) {
        jc = ipp2 - j;
        for (k = 1; k <= l1_v; k++) {
            for (i = 1; i <= ido_v; i++) {
                ch[IDX3(i,k,j, ido_v,l1_v)] = cc[IDX3(i,j,k, ido_v,ip_v)] + cc[IDX3(i,jc,k, ido_v,ip_v)];
                ch[IDX3(i,k,jc, ido_v,l1_v)] = cc[IDX3(i,j,k, ido_v,ip_v)] - cc[IDX3(i,jc,k, ido_v,ip_v)];
            }
        }
    }
    for (k = 1; k <= l1_v; k++) {
        for (i = 1; i <= ido_v; i++) {
            ch[IDX3(i,k,1, ido_v,l1_v)] = cc[IDX3(i,1,k, ido_v,ip_v)];
        }
    }
    goto L112;

L106:
    for (j = 2; j <= ipph; j++) {
        jc = ipp2 - j;
        for (i = 1; i <= ido_v; i++) {
            for (k = 1; k <= l1_v; k++) {
                ch[IDX3(i,k,j, ido_v,l1_v)] = cc[IDX3(i,j,k, ido_v,ip_v)] + cc[IDX3(i,jc,k, ido_v,ip_v)];
                ch[IDX3(i,k,jc, ido_v,l1_v)] = cc[IDX3(i,j,k, ido_v,ip_v)] - cc[IDX3(i,jc,k, ido_v,ip_v)];
            }
        }
    }
    for (i = 1; i <= ido_v; i++) {
        for (k = 1; k <= l1_v; k++) {
            ch[IDX3(i,k,1, ido_v,l1_v)] = cc[IDX3(i,1,k, ido_v,ip_v)];
        }
    }

L112:
    idl = 2 - ido_v;
    inc = 0;
    for (l = 2; l <= ipph; l++) {
        lc = ipp2 - l;
        idl = idl + ido_v;
        for (ik = 1; ik <= idl1_v; ik++) {
            c2[IDX2(ik,l, idl1_v)] = ch2[IDX2(ik,1, idl1_v)] + wa[idl-1-1]*ch2[IDX2(ik,2, idl1_v)];
            c2[IDX2(ik,lc, idl1_v)] = wa[idl-1]*ch2[IDX2(ik,ip_v, idl1_v)];
        }
        idlj = idl;
        inc = inc + ido_v;
        for (j = 3; j <= ipph; j++) {
            jc = ipp2 - j;
            idlj = idlj + inc;
            if (idlj > idp) idlj = idlj - idp;
            war = wa[idlj-1-1];
            wai = wa[idlj-1];
            for (ik = 1; ik <= idl1_v; ik++) {
                c2[IDX2(ik,l, idl1_v)] = c2[IDX2(ik,l, idl1_v)] + war*ch2[IDX2(ik,j, idl1_v)];
                c2[IDX2(ik,lc, idl1_v)] = c2[IDX2(ik,lc, idl1_v)] + wai*ch2[IDX2(ik,jc, idl1_v)];
            }
        }
    }
    for (j = 2; j <= ipph; j++) {
        for (ik = 1; ik <= idl1_v; ik++) {
            ch2[IDX2(ik,1, idl1_v)] = ch2[IDX2(ik,1, idl1_v)] + ch2[IDX2(ik,j, idl1_v)];
        }
    }
    for (j = 2; j <= ipph; j++) {
        jc = ipp2 - j;
        for (ik = 2; ik <= idl1_v; ik += 2) {
            ch2[IDX2(ik-1,j, idl1_v)] = c2[IDX2(ik-1,j, idl1_v)] - c2[IDX2(ik,jc, idl1_v)];
            ch2[IDX2(ik-1,jc, idl1_v)] = c2[IDX2(ik-1,j, idl1_v)] + c2[IDX2(ik,jc, idl1_v)];
            ch2[IDX2(ik,j, idl1_v)] = c2[IDX2(ik,j, idl1_v)] + c2[IDX2(ik-1,jc, idl1_v)];
            ch2[IDX2(ik,jc, idl1_v)] = c2[IDX2(ik,j, idl1_v)] - c2[IDX2(ik-1,jc, idl1_v)];
        }
    }
    *nac_p = 1;
    if (ido_v == 2) return;
    *nac_p = 0;
    for (ik = 1; ik <= idl1_v; ik++) {
        c2[IDX2(ik,1, idl1_v)] = ch2[IDX2(ik,1, idl1_v)];
    }
    for (j = 2; j <= ip_v; j++) {
        for (k = 1; k <= l1_v; k++) {
            c1[IDX3(1,k,j, ido_v,l1_v)] = ch[IDX3(1,k,j, ido_v,l1_v)];
            c1[IDX3(2,k,j, ido_v,l1_v)] = ch[IDX3(2,k,j, ido_v,l1_v)];
        }
    }
    if (idot > l1_v) goto L127;
    idij = 0;
    for (j = 2; j <= ip_v; j++) {
        idij = idij + 2;
        for (i = 4; i <= ido_v; i += 2) {
            idij = idij + 2;
            for (k = 1; k <= l1_v; k++) {
                c1[IDX3(i-1,k,j, ido_v,l1_v)] = wa[idij-1-1]*ch[IDX3(i-1,k,j, ido_v,l1_v)] - wa[idij-1]*ch[IDX3(i,k,j, ido_v,l1_v)];
                c1[IDX3(i,k,j, ido_v,l1_v)] = wa[idij-1-1]*ch[IDX3(i,k,j, ido_v,l1_v)] + wa[idij-1]*ch[IDX3(i-1,k,j, ido_v,l1_v)];
            }
        }
    }
    return;
L127:
    idj = 2 - ido_v;
    for (j = 2; j <= ip_v; j++) {
        idj = idj + ido_v;
        for (k = 1; k <= l1_v; k++) {
            idij = idj;
            for (i = 4; i <= ido_v; i += 2) {
                idij = idij + 2;
                c1[IDX3(i-1,k,j, ido_v,l1_v)] = wa[idij-1-1]*ch[IDX3(i-1,k,j, ido_v,l1_v)] - wa[idij-1]*ch[IDX3(i,k,j, ido_v,l1_v)];
                c1[IDX3(i,k,j, ido_v,l1_v)] = wa[idij-1-1]*ch[IDX3(i,k,j, ido_v,l1_v)] + wa[idij-1]*ch[IDX3(i-1,k,j, ido_v,l1_v)];
            }
        }
    }
}

/* ================================================================
 * dpassf2 - forward butterfly, radix 2
 * ================================================================ */
void FNAME(dpassf2)(const fint *ido_p, const fint *l1_p,
                    const double *cc, double *ch, const double *wa1)
{
    fint ido_v = *ido_p;
    fint l1_v = *l1_p;
    fint k, i;
    double tr2, ti2;

    if (ido_v > 2) goto L102;
    for (k = 1; k <= l1_v; k++) {
        ch[IDX3(1,k,1, ido_v,l1_v)] = cc[IDX3(1,1,k, ido_v,2)] + cc[IDX3(1,2,k, ido_v,2)];
        ch[IDX3(1,k,2, ido_v,l1_v)] = cc[IDX3(1,1,k, ido_v,2)] - cc[IDX3(1,2,k, ido_v,2)];
        ch[IDX3(2,k,1, ido_v,l1_v)] = cc[IDX3(2,1,k, ido_v,2)] + cc[IDX3(2,2,k, ido_v,2)];
        ch[IDX3(2,k,2, ido_v,l1_v)] = cc[IDX3(2,1,k, ido_v,2)] - cc[IDX3(2,2,k, ido_v,2)];
    }
    return;
L102:
    for (k = 1; k <= l1_v; k++) {
        for (i = 2; i <= ido_v; i += 2) {
            ch[IDX3(i-1,k,1, ido_v,l1_v)] = cc[IDX3(i-1,1,k, ido_v,2)] + cc[IDX3(i-1,2,k, ido_v,2)];
            tr2 = cc[IDX3(i-1,1,k, ido_v,2)] - cc[IDX3(i-1,2,k, ido_v,2)];
            ch[IDX3(i,k,1, ido_v,l1_v)] = cc[IDX3(i,1,k, ido_v,2)] + cc[IDX3(i,2,k, ido_v,2)];
            ti2 = cc[IDX3(i,1,k, ido_v,2)] - cc[IDX3(i,2,k, ido_v,2)];
            ch[IDX3(i,k,2, ido_v,l1_v)] = wa1[i-1-1]*ti2 - wa1[i-1]*tr2;
            ch[IDX3(i-1,k,2, ido_v,l1_v)] = wa1[i-1-1]*tr2 + wa1[i-1]*ti2;
        }
    }
}

/* ================================================================
 * dpassf3 - forward butterfly, radix 3
 * ================================================================ */
void FNAME(dpassf3)(const fint *ido_p, const fint *l1_p,
                    const double *cc, double *ch,
                    const double *wa1, const double *wa2)
{
    fint ido_v = *ido_p;
    fint l1_v = *l1_p;
    fint k, i;
    double tr2, cr2, ti2, ci2, cr3, ci3, dr2, dr3, di2, di3;

    static const double taur = -0.5;
    static const double taui = -0.86602540378443864676372317075293618;

    if (ido_v != 2) goto L102;
    for (k = 1; k <= l1_v; k++) {
        tr2 = cc[IDX3(1,2,k, ido_v,3)] + cc[IDX3(1,3,k, ido_v,3)];
        cr2 = cc[IDX3(1,1,k, ido_v,3)] + taur*tr2;
        ch[IDX3(1,k,1, ido_v,l1_v)] = cc[IDX3(1,1,k, ido_v,3)] + tr2;
        ti2 = cc[IDX3(2,2,k, ido_v,3)] + cc[IDX3(2,3,k, ido_v,3)];
        ci2 = cc[IDX3(2,1,k, ido_v,3)] + taur*ti2;
        ch[IDX3(2,k,1, ido_v,l1_v)] = cc[IDX3(2,1,k, ido_v,3)] + ti2;
        cr3 = taui*(cc[IDX3(1,2,k, ido_v,3)] - cc[IDX3(1,3,k, ido_v,3)]);
        ci3 = taui*(cc[IDX3(2,2,k, ido_v,3)] - cc[IDX3(2,3,k, ido_v,3)]);
        ch[IDX3(1,k,2, ido_v,l1_v)] = cr2 - ci3;
        ch[IDX3(1,k,3, ido_v,l1_v)] = cr2 + ci3;
        ch[IDX3(2,k,2, ido_v,l1_v)] = ci2 + cr3;
        ch[IDX3(2,k,3, ido_v,l1_v)] = ci2 - cr3;
    }
    return;
L102:
    for (k = 1; k <= l1_v; k++) {
        for (i = 2; i <= ido_v; i += 2) {
            tr2 = cc[IDX3(i-1,2,k, ido_v,3)] + cc[IDX3(i-1,3,k, ido_v,3)];
            cr2 = cc[IDX3(i-1,1,k, ido_v,3)] + taur*tr2;
            ch[IDX3(i-1,k,1, ido_v,l1_v)] = cc[IDX3(i-1,1,k, ido_v,3)] + tr2;
            ti2 = cc[IDX3(i,2,k, ido_v,3)] + cc[IDX3(i,3,k, ido_v,3)];
            ci2 = cc[IDX3(i,1,k, ido_v,3)] + taur*ti2;
            ch[IDX3(i,k,1, ido_v,l1_v)] = cc[IDX3(i,1,k, ido_v,3)] + ti2;
            cr3 = taui*(cc[IDX3(i-1,2,k, ido_v,3)] - cc[IDX3(i-1,3,k, ido_v,3)]);
            ci3 = taui*(cc[IDX3(i,2,k, ido_v,3)] - cc[IDX3(i,3,k, ido_v,3)]);
            dr2 = cr2 - ci3;
            dr3 = cr2 + ci3;
            di2 = ci2 + cr3;
            di3 = ci2 - cr3;
            ch[IDX3(i,k,2, ido_v,l1_v)] = wa1[i-1-1]*di2 - wa1[i-1]*dr2;
            ch[IDX3(i-1,k,2, ido_v,l1_v)] = wa1[i-1-1]*dr2 + wa1[i-1]*di2;
            ch[IDX3(i,k,3, ido_v,l1_v)] = wa2[i-1-1]*di3 - wa2[i-1]*dr3;
            ch[IDX3(i-1,k,3, ido_v,l1_v)] = wa2[i-1-1]*dr3 + wa2[i-1]*di3;
        }
    }
}

/* ================================================================
 * dpassf4 - forward butterfly, radix 4
 * ================================================================ */
void FNAME(dpassf4)(const fint *ido_p, const fint *l1_p,
                    const double *cc, double *ch,
                    const double *wa1, const double *wa2,
                    const double *wa3)
{
    fint ido_v = *ido_p;
    fint l1_v = *l1_p;
    fint k, i;
    double ti1, ti2, tr4, ti3, tr1, tr2, ti4, tr3;
    double cr2, cr3, cr4, ci2, ci3, ci4;

    if (ido_v != 2) goto L102;
    for (k = 1; k <= l1_v; k++) {
        ti1 = cc[IDX3(2,1,k, ido_v,4)] - cc[IDX3(2,3,k, ido_v,4)];
        ti2 = cc[IDX3(2,1,k, ido_v,4)] + cc[IDX3(2,3,k, ido_v,4)];
        tr4 = cc[IDX3(2,2,k, ido_v,4)] - cc[IDX3(2,4,k, ido_v,4)];
        ti3 = cc[IDX3(2,2,k, ido_v,4)] + cc[IDX3(2,4,k, ido_v,4)];
        tr1 = cc[IDX3(1,1,k, ido_v,4)] - cc[IDX3(1,3,k, ido_v,4)];
        tr2 = cc[IDX3(1,1,k, ido_v,4)] + cc[IDX3(1,3,k, ido_v,4)];
        ti4 = cc[IDX3(1,4,k, ido_v,4)] - cc[IDX3(1,2,k, ido_v,4)];
        tr3 = cc[IDX3(1,2,k, ido_v,4)] + cc[IDX3(1,4,k, ido_v,4)];
        ch[IDX3(1,k,1, ido_v,l1_v)] = tr2 + tr3;
        ch[IDX3(1,k,3, ido_v,l1_v)] = tr2 - tr3;
        ch[IDX3(2,k,1, ido_v,l1_v)] = ti2 + ti3;
        ch[IDX3(2,k,3, ido_v,l1_v)] = ti2 - ti3;
        ch[IDX3(1,k,2, ido_v,l1_v)] = tr1 + tr4;
        ch[IDX3(1,k,4, ido_v,l1_v)] = tr1 - tr4;
        ch[IDX3(2,k,2, ido_v,l1_v)] = ti1 + ti4;
        ch[IDX3(2,k,4, ido_v,l1_v)] = ti1 - ti4;
    }
    return;
L102:
    for (k = 1; k <= l1_v; k++) {
        for (i = 2; i <= ido_v; i += 2) {
            ti1 = cc[IDX3(i,1,k, ido_v,4)] - cc[IDX3(i,3,k, ido_v,4)];
            ti2 = cc[IDX3(i,1,k, ido_v,4)] + cc[IDX3(i,3,k, ido_v,4)];
            ti3 = cc[IDX3(i,2,k, ido_v,4)] + cc[IDX3(i,4,k, ido_v,4)];
            tr4 = cc[IDX3(i,2,k, ido_v,4)] - cc[IDX3(i,4,k, ido_v,4)];
            tr1 = cc[IDX3(i-1,1,k, ido_v,4)] - cc[IDX3(i-1,3,k, ido_v,4)];
            tr2 = cc[IDX3(i-1,1,k, ido_v,4)] + cc[IDX3(i-1,3,k, ido_v,4)];
            ti4 = cc[IDX3(i-1,4,k, ido_v,4)] - cc[IDX3(i-1,2,k, ido_v,4)];
            tr3 = cc[IDX3(i-1,2,k, ido_v,4)] + cc[IDX3(i-1,4,k, ido_v,4)];
            ch[IDX3(i-1,k,1, ido_v,l1_v)] = tr2 + tr3;
            cr3 = tr2 - tr3;
            ch[IDX3(i,k,1, ido_v,l1_v)] = ti2 + ti3;
            ci3 = ti2 - ti3;
            cr2 = tr1 + tr4;
            cr4 = tr1 - tr4;
            ci2 = ti1 + ti4;
            ci4 = ti1 - ti4;
            ch[IDX3(i-1,k,2, ido_v,l1_v)] = wa1[i-1-1]*cr2 + wa1[i-1]*ci2;
            ch[IDX3(i,k,2, ido_v,l1_v)] = wa1[i-1-1]*ci2 - wa1[i-1]*cr2;
            ch[IDX3(i-1,k,3, ido_v,l1_v)] = wa2[i-1-1]*cr3 + wa2[i-1]*ci3;
            ch[IDX3(i,k,3, ido_v,l1_v)] = wa2[i-1-1]*ci3 - wa2[i-1]*cr3;
            ch[IDX3(i-1,k,4, ido_v,l1_v)] = wa3[i-1-1]*cr4 + wa3[i-1]*ci4;
            ch[IDX3(i,k,4, ido_v,l1_v)] = wa3[i-1-1]*ci4 - wa3[i-1]*cr4;
        }
    }
}

/* ================================================================
 * dpassf5 - forward butterfly, radix 5
 * ================================================================ */
void FNAME(dpassf5)(const fint *ido_p, const fint *l1_p,
                    const double *cc, double *ch,
                    const double *wa1, const double *wa2,
                    const double *wa3, const double *wa4)
{
    fint ido_v = *ido_p;
    fint l1_v = *l1_p;
    fint k, i;
    double ti5, ti2, ti4, ti3, tr5, tr2, tr4, tr3;
    double cr2, ci2, cr3, ci3, cr5, ci5, cr4, ci4;
    double dr3, dr4, di3, di4, dr5, dr2, di5, di2;

    static const double tr11 =  0.30901699437494742410229341718281905;
    static const double ti11 = -0.95105651629515357211643933337938214;
    static const double tr12 = -0.80901699437494742410229341718281906;
    static const double ti12 = -0.58778525229247312916870595463907276;

    if (ido_v != 2) goto L102;
    for (k = 1; k <= l1_v; k++) {
        ti5 = cc[IDX3(2,2,k, ido_v,5)] - cc[IDX3(2,5,k, ido_v,5)];
        ti2 = cc[IDX3(2,2,k, ido_v,5)] + cc[IDX3(2,5,k, ido_v,5)];
        ti4 = cc[IDX3(2,3,k, ido_v,5)] - cc[IDX3(2,4,k, ido_v,5)];
        ti3 = cc[IDX3(2,3,k, ido_v,5)] + cc[IDX3(2,4,k, ido_v,5)];
        tr5 = cc[IDX3(1,2,k, ido_v,5)] - cc[IDX3(1,5,k, ido_v,5)];
        tr2 = cc[IDX3(1,2,k, ido_v,5)] + cc[IDX3(1,5,k, ido_v,5)];
        tr4 = cc[IDX3(1,3,k, ido_v,5)] - cc[IDX3(1,4,k, ido_v,5)];
        tr3 = cc[IDX3(1,3,k, ido_v,5)] + cc[IDX3(1,4,k, ido_v,5)];
        ch[IDX3(1,k,1, ido_v,l1_v)] = cc[IDX3(1,1,k, ido_v,5)] + tr2 + tr3;
        ch[IDX3(2,k,1, ido_v,l1_v)] = cc[IDX3(2,1,k, ido_v,5)] + ti2 + ti3;
        cr2 = cc[IDX3(1,1,k, ido_v,5)] + tr11*tr2 + tr12*tr3;
        ci2 = cc[IDX3(2,1,k, ido_v,5)] + tr11*ti2 + tr12*ti3;
        cr3 = cc[IDX3(1,1,k, ido_v,5)] + tr12*tr2 + tr11*tr3;
        ci3 = cc[IDX3(2,1,k, ido_v,5)] + tr12*ti2 + tr11*ti3;
        cr5 = ti11*tr5 + ti12*tr4;
        ci5 = ti11*ti5 + ti12*ti4;
        cr4 = ti12*tr5 - ti11*tr4;
        ci4 = ti12*ti5 - ti11*ti4;
        ch[IDX3(1,k,2, ido_v,l1_v)] = cr2 - ci5;
        ch[IDX3(1,k,5, ido_v,l1_v)] = cr2 + ci5;
        ch[IDX3(2,k,2, ido_v,l1_v)] = ci2 + cr5;
        ch[IDX3(2,k,3, ido_v,l1_v)] = ci3 + cr4;
        ch[IDX3(1,k,3, ido_v,l1_v)] = cr3 - ci4;
        ch[IDX3(1,k,4, ido_v,l1_v)] = cr3 + ci4;
        ch[IDX3(2,k,4, ido_v,l1_v)] = ci3 - cr4;
        ch[IDX3(2,k,5, ido_v,l1_v)] = ci2 - cr5;
    }
    return;
L102:
    for (k = 1; k <= l1_v; k++) {
        for (i = 2; i <= ido_v; i += 2) {
            ti5 = cc[IDX3(i,2,k, ido_v,5)] - cc[IDX3(i,5,k, ido_v,5)];
            ti2 = cc[IDX3(i,2,k, ido_v,5)] + cc[IDX3(i,5,k, ido_v,5)];
            ti4 = cc[IDX3(i,3,k, ido_v,5)] - cc[IDX3(i,4,k, ido_v,5)];
            ti3 = cc[IDX3(i,3,k, ido_v,5)] + cc[IDX3(i,4,k, ido_v,5)];
            tr5 = cc[IDX3(i-1,2,k, ido_v,5)] - cc[IDX3(i-1,5,k, ido_v,5)];
            tr2 = cc[IDX3(i-1,2,k, ido_v,5)] + cc[IDX3(i-1,5,k, ido_v,5)];
            tr4 = cc[IDX3(i-1,3,k, ido_v,5)] - cc[IDX3(i-1,4,k, ido_v,5)];
            tr3 = cc[IDX3(i-1,3,k, ido_v,5)] + cc[IDX3(i-1,4,k, ido_v,5)];
            ch[IDX3(i-1,k,1, ido_v,l1_v)] = cc[IDX3(i-1,1,k, ido_v,5)] + tr2 + tr3;
            ch[IDX3(i,k,1, ido_v,l1_v)] = cc[IDX3(i,1,k, ido_v,5)] + ti2 + ti3;
            cr2 = cc[IDX3(i-1,1,k, ido_v,5)] + tr11*tr2 + tr12*tr3;
            ci2 = cc[IDX3(i,1,k, ido_v,5)] + tr11*ti2 + tr12*ti3;
            cr3 = cc[IDX3(i-1,1,k, ido_v,5)] + tr12*tr2 + tr11*tr3;
            ci3 = cc[IDX3(i,1,k, ido_v,5)] + tr12*ti2 + tr11*ti3;
            cr5 = ti11*tr5 + ti12*tr4;
            ci5 = ti11*ti5 + ti12*ti4;
            cr4 = ti12*tr5 - ti11*tr4;
            ci4 = ti12*ti5 - ti11*ti4;
            dr3 = cr3 - ci4;
            dr4 = cr3 + ci4;
            di3 = ci3 + cr4;
            di4 = ci3 - cr4;
            dr5 = cr2 + ci5;
            dr2 = cr2 - ci5;
            di5 = ci2 - cr5;
            di2 = ci2 + cr5;
            ch[IDX3(i-1,k,2, ido_v,l1_v)] = wa1[i-1-1]*dr2 + wa1[i-1]*di2;
            ch[IDX3(i,k,2, ido_v,l1_v)] = wa1[i-1-1]*di2 - wa1[i-1]*dr2;
            ch[IDX3(i-1,k,3, ido_v,l1_v)] = wa2[i-1-1]*dr3 + wa2[i-1]*di3;
            ch[IDX3(i,k,3, ido_v,l1_v)] = wa2[i-1-1]*di3 - wa2[i-1]*dr3;
            ch[IDX3(i-1,k,4, ido_v,l1_v)] = wa3[i-1-1]*dr4 + wa3[i-1]*di4;
            ch[IDX3(i,k,4, ido_v,l1_v)] = wa3[i-1-1]*di4 - wa3[i-1]*dr4;
            ch[IDX3(i-1,k,5, ido_v,l1_v)] = wa4[i-1-1]*dr5 + wa4[i-1]*di5;
            ch[IDX3(i,k,5, ido_v,l1_v)] = wa4[i-1-1]*di5 - wa4[i-1]*dr5;
        }
    }
}

/* ================================================================
 * dpassf - forward butterfly, general radix
 *
 * Note: cc, c1, c2 may alias the same buffer (3 views).
 *       ch, ch2 may alias the same buffer (2 views).
 * ================================================================ */
void FNAME(dpassf)(fint *nac_p, const fint *ido_p, const fint *ip_p,
                   const fint *l1_p, const fint *idl1_p,
                   double *cc, double *c1, double *c2,
                   double *ch, double *ch2, const double *wa)
{
    fint ido_v = *ido_p;
    fint ip_v = *ip_p;
    fint l1_v = *l1_p;
    fint idl1_v = *idl1_p;
    fint idot, ipp2, ipph, idp;
    fint j, jc, k, i, l, lc, ik, idl, inc, idlj, idij, idj;
    double war, wai;

    idot = ido_v / 2;
    ipp2 = ip_v + 2;
    ipph = (ip_v + 1) / 2;
    idp = ip_v * ido_v;

    if (ido_v < l1_v) goto L106;
    for (j = 2; j <= ipph; j++) {
        jc = ipp2 - j;
        for (k = 1; k <= l1_v; k++) {
            for (i = 1; i <= ido_v; i++) {
                ch[IDX3(i,k,j, ido_v,l1_v)] = cc[IDX3(i,j,k, ido_v,ip_v)] + cc[IDX3(i,jc,k, ido_v,ip_v)];
                ch[IDX3(i,k,jc, ido_v,l1_v)] = cc[IDX3(i,j,k, ido_v,ip_v)] - cc[IDX3(i,jc,k, ido_v,ip_v)];
            }
        }
    }
    for (k = 1; k <= l1_v; k++) {
        for (i = 1; i <= ido_v; i++) {
            ch[IDX3(i,k,1, ido_v,l1_v)] = cc[IDX3(i,1,k, ido_v,ip_v)];
        }
    }
    goto L112;

L106:
    for (j = 2; j <= ipph; j++) {
        jc = ipp2 - j;
        for (i = 1; i <= ido_v; i++) {
            for (k = 1; k <= l1_v; k++) {
                ch[IDX3(i,k,j, ido_v,l1_v)] = cc[IDX3(i,j,k, ido_v,ip_v)] + cc[IDX3(i,jc,k, ido_v,ip_v)];
                ch[IDX3(i,k,jc, ido_v,l1_v)] = cc[IDX3(i,j,k, ido_v,ip_v)] - cc[IDX3(i,jc,k, ido_v,ip_v)];
            }
        }
    }
    for (i = 1; i <= ido_v; i++) {
        for (k = 1; k <= l1_v; k++) {
            ch[IDX3(i,k,1, ido_v,l1_v)] = cc[IDX3(i,1,k, ido_v,ip_v)];
        }
    }

L112:
    idl = 2 - ido_v;
    inc = 0;
    for (l = 2; l <= ipph; l++) {
        lc = ipp2 - l;
        idl = idl + ido_v;
        for (ik = 1; ik <= idl1_v; ik++) {
            c2[IDX2(ik,l, idl1_v)] = ch2[IDX2(ik,1, idl1_v)] + wa[idl-1-1]*ch2[IDX2(ik,2, idl1_v)];
            c2[IDX2(ik,lc, idl1_v)] = -wa[idl-1]*ch2[IDX2(ik,ip_v, idl1_v)];
        }
        idlj = idl;
        inc = inc + ido_v;
        for (j = 3; j <= ipph; j++) {
            jc = ipp2 - j;
            idlj = idlj + inc;
            if (idlj > idp) idlj = idlj - idp;
            war = wa[idlj-1-1];
            wai = wa[idlj-1];
            for (ik = 1; ik <= idl1_v; ik++) {
                c2[IDX2(ik,l, idl1_v)] = c2[IDX2(ik,l, idl1_v)] + war*ch2[IDX2(ik,j, idl1_v)];
                c2[IDX2(ik,lc, idl1_v)] = c2[IDX2(ik,lc, idl1_v)] - wai*ch2[IDX2(ik,jc, idl1_v)];
            }
        }
    }
    for (j = 2; j <= ipph; j++) {
        for (ik = 1; ik <= idl1_v; ik++) {
            ch2[IDX2(ik,1, idl1_v)] = ch2[IDX2(ik,1, idl1_v)] + ch2[IDX2(ik,j, idl1_v)];
        }
    }
    for (j = 2; j <= ipph; j++) {
        jc = ipp2 - j;
        for (ik = 2; ik <= idl1_v; ik += 2) {
            ch2[IDX2(ik-1,j, idl1_v)] = c2[IDX2(ik-1,j, idl1_v)] - c2[IDX2(ik,jc, idl1_v)];
            ch2[IDX2(ik-1,jc, idl1_v)] = c2[IDX2(ik-1,j, idl1_v)] + c2[IDX2(ik,jc, idl1_v)];
            ch2[IDX2(ik,j, idl1_v)] = c2[IDX2(ik,j, idl1_v)] + c2[IDX2(ik-1,jc, idl1_v)];
            ch2[IDX2(ik,jc, idl1_v)] = c2[IDX2(ik,j, idl1_v)] - c2[IDX2(ik-1,jc, idl1_v)];
        }
    }
    *nac_p = 1;
    if (ido_v == 2) return;
    *nac_p = 0;
    for (ik = 1; ik <= idl1_v; ik++) {
        c2[IDX2(ik,1, idl1_v)] = ch2[IDX2(ik,1, idl1_v)];
    }
    for (j = 2; j <= ip_v; j++) {
        for (k = 1; k <= l1_v; k++) {
            c1[IDX3(1,k,j, ido_v,l1_v)] = ch[IDX3(1,k,j, ido_v,l1_v)];
            c1[IDX3(2,k,j, ido_v,l1_v)] = ch[IDX3(2,k,j, ido_v,l1_v)];
        }
    }
    if (idot > l1_v) goto L127;
    idij = 0;
    for (j = 2; j <= ip_v; j++) {
        idij = idij + 2;
        for (i = 4; i <= ido_v; i += 2) {
            idij = idij + 2;
            for (k = 1; k <= l1_v; k++) {
                c1[IDX3(i-1,k,j, ido_v,l1_v)] = wa[idij-1-1]*ch[IDX3(i-1,k,j, ido_v,l1_v)] + wa[idij-1]*ch[IDX3(i,k,j, ido_v,l1_v)];
                c1[IDX3(i,k,j, ido_v,l1_v)] = wa[idij-1-1]*ch[IDX3(i,k,j, ido_v,l1_v)] - wa[idij-1]*ch[IDX3(i-1,k,j, ido_v,l1_v)];
            }
        }
    }
    return;
L127:
    idj = 2 - ido_v;
    for (j = 2; j <= ip_v; j++) {
        idj = idj + ido_v;
        for (k = 1; k <= l1_v; k++) {
            idij = idj;
            for (i = 4; i <= ido_v; i += 2) {
                idij = idij + 2;
                c1[IDX3(i-1,k,j, ido_v,l1_v)] = wa[idij-1-1]*ch[IDX3(i-1,k,j, ido_v,l1_v)] + wa[idij-1]*ch[IDX3(i,k,j, ido_v,l1_v)];
                c1[IDX3(i,k,j, ido_v,l1_v)] = wa[idij-1-1]*ch[IDX3(i,k,j, ido_v,l1_v)] - wa[idij-1]*ch[IDX3(i-1,k,j, ido_v,l1_v)];
            }
        }
    }
}

/* ================================================================
 * zfftb1 - backward FFT dispatch (called by zfftb)
 * ================================================================ */
void FNAME(zfftb1)(const fint *n_p, double *c, double *ch,
                   double *wa, fint *ifac)
{
    fint n_v = *n_p;
    fint nf, na, l1, iw, k1, ip, l2, ido, idot, idl1;
    fint ix2, ix3, ix4;
    fint nac, n2, i;

    nf = ifac[1];
    na = 0;
    l1 = 1;
    iw = 1;
    for (k1 = 1; k1 <= nf; k1++) {
        ip = ifac[k1 + 2 - 1];
        l2 = ip * l1;
        ido = n_v / l2;
        idot = ido + ido;
        idl1 = idot * l1;
        if (ip != 4) goto L103;
        ix2 = iw + idot;
        ix3 = ix2 + idot;
        if (na != 0) goto L101;
        FNAME(dpassb4)(&idot, &l1, c, ch, &wa[iw-1], &wa[ix2-1], &wa[ix3-1]);
        goto L102;
    L101:
        FNAME(dpassb4)(&idot, &l1, ch, c, &wa[iw-1], &wa[ix2-1], &wa[ix3-1]);
    L102:
        na = 1 - na;
        goto L115;
    L103:
        if (ip != 2) goto L106;
        if (na != 0) goto L104;
        FNAME(dpassb2)(&idot, &l1, c, ch, &wa[iw-1]);
        goto L105;
    L104:
        FNAME(dpassb2)(&idot, &l1, ch, c, &wa[iw-1]);
    L105:
        na = 1 - na;
        goto L115;
    L106:
        if (ip != 3) goto L109;
        ix2 = iw + idot;
        if (na != 0) goto L107;
        FNAME(dpassb3)(&idot, &l1, c, ch, &wa[iw-1], &wa[ix2-1]);
        goto L108;
    L107:
        FNAME(dpassb3)(&idot, &l1, ch, c, &wa[iw-1], &wa[ix2-1]);
    L108:
        na = 1 - na;
        goto L115;
    L109:
        if (ip != 5) goto L112;
        ix2 = iw + idot;
        ix3 = ix2 + idot;
        ix4 = ix3 + idot;
        if (na != 0) goto L110;
        FNAME(dpassb5)(&idot, &l1, c, ch, &wa[iw-1], &wa[ix2-1], &wa[ix3-1], &wa[ix4-1]);
        goto L111;
    L110:
        FNAME(dpassb5)(&idot, &l1, ch, c, &wa[iw-1], &wa[ix2-1], &wa[ix3-1], &wa[ix4-1]);
    L111:
        na = 1 - na;
        goto L115;
    L112:
        if (na != 0) goto L113;
        FNAME(dpassb)(&nac, &idot, &ip, &l1, &idl1, c, c, c, ch, ch, &wa[iw-1]);
        goto L114;
    L113:
        FNAME(dpassb)(&nac, &idot, &ip, &l1, &idl1, ch, ch, ch, c, c, &wa[iw-1]);
    L114:
        if (nac != 0) na = 1 - na;
    L115:
        l1 = l2;
        iw = iw + (ip - 1) * idot;
    }
    if (na == 0) return;
    n2 = n_v + n_v;
    for (i = 1; i <= n2; i++) {
        c[i-1] = ch[i-1];
    }
}

/* ================================================================
 * zfftb - backward (inverse) FFT of a complex periodic sequence
 *
 * Thread-safe: allocates local workspace.
 * ================================================================ */
void FNAME(zfftb)(const fint *n_p, double *c, double *wsave)
{
    fint n_v = *n_p;
    fint iw1, iw2, wsz, i;
    double *wsavep;

    if (n_v == 1) return;

    wsz = 4 * n_v + 100;
    wsavep = (double *)malloc((size_t)wsz * sizeof(double));
    for (i = 0; i < wsz; i++) {
        wsavep[i] = wsave[i];
    }

    iw1 = n_v + n_v + 1;
    iw2 = iw1 + n_v + n_v;
    FNAME(zfftb1)(&n_v, c, wsavep, &wsavep[iw1-1],
                  (fint *)((void *)&wsavep[iw2-1]));

    free(wsavep);
}

/* ================================================================
 * zfftf1 - forward FFT dispatch (called by zfftf)
 * ================================================================ */
void FNAME(zfftf1)(const fint *n_p, double *c, double *ch,
                   double *wa, fint *ifac)
{
    fint n_v = *n_p;
    fint nf, na, l1, iw, k1, ip, l2, ido, idot, idl1;
    fint ix2, ix3, ix4;
    fint nac, n2, i;

    nf = ifac[1];
    na = 0;
    l1 = 1;
    iw = 1;
    for (k1 = 1; k1 <= nf; k1++) {
        ip = ifac[k1 + 2 - 1];
        l2 = ip * l1;
        ido = n_v / l2;
        idot = ido + ido;
        idl1 = idot * l1;
        if (ip != 4) goto L103;
        ix2 = iw + idot;
        ix3 = ix2 + idot;
        if (na != 0) goto L101;
        FNAME(dpassf4)(&idot, &l1, c, ch, &wa[iw-1], &wa[ix2-1], &wa[ix3-1]);
        goto L102;
    L101:
        FNAME(dpassf4)(&idot, &l1, ch, c, &wa[iw-1], &wa[ix2-1], &wa[ix3-1]);
    L102:
        na = 1 - na;
        goto L115;
    L103:
        if (ip != 2) goto L106;
        if (na != 0) goto L104;
        FNAME(dpassf2)(&idot, &l1, c, ch, &wa[iw-1]);
        goto L105;
    L104:
        FNAME(dpassf2)(&idot, &l1, ch, c, &wa[iw-1]);
    L105:
        na = 1 - na;
        goto L115;
    L106:
        if (ip != 3) goto L109;
        ix2 = iw + idot;
        if (na != 0) goto L107;
        FNAME(dpassf3)(&idot, &l1, c, ch, &wa[iw-1], &wa[ix2-1]);
        goto L108;
    L107:
        FNAME(dpassf3)(&idot, &l1, ch, c, &wa[iw-1], &wa[ix2-1]);
    L108:
        na = 1 - na;
        goto L115;
    L109:
        if (ip != 5) goto L112;
        ix2 = iw + idot;
        ix3 = ix2 + idot;
        ix4 = ix3 + idot;
        if (na != 0) goto L110;
        FNAME(dpassf5)(&idot, &l1, c, ch, &wa[iw-1], &wa[ix2-1], &wa[ix3-1], &wa[ix4-1]);
        goto L111;
    L110:
        FNAME(dpassf5)(&idot, &l1, ch, c, &wa[iw-1], &wa[ix2-1], &wa[ix3-1], &wa[ix4-1]);
    L111:
        na = 1 - na;
        goto L115;
    L112:
        if (na != 0) goto L113;
        FNAME(dpassf)(&nac, &idot, &ip, &l1, &idl1, c, c, c, ch, ch, &wa[iw-1]);
        goto L114;
    L113:
        FNAME(dpassf)(&nac, &idot, &ip, &l1, &idl1, ch, ch, ch, c, c, &wa[iw-1]);
    L114:
        if (nac != 0) na = 1 - na;
    L115:
        l1 = l2;
        iw = iw + (ip - 1) * idot;
    }
    if (na == 0) return;
    n2 = n_v + n_v;
    for (i = 1; i <= n2; i++) {
        c[i-1] = ch[i-1];
    }
}

/* ================================================================
 * zfftf - forward FFT of a complex periodic sequence
 *
 * Thread-safe: allocates local workspace.
 * ================================================================ */
void FNAME(zfftf)(const fint *n_p, double *c, double *wsave)
{
    fint n_v = *n_p;
    fint iw1, iw2, wsz, i;
    double *wsavep;

    if (n_v == 1) return;

    wsz = 4 * n_v + 100;
    wsavep = (double *)malloc((size_t)wsz * sizeof(double));
    for (i = 0; i < wsz; i++) {
        wsavep[i] = wsave[i];
    }

    iw1 = n_v + n_v + 1;
    iw2 = iw1 + n_v + n_v;
    FNAME(zfftf1)(&n_v, c, wsavep, &wsavep[iw1-1],
                  (fint *)((void *)&wsavep[iw2-1]));

    free(wsavep);
}

/* ================================================================
 * zffti1 - initialize twiddle factors and factorization
 * ================================================================ */
void FNAME(zffti1)(const fint *n_p, double *wa, fint *ifac)
{
    fint n_v = *n_p;
    fint ntryh[4] = {3, 4, 2, 5};
    fint nl, nf, j, ntry, nq, nr, ib;
    fint i, l1, k1, ip, ld, l2, ido, idot, ipm, i1, ii;
    double tpi, argh, fi, argld, arg;

    nl = n_v;
    nf = 0;
    j = 0;
    ntry = 0;

L101:
    j = j + 1;
    if (j <= 4)
        ntry = ntryh[j-1];
    else
        ntry = ntry + 2;

L104:
    nq = nl / ntry;
    nr = nl - ntry * nq;
    if (nr != 0) goto L101;

    /* factor found */
    nf = nf + 1;
    ifac[nf + 2 - 1] = ntry;
    nl = nq;
    if (ntry != 2) goto L107;
    if (nf == 1) goto L107;
    for (i = 2; i <= nf; i++) {
        ib = nf - i + 2;
        ifac[ib + 2 - 1] = ifac[ib + 1 - 1];
    }
    ifac[3 - 1] = 2;

L107:
    if (nl != 1) goto L104;
    ifac[1 - 1] = n_v;
    ifac[2 - 1] = nf;

    tpi = 6.2831853071795864769252867665590057;
    argh = tpi / (double)n_v;
    i = 2;
    l1 = 1;
    for (k1 = 1; k1 <= nf; k1++) {
        ip = ifac[k1 + 2 - 1];
        ld = 0;
        l2 = l1 * ip;
        ido = n_v / l2;
        idot = ido + ido + 2;
        ipm = ip - 1;
        for (j = 1; j <= ipm; j++) {
            i1 = i;
            wa[i-1-1] = 1.0;
            wa[i-1] = 0.0;
            ld = ld + l1;
            fi = 0.0;
            argld = (double)ld * argh;
            for (ii = 4; ii <= idot; ii += 2) {
                i = i + 2;
                fi = fi + 1.0;
                arg = fi * argld;
                wa[i-1-1] = cos(arg);
                wa[i-1] = sin(arg);
            }
            if (ip > 5) {
                wa[i1-1-1] = wa[i-1-1];
                wa[i1-1] = wa[i-1];
            }
        }
        l1 = l2;
    }
}

/* ================================================================
 * zffti - initialize zfftf and zfftb
 * ================================================================ */
void FNAME(zffti)(const fint *n_p, double *wsave)
{
    fint n_v = *n_p;
    fint iw1, iw2;

    if (n_v == 1) return;
    iw1 = n_v + n_v + 1;
    iw2 = iw1 + n_v + n_v;
    FNAME(zffti1)(&n_v, &wsave[iw1-1],
                  (fint *)((void *)&wsave[iw2-1]));
}
