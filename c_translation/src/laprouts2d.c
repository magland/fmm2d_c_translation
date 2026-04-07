/*
 * laprouts2d.c - C translation of src/laplace/laprouts2d.f
 *
 * Multipole and local expansion form/eval/translate routines for the
 * 2D Laplace FMM (Cauchy form). The expansion convention follows the
 * Fortran library: Re(log(z)) = log(|z|) for the n=0 multipole/local
 * term, complex powers of z (or rscale/z) for higher orders. All
 * routines INCREMENT their output buffers; they do not overwrite.
 *
 * Floating-point operation order is preserved exactly from the
 * Fortran source so that bit-for-bit equality with the Fortran
 * reference is achievable when both are compiled at -O0. The most
 * important consequence is: Fortran's left-to-right associativity for
 * `x = x + a + b` evaluates `(x + a) + b`, whereas the C shorthand
 * `x += a + b` parses as `x + (a + b)`. Wherever the Fortran source
 * has a multi-term right-hand side, the translation splits it into
 * one `+=` per term so the rounding pattern matches.
 */

#include "laprouts2d.h"

/* mpole(ii, n) where ii is 1-based and n is 0-based, leading dim nd */
#define MIDX(ii, n, nd) ((n) * (nd) + ((ii) - 1))

/* carray(l, m) where both are 0-based, leading dim is ldc+1 */
#define CIDX(l, m, ldc) ((m) * ((ldc) + 1) + (l))


/*
 * l2dformmpc: form multipole expansion from charges
 *
 *   mpole_0 += sum_j charge_j
 *   mpole_n += sum_j -charge_j (z0/rscale)^n / n         (n >= 1)
 *
 * The Fortran routine populates a local zpow(0:nterms) array per
 * source: zpow(0) = -1, zpow(n) = zpow(n-1)*ztemp1, then divides by
 * n, then sets zpow(0) = 1 so that the final accumulation
 * mpole(ii,n) += charge(ii,j)*zpow(n) handles the n=0 case (which is
 * just a sum of charges).
 */
void FNAME(l2dformmpc)(const fint *nd, const double *rscale, const double *source,
                       const fint *ns, const fcomplex *charge, const double *center,
                       const fint *nterms, fcomplex *mpole)
{
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nterms_v = *nterms;
    double rscale_v = *rscale;
    fint j, n, ii;
    double zdiff1, zdiff2;
    fcomplex z0, ztemp1;
    fcomplex *zpow = (fcomplex *)malloc((nterms_v + 1) * sizeof(fcomplex));

    for (j = 1; j <= ns_v; j++) {
        zdiff1 = source[FA2(1, j, 2)] - center[0];
        zdiff2 = source[FA2(2, j, 2)] - center[1];

        z0 = zdiff1 + zdiff2 * I;
        ztemp1 = z0 / rscale_v;

        zpow[0] = -1.0;
        for (n = 1; n <= nterms_v; n++) {
            zpow[n] = zpow[n - 1] * ztemp1;
        }
        for (n = 1; n <= nterms_v; n++) {
            zpow[n] = zpow[n] / ((double)(n));
        }
        zpow[0] = 1.0;
        for (n = 0; n <= nterms_v; n++) {
            for (ii = 1; ii <= nd_v; ii++) {
                mpole[MIDX(ii, n, nd_v)] += charge[FA2(ii, j, nd_v)] * zpow[n];
            }
        }
    }

    free(zpow);
}


/*
 * l2dformmpd: form multipole expansion from dipoles
 *
 *   mpole_n += sum_j dipstr_j z0^(n-1) / rscale^n          (n >= 1)
 *
 * Fortran source declares zpow(nterms) — a 1-indexed array of length
 * nterms — and only writes elements 1..nterms (the n=0 multipole is
 * untouched). We mirror that with a length-nterms allocation and
 * shift the index by one.
 */
void FNAME(l2dformmpd)(const fint *nd, const double *rscale, const double *source,
                       const fint *ns, const fcomplex *dipstr, const double *center,
                       const fint *nterms, fcomplex *mpole)
{
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nterms_v = *nterms;
    double rscale_v = *rscale;
    fint j, n, ii;
    double zdiff1, zdiff2;
    fcomplex z0, ztemp1;
    /* zpow is 1-indexed in Fortran (length nterms). We allocate one
     * extra slot so that zpow[1..nterms] corresponds to Fortran
     * zpow(1..nterms). */
    fcomplex *zpow = (fcomplex *)malloc((nterms_v + 1) * sizeof(fcomplex));

    for (j = 1; j <= ns_v; j++) {
        zdiff1 = source[FA2(1, j, 2)] - center[0];
        zdiff2 = source[FA2(2, j, 2)] - center[1];
        z0 = zdiff1 + zdiff2 * I;
        ztemp1 = z0 / rscale_v;
        zpow[1] = 1.0 / rscale_v;
        for (n = 2; n <= nterms_v; n++) {
            zpow[n] = zpow[n - 1] * ztemp1;
        }
        for (n = 1; n <= nterms_v; n++) {
            for (ii = 1; ii <= nd_v; ii++) {
                mpole[MIDX(ii, n, nd_v)] += dipstr[FA2(ii, j, nd_v)] * zpow[n];
            }
        }
    }

    free(zpow);
}


/*
 * l2dformmpcd: form multipole expansion from charges + dipoles
 *
 *   mpole_0 += sum_j charge_j
 *   mpole_n += sum_j (dipstr_j z0^(n-1) - charge_j z0^n / n) / rscale^n
 *
 * Per the comment in cauchykernels2d.c, the Fortran statement
 *   mpole(ii,n) = mpole(ii,n) + dipstr*zpowd + charge*zpowc
 * evaluates left-to-right, i.e. (mpole + dipstr*zpowd) + charge*zpowc.
 * To preserve this rounding order in C we use two += statements.
 */
void FNAME(l2dformmpcd)(const fint *nd, const double *rscale, const double *source,
                        const fint *ns, const fcomplex *charge, const fcomplex *dipstr,
                        const double *center, const fint *nterms, fcomplex *mpole)
{
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nterms_v = *nterms;
    double rscale_v = *rscale;
    fint j, n, ii;
    double zdiff1, zdiff2;
    fcomplex z0, ztemp1;
    fcomplex *zpowd = (fcomplex *)malloc((nterms_v + 1) * sizeof(fcomplex));
    fcomplex *zpowc = (fcomplex *)malloc((nterms_v + 1) * sizeof(fcomplex));

    for (j = 1; j <= ns_v; j++) {
        zdiff1 = source[FA2(1, j, 2)] - center[0];
        zdiff2 = source[FA2(2, j, 2)] - center[1];
        z0 = zdiff1 + zdiff2 * I;
        ztemp1 = z0 / rscale_v;
        zpowd[0] = 0.0;
        zpowd[1] = 1.0 / rscale_v;
        zpowc[0] = 1.0;
        zpowc[1] = -ztemp1;
        for (n = 2; n <= nterms_v; n++) {
            zpowd[n] = zpowd[n - 1] * ztemp1;
            zpowc[n] = zpowc[n - 1] * ztemp1;
        }
        for (n = 1; n <= nterms_v; n++) {
            zpowc[n] = zpowc[n] / ((double)(n));
        }
        for (n = 0; n <= nterms_v; n++) {
            for (ii = 1; ii <= nd_v; ii++) {
                /* Fortran: mpole(ii,n) = mpole(ii,n)
                 *                       + dipstr(ii,j)*zpowd(n)
                 *                       + charge(ii,j)*zpowc(n)
                 * left-to-right: ((mpole + dipstr*zpowd) + charge*zpowc).
                 * Splitting into two += statements preserves that. */
                mpole[MIDX(ii, n, nd_v)] += dipstr[FA2(ii, j, nd_v)] * zpowd[n];
                mpole[MIDX(ii, n, nd_v)] += charge[FA2(ii, j, nd_v)] * zpowc[n];
            }
        }
    }

    free(zpowd);
    free(zpowc);
}


/*
 * l2dmpevalp: evaluate multipole expansion potential at targets
 *
 *   pot += mpole_0 log(|z|) + sum_{n>=1} mpole_n (rscale/z)^n
 *
 * The single-term Fortran update mpole(ii,n)*zpow(n) translates to
 * one += per element.
 */
void FNAME(l2dmpevalp)(const fint *nd, const double *rscale, const double *center,
                       const fcomplex *mpole, const fint *nterms,
                       const double *ztarg, const fint *ntarg, fcomplex *pot)
{
    fint nd_v = *nd;
    fint nterms_v = *nterms;
    fint ntarg_v = *ntarg;
    double rscale_v = *rscale;
    fint k, n, ii;
    double zdiff1, zdiff2;
    fcomplex z0, ztemp1;
    fcomplex *zpow = (fcomplex *)malloc((nterms_v + 1) * sizeof(fcomplex));

    for (k = 1; k <= ntarg_v; k++) {
        zdiff1 = ztarg[FA2(1, k, 2)] - center[0];
        zdiff2 = ztarg[FA2(2, k, 2)] - center[1];
        z0 = zdiff1 + zdiff2 * I;
        zpow[0] = log(cabs(z0));
        ztemp1 = rscale_v / z0;
        zpow[1] = ztemp1;
        for (n = 2; n <= nterms_v; n++) {
            zpow[n] = zpow[n - 1] * ztemp1;
        }
        for (n = 0; n <= nterms_v; n++) {
            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii, k, nd_v)] += mpole[MIDX(ii, n, nd_v)] * zpow[n];
            }
        }
    }

    free(zpow);
}


/*
 * l2dmpevalg: evaluate multipole potential + gradient at targets
 *
 *   pot  += mpole_0 log(|z|) + sum_{n>=1} mpole_n (rscale/z)^n
 *   grad += d/dz pot
 *
 * Two single-term Fortran statements -> two single += C statements.
 * zpow has length nterms+2 in Fortran (declared 0:nterms+3 but only
 * filled through n=nterms+1); we allocate the same.
 */
void FNAME(l2dmpevalg)(const fint *nd, const double *rscale, const double *center,
                       const fcomplex *mpole, const fint *nterms,
                       const double *ztarg, const fint *ntarg,
                       fcomplex *pot, fcomplex *grad)
{
    fint nd_v = *nd;
    fint nterms_v = *nterms;
    fint ntarg_v = *ntarg;
    double rscale_v = *rscale;
    fint k, n, ii;
    double rinv;
    double zdiff1, zdiff2;
    fcomplex z, ztemp1;
    fcomplex *zpow = (fcomplex *)malloc((nterms_v + 4) * sizeof(fcomplex));
    fcomplex *zpowg = (fcomplex *)malloc((nterms_v + 4) * sizeof(fcomplex));

    rinv = 1.0 / rscale_v;
    for (k = 1; k <= ntarg_v; k++) {
        zdiff1 = ztarg[FA2(1, k, 2)] - center[0];
        zdiff2 = ztarg[FA2(2, k, 2)] - center[1];
        z = zdiff1 + zdiff2 * I;
        zpow[0] = log(cabs(z));
        ztemp1 = rscale_v / z;
        zpow[1] = ztemp1;
        for (n = 2; n <= nterms_v + 1; n++) {
            zpow[n] = zpow[n - 1] * ztemp1;
        }
        zpowg[0] = zpow[1] * rinv;
        for (n = 1; n <= nterms_v; n++) {
            zpowg[n] = -zpow[n + 1] * rinv * n;
        }
        for (n = 0; n <= nterms_v; n++) {
            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii, k, nd_v)]  += mpole[MIDX(ii, n, nd_v)] * zpow[n];
                grad[FA2(ii, k, nd_v)] += mpole[MIDX(ii, n, nd_v)] * zpowg[n];
            }
        }
    }

    free(zpow);
    free(zpowg);
}


/*
 * l2dmpevalh: evaluate multipole potential + gradient + Hessian
 *
 *   pot  += mpole_0 log(|z|) + sum_{n>=1} mpole_n (rscale/z)^n
 *   grad += d/dz pot
 *   hess += d/dz grad
 */
void FNAME(l2dmpevalh)(const fint *nd, const double *rscale, const double *center,
                       const fcomplex *mpole, const fint *nterms,
                       const double *ztarg, const fint *ntarg,
                       fcomplex *pot, fcomplex *grad, fcomplex *hess)
{
    fint nd_v = *nd;
    fint nterms_v = *nterms;
    fint ntarg_v = *ntarg;
    double rscale_v = *rscale;
    fint k, n, ii;
    double rinv, rinv2;
    double zdiff1, zdiff2;
    fcomplex z, ztemp1;
    fcomplex *zpow  = (fcomplex *)malloc((nterms_v + 4) * sizeof(fcomplex));
    fcomplex *zpowg = (fcomplex *)malloc((nterms_v + 4) * sizeof(fcomplex));
    fcomplex *zpowh = (fcomplex *)malloc((nterms_v + 4) * sizeof(fcomplex));

    rinv = 1.0 / rscale_v;
    rinv2 = rinv * rinv;
    zpowg[0] = 0.0;
    zpowh[0] = 0.0;
    zpowh[1] = 0.0;

    for (k = 1; k <= ntarg_v; k++) {
        zdiff1 = ztarg[FA2(1, k, 2)] - center[0];
        zdiff2 = ztarg[FA2(2, k, 2)] - center[1];
        z = zdiff1 + zdiff2 * I;
        zpow[0] = log(cabs(z));
        ztemp1 = rscale_v / z;
        zpow[1] = ztemp1;
        for (n = 2; n <= nterms_v + 2; n++) {
            zpow[n] = zpow[n - 1] * ztemp1;
        }

        zpowg[0] = zpow[1] * rinv;
        for (n = 1; n <= nterms_v; n++) {
            zpowg[n] = -zpow[n + 1] * rinv * n;
        }

        zpowh[0] = -zpow[2] * rinv2;
        for (n = 1; n <= nterms_v; n++) {
            zpowh[n] = zpow[n + 2] * n * (n + 1) * rinv2;
        }

        for (n = 0; n <= nterms_v; n++) {
            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii, k, nd_v)]  += mpole[MIDX(ii, n, nd_v)] * zpow[n];
                grad[FA2(ii, k, nd_v)] += mpole[MIDX(ii, n, nd_v)] * zpowg[n];
                hess[FA2(ii, k, nd_v)] += mpole[MIDX(ii, n, nd_v)] * zpowh[n];
            }
        }
    }

    free(zpow);
    free(zpowg);
    free(zpowh);
}


/*
 * l2dformtac: form local expansion from charges
 *
 *   local_0 += sum_j charge_j log(|z0|)
 *   local_n -= sum_j charge_j (rscale/z0)^n / n            (n >= 1)
 *
 * Same zpow(0) trick as l2dformmpc: zpow(0) starts as -1 to seed the
 * recurrence, then is overwritten with log(|z0|) before the final
 * accumulation. (zpow(n) for n>=1 already carries the correct -1
 * sign from the seed.)
 */
void FNAME(l2dformtac)(const fint *nd, const double *rscale, const double *source,
                       const fint *ns, const fcomplex *charge, const double *center,
                       const fint *nterms, fcomplex *local)
{
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nterms_v = *nterms;
    double rscale_v = *rscale;
    fint j, n, ii;
    double zdiff1, zdiff2;
    fcomplex z0, ztemp1;
    fcomplex *zpow = (fcomplex *)malloc((nterms_v + 1) * sizeof(fcomplex));

    for (j = 1; j <= ns_v; j++) {
        zdiff1 = source[FA2(1, j, 2)] - center[0];
        zdiff2 = source[FA2(2, j, 2)] - center[1];

        z0 = zdiff1 + zdiff2 * I;
        ztemp1 = rscale_v / z0;
        zpow[0] = -1.0;
        for (n = 1; n <= nterms_v; n++) {
            zpow[n] = zpow[n - 1] * ztemp1;
        }
        for (n = 1; n <= nterms_v; n++) {
            zpow[n] = zpow[n] / ((double)(n));
        }
        zpow[0] = log(cabs(z0));
        for (n = 0; n <= nterms_v; n++) {
            for (ii = 1; ii <= nd_v; ii++) {
                local[MIDX(ii, n, nd_v)] += charge[FA2(ii, j, nd_v)] * zpow[n];
            }
        }
    }

    free(zpow);
}


/*
 * l2dformtad: form local expansion from dipoles
 *
 *   local_n -= sum_j dipstr_j rscale^n / z0^(n+1)
 *
 * The seed zpow(0) = -1/z0, recurrence zpow(n) = zpow(n-1)*ztemp1
 * with ztemp1 = rscale/z0.
 */
void FNAME(l2dformtad)(const fint *nd, const double *rscale, const double *source,
                       const fint *ns, const fcomplex *dipstr, const double *center,
                       const fint *nterms, fcomplex *local)
{
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nterms_v = *nterms;
    double rscale_v = *rscale;
    fint j, n, ii;
    double zdiff1, zdiff2;
    fcomplex z0, ztemp1;
    fcomplex *zpow = (fcomplex *)malloc((nterms_v + 1) * sizeof(fcomplex));

    for (j = 1; j <= ns_v; j++) {
        zdiff1 = source[FA2(1, j, 2)] - center[0];
        zdiff2 = source[FA2(2, j, 2)] - center[1];
        z0 = zdiff1 + zdiff2 * I;
        ztemp1 = rscale_v / z0;
        zpow[0] = -1.0 / z0;
        for (n = 1; n <= nterms_v; n++) {
            zpow[n] = zpow[n - 1] * ztemp1;
        }
        for (n = 0; n <= nterms_v; n++) {
            for (ii = 1; ii <= nd_v; ii++) {
                local[MIDX(ii, n, nd_v)] += dipstr[FA2(ii, j, nd_v)] * zpow[n];
            }
        }
    }

    free(zpow);
}


/*
 * l2dformtacd: form local expansion from charges + dipoles
 *
 *   local_0 += sum_j charge_j log(|z0|)
 *   local_n -= sum_j (dipstr_j/z0^(n+1) + charge_j/(n z0^n)) rscale^n
 *
 * Two-term Fortran update -> two += C statements (left-to-right
 * preservation: charge contribution first, then dipstr).
 *
 * Note: the Fortran source assigns ztemp1 = rscale/z0 a second time
 * before computing zpowd. The value is the same as the first
 * assignment, but we keep the same statements in the same order to
 * mirror the source.
 */
void FNAME(l2dformtacd)(const fint *nd, const double *rscale, const double *source,
                        const fint *ns, const fcomplex *charge, const fcomplex *dipstr,
                        const double *center, const fint *nterms, fcomplex *local)
{
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nterms_v = *nterms;
    double rscale_v = *rscale;
    fint j, n, ii;
    double zdiff1, zdiff2;
    fcomplex z0, ztemp1;
    fcomplex *zpow  = (fcomplex *)malloc((nterms_v + 1) * sizeof(fcomplex));
    fcomplex *zpowd = (fcomplex *)malloc((nterms_v + 1) * sizeof(fcomplex));

    for (j = 1; j <= ns_v; j++) {
        zdiff1 = source[FA2(1, j, 2)] - center[0];
        zdiff2 = source[FA2(2, j, 2)] - center[1];

        z0 = zdiff1 + zdiff2 * I;
        ztemp1 = rscale_v / z0;
        zpow[0] = -1.0;
        for (n = 1; n <= nterms_v; n++) {
            zpow[n] = zpow[n - 1] * ztemp1;
        }
        for (n = 1; n <= nterms_v; n++) {
            zpow[n] = zpow[n] / ((double)(n));
        }
        zpow[0] = log(cabs(z0));

        ztemp1 = rscale_v / z0;
        zpowd[0] = -1.0 / z0;
        for (n = 1; n <= nterms_v; n++) {
            zpowd[n] = zpowd[n - 1] * ztemp1;
        }

        for (n = 0; n <= nterms_v; n++) {
            for (ii = 1; ii <= nd_v; ii++) {
                /* Fortran: local(ii,n) = local(ii,n)
                 *                       + charge(ii,j)*zpow(n)
                 *                       + dipstr(ii,j)*zpowd(n)
                 * is ((local + charge*zpow) + dipstr*zpowd). */
                local[MIDX(ii, n, nd_v)] += charge[FA2(ii, j, nd_v)] * zpow[n];
                local[MIDX(ii, n, nd_v)] += dipstr[FA2(ii, j, nd_v)] * zpowd[n];
            }
        }
    }

    free(zpow);
    free(zpowd);
}


/*
 * l2dtaevalp: evaluate local expansion potential at targets
 *
 *   pot += sum_{n=0..nterms} local_n (z/rscale)^n
 */
void FNAME(l2dtaevalp)(const fint *nd, const double *rscale, const double *center,
                       const fcomplex *local, const fint *nterms,
                       const double *ztarg, const fint *ntarg, fcomplex *pot)
{
    fint nd_v = *nd;
    fint nterms_v = *nterms;
    fint ntarg_v = *ntarg;
    double rscale_v = *rscale;
    fint k, i, ii;
    double zdiff1, zdiff2;
    fcomplex z;
    fcomplex *zpow = (fcomplex *)malloc((nterms_v + 1) * sizeof(fcomplex));

    for (k = 1; k <= ntarg_v; k++) {
        zdiff1 = ztarg[FA2(1, k, 2)] - center[0];
        zdiff2 = ztarg[FA2(2, k, 2)] - center[1];
        z = (zdiff1 + zdiff2 * I) / rscale_v;
        zpow[0] = 1.0;
        for (i = 1; i <= nterms_v; i++) {
            zpow[i] = zpow[i - 1] * z;
        }

        for (i = 0; i <= nterms_v; i++) {
            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii, k, nd_v)] += local[MIDX(ii, i, nd_v)] * zpow[i];
            }
        }
    }

    free(zpow);
}


/*
 * l2dtaevalg: evaluate local expansion potential + gradient at targets
 *
 *   pot  += sum local_n (z/rscale)^n
 *   grad += sum local_n n (z/rscale)^(n-1) / rscale
 */
void FNAME(l2dtaevalg)(const fint *nd, const double *rscale, const double *center,
                       const fcomplex *local, const fint *nterms,
                       const double *ztarg, const fint *ntarg,
                       fcomplex *pot, fcomplex *grad)
{
    fint nd_v = *nd;
    fint nterms_v = *nterms;
    fint ntarg_v = *ntarg;
    double rscale_v = *rscale;
    fint k, n, ii;
    double rinv;
    double zdiff1, zdiff2;
    fcomplex z;
    fcomplex *zpow  = (fcomplex *)malloc((nterms_v + 1) * sizeof(fcomplex));
    fcomplex *zpowg = (fcomplex *)malloc((nterms_v + 1) * sizeof(fcomplex));

    rinv = 1.0 / rscale_v;
    for (k = 1; k <= ntarg_v; k++) {
        zdiff1 = ztarg[FA2(1, k, 2)] - center[0];
        zdiff2 = ztarg[FA2(2, k, 2)] - center[1];
        z = (zdiff1 + zdiff2 * I) / rscale_v;
        zpow[0] = 1.0;
        zpowg[0] = 0.0;
        for (n = 1; n <= nterms_v; n++) {
            zpow[n] = zpow[n - 1] * z;
        }
        for (n = 1; n <= nterms_v; n++) {
            zpowg[n] = zpow[n - 1] * n * rinv;
        }

        for (n = 0; n <= nterms_v; n++) {
            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii, k, nd_v)]  += local[MIDX(ii, n, nd_v)] * zpow[n];
                grad[FA2(ii, k, nd_v)] += local[MIDX(ii, n, nd_v)] * zpowg[n];
            }
        }
    }

    free(zpow);
    free(zpowg);
}


/*
 * l2dtaevalh: evaluate local expansion potential + gradient + Hessian
 */
void FNAME(l2dtaevalh)(const fint *nd, const double *rscale, const double *center,
                       const fcomplex *local, const fint *nterms,
                       const double *ztarg, const fint *ntarg,
                       fcomplex *pot, fcomplex *grad, fcomplex *hess)
{
    fint nd_v = *nd;
    fint nterms_v = *nterms;
    fint ntarg_v = *ntarg;
    double rscale_v = *rscale;
    fint k, n, ii;
    double rinv, rinv2;
    double zdiff1, zdiff2;
    fcomplex z;
    fcomplex *zpow  = (fcomplex *)malloc((nterms_v + 1) * sizeof(fcomplex));
    fcomplex *zpowg = (fcomplex *)malloc((nterms_v + 1) * sizeof(fcomplex));
    fcomplex *zpowh = (fcomplex *)malloc((nterms_v + 1) * sizeof(fcomplex));

    rinv = 1.0 / rscale_v;
    rinv2 = rinv * rinv;
    for (k = 1; k <= ntarg_v; k++) {
        zdiff1 = ztarg[FA2(1, k, 2)] - center[0];
        zdiff2 = ztarg[FA2(2, k, 2)] - center[1];
        z = (zdiff1 + zdiff2 * I) / rscale_v;
        zpow[0] = 1.0;
        zpowg[0] = 0.0;
        zpowh[0] = 0.0;
        zpowh[1] = 0.0;
        for (n = 1; n <= nterms_v; n++) {
            zpow[n] = zpow[n - 1] * z;
        }
        for (n = 1; n <= nterms_v; n++) {
            zpowg[n] = zpow[n - 1] * n * rinv;
        }
        for (n = 2; n <= nterms_v; n++) {
            zpowh[n] = zpow[n - 2] * n * (n - 1) * rinv2;
        }

        for (n = 0; n <= nterms_v; n++) {
            for (ii = 1; ii <= nd_v; ii++) {
                pot[FA2(ii, k, nd_v)]  += local[MIDX(ii, n, nd_v)] * zpow[n];
                grad[FA2(ii, k, nd_v)] += local[MIDX(ii, n, nd_v)] * zpowg[n];
                hess[FA2(ii, k, nd_v)] += local[MIDX(ii, n, nd_v)] * zpowh[n];
            }
        }
    }

    free(zpow);
    free(zpowg);
    free(zpowh);
}


/*
 * l2dmpmp: M2M translation. Shifts multipole expansion HEXP1 (about
 * CENTER1) to a new center CENTER2 and INCREMENTS HEXP2.
 *
 * The Fortran source allocates four temporary arrays
 * (z0pow1, z0pow2, hexp1tmp, hexp2tmp), zeros the two hexp temporaries,
 * scales hexp1 elementwise into hexp1tmp, builds hexp2tmp by
 * binomial-weighted accumulation, scales by z0pow2, and finally
 * increments hexp2.
 *
 * Note the order in `hexp2tmp(ii,i) = hexp2tmp(ii,i) - hexp1tmp(ii,0)/i`
 * which becomes a single -= statement (single-term right-hand side).
 *
 * The double-loop accumulation `hexp2tmp(ii,i) = hexp2tmp(ii,i) +
 * hexp1tmp(ii,j)*carray(i-1,j-1)` is also single-term so a single +=
 * is correct.
 *
 * Crucial subtlety: `ztemp1 = 1/(z0/rscale1)` is NOT the same as
 * `rscale1/z0` in floating point. We translate it literally as
 * `1.0 / (z0 / rscale1_v)` so the operation order matches.
 */
void FNAME(l2dmpmp)(const fint *nd,
                    const double *rscale1, const double *center1,
                    const fcomplex *hexp1, const fint *nterms1,
                    const double *rscale2, const double *center2,
                    fcomplex *hexp2, const fint *nterms2,
                    const double *carray, const fint *ldc)
{
    fint nd_v = *nd;
    fint nterms1_v = *nterms1;
    fint nterms2_v = *nterms2;
    fint ldc_v = *ldc;
    double rscale1_v = *rscale1;
    double rscale2_v = *rscale2;
    fint nmax;
    fint i, j, ii;
    double zdiff1, zdiff2;
    fcomplex z0, ztemp1, ztemp2;
    fcomplex *z0pow1;
    fcomplex *z0pow2;
    fcomplex *hexp1tmp;
    fcomplex *hexp2tmp;

    nmax = (nterms1_v > nterms2_v) ? nterms1_v : nterms2_v;

    zdiff1 = center2[0] - center1[0];
    zdiff2 = center2[1] - center1[1];
    z0 = -(zdiff1 + zdiff2 * I);

    z0pow1 = (fcomplex *)malloc((nmax + 1) * sizeof(fcomplex));
    ztemp1 = 1.0 / (z0 / rscale1_v);
    ztemp2 = ztemp1;
    z0pow1[0] = 1.0;
    for (i = 1; i <= nmax; i++) {
        z0pow1[i] = ztemp1;
        ztemp1 = ztemp1 * ztemp2;
    }

    z0pow2 = (fcomplex *)malloc((nmax + 1) * sizeof(fcomplex));
    ztemp1 = z0 / rscale2_v;
    ztemp2 = ztemp1;
    z0pow2[0] = 1.0;
    for (i = 1; i <= nmax; i++) {
        z0pow2[i] = ztemp1;
        ztemp1 = ztemp1 * ztemp2;
    }

    hexp1tmp = (fcomplex *)malloc(nd_v * (nterms1_v + 1) * sizeof(fcomplex));
    hexp2tmp = (fcomplex *)malloc(nd_v * (nterms2_v + 1) * sizeof(fcomplex));

    for (i = 0; i <= nterms1_v; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            hexp1tmp[MIDX(ii, i, nd_v)] = 0.0;
        }
    }

    for (i = 0; i <= nterms2_v; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            hexp2tmp[MIDX(ii, i, nd_v)] = 0.0;
        }
    }

    for (i = 0; i <= nterms1_v; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            hexp1tmp[MIDX(ii, i, nd_v)] = hexp1[MIDX(ii, i, nd_v)] * z0pow1[i];
        }
    }

    for (ii = 1; ii <= nd_v; ii++) {
        hexp2tmp[MIDX(ii, 0, nd_v)] += hexp1[MIDX(ii, 0, nd_v)];
    }

    for (i = 1; i <= nterms2_v; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            hexp2tmp[MIDX(ii, i, nd_v)] -= hexp1tmp[MIDX(ii, 0, nd_v)] / i;
        }
        {
            fint jmax = (i < nterms1_v) ? i : nterms1_v;
            for (j = 1; j <= jmax; j++) {
                for (ii = 1; ii <= nd_v; ii++) {
                    hexp2tmp[MIDX(ii, i, nd_v)] +=
                        hexp1tmp[MIDX(ii, j, nd_v)] * carray[CIDX(i - 1, j - 1, ldc_v)];
                }
            }
        }

        for (ii = 1; ii <= nd_v; ii++) {
            hexp2tmp[MIDX(ii, i, nd_v)] = hexp2tmp[MIDX(ii, i, nd_v)] * z0pow2[i];
        }
    }

    for (i = 0; i <= nterms2_v; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            hexp2[MIDX(ii, i, nd_v)] += hexp2tmp[MIDX(ii, i, nd_v)];
        }
    }

    free(z0pow1);
    free(z0pow2);
    free(hexp1tmp);
    free(hexp2tmp);
}


/*
 * l2dlocloc: L2L translation. Shifts local expansion JEXP1 to a new
 * center and INCREMENTS JEXP2.
 *
 * Note: in this routine z0 has the OPPOSITE sign convention from
 * l2dmpmp / l2dmploc — Fortran source has `z0 = dcmplx(zdiff(1),
 * zdiff(2))` (no leading minus). The translation matches.
 */
void FNAME(l2dlocloc)(const fint *nd,
                      const double *rscale1, const double *center1,
                      const fcomplex *jexp1, const fint *nterms1,
                      const double *rscale2, const double *center2,
                      fcomplex *jexp2, const fint *nterms2,
                      const double *carray, const fint *ldc)
{
    fint nd_v = *nd;
    fint nterms1_v = *nterms1;
    fint nterms2_v = *nterms2;
    fint ldc_v = *ldc;
    double rscale1_v = *rscale1;
    double rscale2_v = *rscale2;
    fint nmax;
    fint i, j, ii;
    double zdiff1, zdiff2;
    fcomplex z0, ztemp1, ztemp2;
    fcomplex *z0pow1;
    fcomplex *z0pow2;
    fcomplex *jexp1tmp;
    fcomplex *jexp2tmp;

    nmax = (nterms1_v > nterms2_v) ? nterms1_v : nterms2_v;

    zdiff1 = center2[0] - center1[0];
    zdiff2 = center2[1] - center1[1];
    z0 = zdiff1 + zdiff2 * I;

    z0pow1 = (fcomplex *)malloc((nmax + 1) * sizeof(fcomplex));
    ztemp1 = (z0 / rscale1_v);
    ztemp2 = ztemp1;
    z0pow1[0] = 1.0;
    for (i = 1; i <= nmax; i++) {
        z0pow1[i] = ztemp1;
        ztemp1 = ztemp1 * ztemp2;
    }

    z0pow2 = (fcomplex *)malloc((nmax + 1) * sizeof(fcomplex));
    ztemp1 = 1.0 / (z0 / rscale2_v);
    ztemp2 = ztemp1;
    z0pow2[0] = 1.0;
    for (i = 1; i <= nmax; i++) {
        z0pow2[i] = ztemp1;
        ztemp1 = ztemp1 * ztemp2;
    }

    jexp1tmp = (fcomplex *)malloc(nd_v * (nterms1_v + 1) * sizeof(fcomplex));
    jexp2tmp = (fcomplex *)malloc(nd_v * (nterms2_v + 1) * sizeof(fcomplex));

    for (i = 0; i <= nterms1_v; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            jexp1tmp[MIDX(ii, i, nd_v)] = 0.0;
        }
    }

    for (i = 0; i <= nterms2_v; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            jexp2tmp[MIDX(ii, i, nd_v)] = 0.0;
        }
    }

    for (i = 0; i <= nterms1_v; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            jexp1tmp[MIDX(ii, i, nd_v)] = jexp1[MIDX(ii, i, nd_v)] * z0pow1[i];
        }
    }

    for (i = 0; i <= nterms2_v; i++) {
        for (j = i; j <= nterms1_v; j++) {
            for (ii = 1; ii <= nd_v; ii++) {
                jexp2tmp[MIDX(ii, i, nd_v)] +=
                    jexp1tmp[MIDX(ii, j, nd_v)] * carray[CIDX(j, i, ldc_v)];
            }
        }

        for (ii = 1; ii <= nd_v; ii++) {
            jexp2tmp[MIDX(ii, i, nd_v)] = jexp2tmp[MIDX(ii, i, nd_v)] * z0pow2[i];
        }
    }

    for (i = 0; i <= nterms2_v; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            jexp2[MIDX(ii, i, nd_v)] += jexp2tmp[MIDX(ii, i, nd_v)];
        }
    }

    free(z0pow1);
    free(z0pow2);
    free(jexp1tmp);
    free(jexp2tmp);
}


/*
 * l2dmploc: M2L translation. Converts multipole expansion HEXP1
 * about CENTER1 to a local expansion JEXP2 about CENTER2 and
 * INCREMENTS JEXP2.
 *
 * Subtleties:
 *  - z0 = -(center2 - center1) here.
 *  - The Fortran builds z0pow1 and z0pow2 in a single loop with two
 *    parallel recurrences (rather than two separate loops). The
 *    translation matches the operation order exactly.
 *  - jexp2tmp(ii,0) is FIRST set to hexp1tmp(ii,0)*rtmp (NOT
 *    incremented — this overwrites the prior zero), then incremented
 *    by hexp1tmp(ii,j) for j=1..nterms1.
 *  - For i>=1: jexp2tmp(ii,i) -= hexp1tmp(ii,0)/i, then accumulate
 *    contributions from j=1..nterms1, then scale by z0pow2(i).
 */
void FNAME(l2dmploc)(const fint *nd,
                     const double *rscale1, const double *center1,
                     const fcomplex *hexp1, const fint *nterms1,
                     const double *rscale2, const double *center2,
                     fcomplex *jexp2, const fint *nterms2,
                     const double *carray, const fint *ldc)
{
    fint nd_v = *nd;
    fint nterms1_v = *nterms1;
    fint nterms2_v = *nterms2;
    fint ldc_v = *ldc;
    double rscale1_v = *rscale1;
    double rscale2_v = *rscale2;
    fint nmax;
    fint i, j, ii;
    double zdiff1, zdiff2;
    double rtmp;
    fcomplex z0, ztemp1, ztemp2, ztemp3;
    fcomplex *z0pow1;
    fcomplex *z0pow2;
    fcomplex *hexp1tmp;
    fcomplex *jexp2tmp;

    nmax = (nterms1_v > nterms2_v) ? nterms1_v : nterms2_v;

    zdiff1 = center2[0] - center1[0];
    zdiff2 = center2[1] - center1[1];
    z0 = -(zdiff1 + zdiff2 * I);

    z0pow1 = (fcomplex *)malloc((nmax + 1) * sizeof(fcomplex));
    z0pow2 = (fcomplex *)malloc((nmax + 1) * sizeof(fcomplex));
    ztemp1 = 1.0 / z0;
    z0pow1[0] = 1.0;
    z0pow2[0] = 1.0;
    ztemp2 = ztemp1 * rscale2_v;
    ztemp3 = -ztemp1 * rscale1_v;
    z0pow1[0] = 1.0;
    for (i = 1; i <= nmax; i++) {
        z0pow2[i] = ztemp2;
        z0pow1[i] = ztemp3;
        ztemp2 = ztemp2 * ztemp1 * rscale2_v;
        ztemp3 = -ztemp3 * ztemp1 * rscale1_v;
    }

    hexp1tmp = (fcomplex *)malloc(nd_v * (nterms1_v + 1) * sizeof(fcomplex));
    jexp2tmp = (fcomplex *)malloc(nd_v * (nterms2_v + 1) * sizeof(fcomplex));

    for (i = 0; i <= nterms1_v; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            hexp1tmp[MIDX(ii, i, nd_v)] = 0.0;
        }
    }

    for (i = 0; i <= nterms2_v; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            jexp2tmp[MIDX(ii, i, nd_v)] = 0.0;
        }
    }

    for (i = 0; i <= nterms1_v; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            hexp1tmp[MIDX(ii, i, nd_v)] = hexp1[MIDX(ii, i, nd_v)] * z0pow1[i];
        }
    }

    rtmp = log(cabs(z0));
    for (ii = 1; ii <= nd_v; ii++) {
        /* Fortran: jexp2tmp(ii,0) = hexp1tmp(ii,0)*rtmp  (assignment,
         * not increment). The prior zero in jexp2tmp is overwritten. */
        jexp2tmp[MIDX(ii, 0, nd_v)] = hexp1tmp[MIDX(ii, 0, nd_v)] * rtmp;
    }

    for (j = 1; j <= nterms1_v; j++) {
        for (ii = 1; ii <= nd_v; ii++) {
            jexp2tmp[MIDX(ii, 0, nd_v)] += hexp1tmp[MIDX(ii, j, nd_v)];
        }
    }

    for (i = 1; i <= nterms2_v; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            jexp2tmp[MIDX(ii, i, nd_v)] -= hexp1tmp[MIDX(ii, 0, nd_v)] / i;
        }
        for (j = 1; j <= nterms1_v; j++) {
            for (ii = 1; ii <= nd_v; ii++) {
                jexp2tmp[MIDX(ii, i, nd_v)] +=
                    hexp1tmp[MIDX(ii, j, nd_v)] * carray[CIDX(i + j - 1, j - 1, ldc_v)];
            }
        }

        for (ii = 1; ii <= nd_v; ii++) {
            jexp2tmp[MIDX(ii, i, nd_v)] = jexp2tmp[MIDX(ii, i, nd_v)] * z0pow2[i];
        }
    }

    for (i = 0; i <= nterms2_v; i++) {
        for (ii = 1; ii <= nd_v; ii++) {
            jexp2[MIDX(ii, i, nd_v)] += jexp2tmp[MIDX(ii, i, nd_v)];
        }
    }

    free(z0pow1);
    free(z0pow2);
    free(hexp1tmp);
    free(jexp2tmp);
}


/*
 * l2dmpzero: zero out a multipole array of shape (nd, 0:nterms).
 */
void FNAME(l2dmpzero)(const fint *nd, fcomplex *mpole, const fint *nterms)
{
    fint nd_v = *nd;
    fint nterms_v = *nterms;
    fint n, ii;

    for (n = 0; n <= nterms_v; n++) {
        for (ii = 1; ii <= nd_v; ii++) {
            mpole[MIDX(ii, n, nd_v)] = 0.0;
        }
    }
}
