/*
 * bhrouts2d.c - C translation of src/biharmonic/bhrouts2d.f
 *
 * Form / eval / translate routines for the 2D biharmonic (Stokes) FMM.
 * Each multipole or local expansion has FIVE sets of coefficients per
 * density slot, indexed in Fortran as `mpole(nd, 5, 0:nterms)`.
 *
 * Translated 1:1 from the Fortran reference. Multi-term Fortran
 * additions are split into separate `+=` statements (left-to-right) to
 * preserve operation order — required for the -O0 bit-for-bit diff
 * test.
 */

#include "bhrouts2d.h"

/*
 * Indexing helpers.
 *
 *   mpole(jj, kk, n)  with jj 1-based, kk 1..5 (1-based), n 0-based.
 *     Leading dimensions are (nd, 5).
 *
 *   h_a(jj, n)  with jj 1-based, n 0-based. Leading dim nd.
 *
 *   carray(l, m)  with both 0-based; leading dim is (ldc+1).
 *
 *   zpow(i)  is 1D; either 0-based (just zpow[i]) or 1-based
 *            (use zpow[i-1] -- handled at the call site).
 */
#define MPIDX(jj, kk, n, nd) (((n) * 5 + ((kk) - 1)) * (nd) + ((jj) - 1))
#define HIDX(jj, n, nd)      ((n) * (nd) + ((jj) - 1))
#define CIDX(l, m, ldc)      ((m) * ((ldc) + 1) + (l))


/* ================================================================ */
/* bh2dmpevalp - eval mpole expansion -> velocity                    */
/* ================================================================ */
void FNAME(bh2dmpevalp)(const fint *nd, const double *rscale,
                        const double *center, const fcomplex *mpole,
                        const fint *nterms, const double *ztarg,
                        const fint *ntarg, fcomplex *vel)
{
    fint nd_v = *nd;
    fint nterms_v = *nterms;
    fint ntarg_v = *ntarg;
    double rscale_v = *rscale;

    fint k, i, idim, nmax;
    fcomplex zc, zt, zdis, zdisinv, ztemp;
    fcomplex *zpow; /* zpow(1:nterms+5), 1-based */

    nmax = nterms_v + 3;
    zc = center[0] + center[1] * I;

    zpow = (fcomplex *)malloc((nterms_v + 5) * sizeof(fcomplex));

    for (k = 1; k <= ntarg_v; k++) {
        zt = ztarg[FA2(1, k, 2)] + ztarg[FA2(2, k, 2)] * I;
        zdis = zt - zc;
        zdisinv = 1.0 / zdis;
        ztemp = rscale_v * zdisinv;
        zpow[0] = ztemp;                          /* zpow(1) */
        for (i = 2; i <= nmax; i++) {
            zpow[i - 1] = zpow[i - 2] * ztemp;    /* zpow(i)=zpow(i-1)*ztemp */
        }

        /* log term */
        for (idim = 1; idim <= nd_v; idim++) {
            vel[FA2(idim, k, nd_v)] +=
                (mpole[MPIDX(idim, 4, 0, nd_v)]
                 + I * mpole[MPIDX(idim, 5, 0, nd_v)])
                * log(cabs(zdis));
        }

        for (i = 1; i <= nterms_v; i++) {
            for (idim = 1; idim <= nd_v; idim++) {
                /* vel = vel + mpole(1,i)*zpow(i) + mpole(2,i)*conj(zpow(i)) */
                vel[FA2(idim, k, nd_v)] +=
                    mpole[MPIDX(idim, 1, i, nd_v)] * zpow[i - 1];
                vel[FA2(idim, k, nd_v)] +=
                    mpole[MPIDX(idim, 2, i, nd_v)] * conj(zpow[i - 1]);
                /* vel = vel + mpole(3,i)*conj(zpow(i))*zdis */
                vel[FA2(idim, k, nd_v)] +=
                    mpole[MPIDX(idim, 3, i, nd_v)] * conj(zpow[i - 1]) * zdis;
                /* vel = vel + Re(mpole(4,i)*zpow(i)) + i*Re(mpole(5,i)*zpow(i)) */
                vel[FA2(idim, k, nd_v)] +=
                    creal(mpole[MPIDX(idim, 4, i, nd_v)] * zpow[i - 1]);
                vel[FA2(idim, k, nd_v)] +=
                    I * creal(mpole[MPIDX(idim, 5, i, nd_v)] * zpow[i - 1]);
            }
        }
    }

    free(zpow);
}


/* ================================================================ */
/* bh2dmpevalg - eval mpole expansion -> velocity + gradient         */
/* ================================================================ */
void FNAME(bh2dmpevalg)(const fint *nd, const double *rscale,
                        const double *center, const fcomplex *mpole,
                        const fint *nterms, const double *ztarg,
                        const fint *ntarg, fcomplex *vel, fcomplex *grad)
{
    fint nd_v = *nd;
    fint nterms_v = *nterms;
    fint ntarg_v = *ntarg;
    double rscale_v = *rscale;
    double rinv = 1.0 / rscale_v;

    fint k, i, idim, nmax;
    fcomplex zc, zt, zdis, zdisinv, ztemp;
    fcomplex *zpow;

    nmax = nterms_v + 3;
    zc = center[0] + center[1] * I;

    zpow = (fcomplex *)malloc((nterms_v + 5) * sizeof(fcomplex));

    for (k = 1; k <= ntarg_v; k++) {
        zt = ztarg[FA2(1, k, 2)] + ztarg[FA2(2, k, 2)] * I;
        zdis = zt - zc;
        zdisinv = 1.0 / zdis;
        ztemp = rscale_v * zdisinv;
        zpow[0] = ztemp;
        for (i = 2; i <= nmax; i++) {
            zpow[i - 1] = zpow[i - 2] * ztemp;
        }

        for (idim = 1; idim <= nd_v; idim++) {
            /* vel += (mp4(0)+i*mp5(0))*log|zdis| */
            vel[FA2(idim, k, nd_v)] +=
                (mpole[MPIDX(idim, 4, 0, nd_v)]
                 + I * mpole[MPIDX(idim, 5, 0, nd_v)])
                * log(cabs(zdis));
            /* grad(1) += 0.5*(mp4(0)+i*mp5(0))*zpow(1)*rinv */
            grad[FA3(idim, 1, k, nd_v, 3)] +=
                0.5 * (mpole[MPIDX(idim, 4, 0, nd_v)]
                       + I * mpole[MPIDX(idim, 5, 0, nd_v)])
                * zpow[0] * rinv;
            /* grad(3) += 0.5*(mp4(0)+i*mp5(0))*conj(zpow(1))*rinv */
            grad[FA3(idim, 3, k, nd_v, 3)] +=
                0.5 * (mpole[MPIDX(idim, 4, 0, nd_v)]
                       + I * mpole[MPIDX(idim, 5, 0, nd_v)])
                * conj(zpow[0]) * rinv;
        }

        for (i = 1; i <= nterms_v; i++) {
            for (idim = 1; idim <= nd_v; idim++) {
                /* vel = vel + mp(1,i)*zpow(i) + mp(2,i)*conj(zpow(i)) */
                vel[FA2(idim, k, nd_v)] +=
                    mpole[MPIDX(idim, 1, i, nd_v)] * zpow[i - 1];
                vel[FA2(idim, k, nd_v)] +=
                    mpole[MPIDX(idim, 2, i, nd_v)] * conj(zpow[i - 1]);
                /* vel = vel + mp(3,i)*conj(zpow(i))*zdis */
                vel[FA2(idim, k, nd_v)] +=
                    mpole[MPIDX(idim, 3, i, nd_v)] * conj(zpow[i - 1]) * zdis;
                /* vel = vel + Re(mp(4,i)*zpow(i)) + i*Re(mp(5,i)*zpow(i)) */
                vel[FA2(idim, k, nd_v)] +=
                    creal(mpole[MPIDX(idim, 4, i, nd_v)] * zpow[i - 1]);
                vel[FA2(idim, k, nd_v)] +=
                    I * creal(mpole[MPIDX(idim, 5, i, nd_v)] * zpow[i - 1]);

                /* grad(1) -= mp(1,i)*zpow(i+1)*i*rinv */
                grad[FA3(idim, 1, k, nd_v, 3)] -=
                    mpole[MPIDX(idim, 1, i, nd_v)] * zpow[i] * i * rinv;
                /* grad(1) -= 0.5*mp(4,i)*i*zpow(i+1)*rinv */
                grad[FA3(idim, 1, k, nd_v, 3)] -=
                    0.5 * mpole[MPIDX(idim, 4, i, nd_v)]
                    * i * zpow[i] * rinv;
                /* grad(1) -= eye*0.5*mp(5,i)*i*zpow(i+1)*rinv */
                grad[FA3(idim, 1, k, nd_v, 3)] -=
                    I * 0.5 * mpole[MPIDX(idim, 5, i, nd_v)]
                    * i * zpow[i] * rinv;

                /* grad(2) += mp(3,i)*conj(zpow(i)) */
                grad[FA3(idim, 2, k, nd_v, 3)] +=
                    mpole[MPIDX(idim, 3, i, nd_v)] * conj(zpow[i - 1]);

                /* grad(3) -= mp(2,i)*conj(zpow(i+1))*i*rinv */
                grad[FA3(idim, 3, k, nd_v, 3)] -=
                    mpole[MPIDX(idim, 2, i, nd_v)]
                    * conj(zpow[i]) * i * rinv;
                /* grad(3) -= zdis*mp(3,i)*conj(zpow(i+1))*i*rinv */
                grad[FA3(idim, 3, k, nd_v, 3)] -=
                    zdis * mpole[MPIDX(idim, 3, i, nd_v)]
                    * conj(zpow[i]) * i * rinv;
                /* grad(3) -= conj(0.5*mp(4,i)*i*zpow(i+1))*rinv */
                grad[FA3(idim, 3, k, nd_v, 3)] -=
                    conj(0.5 * mpole[MPIDX(idim, 4, i, nd_v)]
                         * i * zpow[i]) * rinv;
                /* grad(3) += conj(eye*0.5*mp(5,i)*i*zpow(i+1))*rinv */
                grad[FA3(idim, 3, k, nd_v, 3)] +=
                    conj(I * 0.5 * mpole[MPIDX(idim, 5, i, nd_v)]
                         * i * zpow[i]) * rinv;
            }
        }
    }

    free(zpow);
}


/* ================================================================ */
/* bh2dtaevalp - eval local expansion -> velocity                    */
/* ================================================================ */
void FNAME(bh2dtaevalp)(const fint *nd, const double *rscale,
                        const double *center, const fcomplex *mpole,
                        const fint *nterms, const double *ztarg,
                        const fint *ntarg, fcomplex *vel)
{
    fint nd_v = *nd;
    fint nterms_v = *nterms;
    fint ntarg_v = *ntarg;
    double rscale_v = *rscale;
    double rinv = 1.0 / rscale_v;

    fint k, i, idim;
    fcomplex zc, zt, zdis, ztemp1;
    fcomplex *zpow; /* zpow(0:nterms), 0-based */

    zc = center[0] + center[1] * I;
    zpow = (fcomplex *)malloc((nterms_v + 1) * sizeof(fcomplex));

    for (k = 1; k <= ntarg_v; k++) {
        zt = ztarg[FA2(1, k, 2)] + ztarg[FA2(2, k, 2)] * I;
        zdis = zt - zc;
        ztemp1 = zdis * rinv;
        zpow[0] = 1.0;
        for (i = 1; i <= nterms_v; i++) {
            zpow[i] = zpow[i - 1] * ztemp1;
        }

        for (i = 0; i <= nterms_v; i++) {
            for (idim = 1; idim <= nd_v; idim++) {
                /* vel = vel + mp(1,i)*zpow(i) + mp(2,i)*conj(zpow(i)) */
                vel[FA2(idim, k, nd_v)] +=
                    mpole[MPIDX(idim, 1, i, nd_v)] * zpow[i];
                vel[FA2(idim, k, nd_v)] +=
                    mpole[MPIDX(idim, 2, i, nd_v)] * conj(zpow[i]);
                /* vel = vel + mp(3,i)*conj(zpow(i))*zdis */
                vel[FA2(idim, k, nd_v)] +=
                    mpole[MPIDX(idim, 3, i, nd_v)] * conj(zpow[i]) * zdis;
                /* vel = vel + Re(mp(4,i)*zpow(i)) + i*Re(mp(5,i)*zpow(i)) */
                vel[FA2(idim, k, nd_v)] +=
                    creal(mpole[MPIDX(idim, 4, i, nd_v)] * zpow[i]);
                vel[FA2(idim, k, nd_v)] +=
                    I * creal(mpole[MPIDX(idim, 5, i, nd_v)] * zpow[i]);
            }
        }
    }

    free(zpow);
}


/* ================================================================ */
/* bh2dtaevalg - eval local expansion -> velocity + gradient         */
/* ================================================================ */
void FNAME(bh2dtaevalg)(const fint *nd, const double *rscale,
                        const double *center, const fcomplex *mpole,
                        const fint *nterms, const double *ztarg,
                        const fint *ntarg, fcomplex *vel, fcomplex *grad)
{
    fint nd_v = *nd;
    fint nterms_v = *nterms;
    fint ntarg_v = *ntarg;
    double rscale_v = *rscale;
    double rinv = 1.0 / rscale_v;

    fint k, i, idim;
    fcomplex zc, zt, zdis, ztemp1;
    fcomplex *zpow;

    zc = center[0] + center[1] * I;
    zpow = (fcomplex *)malloc((nterms_v + 1) * sizeof(fcomplex));

    for (k = 1; k <= ntarg_v; k++) {
        zt = ztarg[FA2(1, k, 2)] + ztarg[FA2(2, k, 2)] * I;
        zdis = zt - zc;
        ztemp1 = zdis * rinv;
        zpow[0] = 1.0;
        for (i = 1; i <= nterms_v; i++) {
            zpow[i] = zpow[i - 1] * ztemp1;
        }

        for (i = 0; i <= nterms_v; i++) {
            for (idim = 1; idim <= nd_v; idim++) {
                vel[FA2(idim, k, nd_v)] +=
                    mpole[MPIDX(idim, 1, i, nd_v)] * zpow[i];
                vel[FA2(idim, k, nd_v)] +=
                    mpole[MPIDX(idim, 2, i, nd_v)] * conj(zpow[i]);
                vel[FA2(idim, k, nd_v)] +=
                    mpole[MPIDX(idim, 3, i, nd_v)] * conj(zpow[i]) * zdis;
                vel[FA2(idim, k, nd_v)] +=
                    creal(mpole[MPIDX(idim, 4, i, nd_v)] * zpow[i]);
                vel[FA2(idim, k, nd_v)] +=
                    I * creal(mpole[MPIDX(idim, 5, i, nd_v)] * zpow[i]);
            }
        }

        /* grad(2,0) += mpole(3,0) (j=0 term, no zpow factor) */
        for (idim = 1; idim <= nd_v; idim++) {
            grad[FA3(idim, 2, k, nd_v, 3)] +=
                mpole[MPIDX(idim, 3, 0, nd_v)];
        }

        for (i = 1; i <= nterms_v; i++) {
            for (idim = 1; idim <= nd_v; idim++) {
                /* grad(1) += mp(1,i)*zpow(i-1)*i*rinv */
                grad[FA3(idim, 1, k, nd_v, 3)] +=
                    mpole[MPIDX(idim, 1, i, nd_v)]
                    * zpow[i - 1] * i * rinv;
                /* grad(1) += 0.5*mp(4,i)*i*zpow(i-1)*rinv */
                grad[FA3(idim, 1, k, nd_v, 3)] +=
                    0.5 * mpole[MPIDX(idim, 4, i, nd_v)]
                    * i * zpow[i - 1] * rinv;
                /* grad(1) += eye*0.5*mp(5,i)*i*zpow(i-1)*rinv */
                grad[FA3(idim, 1, k, nd_v, 3)] +=
                    I * 0.5 * mpole[MPIDX(idim, 5, i, nd_v)]
                    * i * zpow[i - 1] * rinv;

                /* grad(2) += mp(3,i)*conj(zpow(i)) */
                grad[FA3(idim, 2, k, nd_v, 3)] +=
                    mpole[MPIDX(idim, 3, i, nd_v)] * conj(zpow[i]);

                /* grad(3) += mp(2,i)*conj(zpow(i-1))*i*rinv */
                grad[FA3(idim, 3, k, nd_v, 3)] +=
                    mpole[MPIDX(idim, 2, i, nd_v)]
                    * conj(zpow[i - 1]) * i * rinv;
                /* grad(3) += mp(3,i)*conj(zpow(i-1))*i*zdis*rinv */
                grad[FA3(idim, 3, k, nd_v, 3)] +=
                    mpole[MPIDX(idim, 3, i, nd_v)]
                    * conj(zpow[i - 1]) * i * zdis * rinv;
                /* grad(3) += conj(0.5*mp(4,i)*i*zpow(i-1))*rinv */
                grad[FA3(idim, 3, k, nd_v, 3)] +=
                    conj(0.5 * mpole[MPIDX(idim, 4, i, nd_v)]
                         * i * zpow[i - 1]) * rinv;
                /* grad(3) -= conj(eye*0.5*mp(5,i)*i*zpow(i-1))*rinv */
                grad[FA3(idim, 3, k, nd_v, 3)] -=
                    conj(I * 0.5 * mpole[MPIDX(idim, 5, i, nd_v)]
                         * i * zpow[i - 1]) * rinv;
            }
        }
    }

    free(zpow);
}


/* ================================================================ */
/* bh2dformmpd - form mpole expansion from dipoles                   */
/* ================================================================ */
void FNAME(bh2dformmpd)(const fint *nd, const double *rscale,
                        const double *sources, const fint *ns,
                        const fcomplex *dip, const double *center,
                        const fint *nterms, fcomplex *mpole)
{
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nterms_v = *nterms;
    double rscale_v = *rscale;
    double rinv = 1.0 / rscale_v;

    fint i, j, idim;
    fcomplex zc, zs, zdis, ztemp, zdisc, ztempc, ztempc2;
    fcomplex zt2, zt3;

    zc = center[0] + center[1] * I;

    for (i = 1; i <= ns_v; i++) {
        zs = sources[FA2(1, i, 2)] + sources[FA2(2, i, 2)] * I;
        zdis = (zs - zc) / rscale_v;
        ztemp = zs - zc;
        zdisc = conj(zdis);

        if (cabs(zdis) <= 1.0e-16) {
            for (idim = 1; idim <= nd_v; idim++) {
                mpole[MPIDX(idim, 1, 1, nd_v)] +=
                    dip[FA3(idim, 1, i, nd_v, 3)] * rinv;
                mpole[MPIDX(idim, 2, 1, nd_v)] +=
                    dip[FA3(idim, 3, i, nd_v, 3)] * rinv;
                mpole[MPIDX(idim, 3, 2, nd_v)] +=
                    dip[FA3(idim, 2, i, nd_v, 3)] * rinv * rinv;
            }
        }

        if (cabs(zdis) > 1.0e-16) {
            ztempc = 1.0 / conj(ztemp);
            ztempc2 = ztempc * ztempc;
            for (j = 1; j <= nterms_v; j++) {
                for (idim = 1; idim <= nd_v; idim++) {
                    zt2 = dip[FA3(idim, 2, i, nd_v, 3)] * ztempc2;
                    zt3 = dip[FA3(idim, 3, i, nd_v, 3)] * ztempc;
                    /* mp(1,j) += dip(1,i)*zdis/ztemp */
                    mpole[MPIDX(idim, 1, j, nd_v)] +=
                        dip[FA3(idim, 1, i, nd_v, 3)] * zdis / ztemp;
                    /* mp(2,j) -= zt2*(j-1)*zdisc*ztemp */
                    mpole[MPIDX(idim, 2, j, nd_v)] -=
                        zt2 * (j - 1) * zdisc * ztemp;
                    /* mp(3,j) += zt2*(j-1)*zdisc */
                    mpole[MPIDX(idim, 3, j, nd_v)] +=
                        zt2 * (j - 1) * zdisc;
                    /* mp(2,j) += zt3*zdisc */
                    mpole[MPIDX(idim, 2, j, nd_v)] +=
                        zt3 * zdisc;
                }
                zdis = zdis * ztemp * rinv;
                zdisc = zdisc / ztempc * rinv;
            }
        }
    }
}


/* ================================================================ */
/* bh2dformmpc - form mpole expansion from charges                   */
/* ================================================================ */
void FNAME(bh2dformmpc)(const fint *nd, const double *rscale,
                        const double *sources, const fint *ns,
                        const fcomplex *c1, const double *center,
                        const fint *nterms, fcomplex *mpole)
{
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nterms_v = *nterms;
    double rscale_v = *rscale;
    double rinv = 1.0 / rscale_v;

    fint i, j, idim;
    fcomplex zc, zs, zdis, ztemp, zdisc, ztempc;
    fcomplex zt1;

    zc = center[0] + center[1] * I;

    for (i = 1; i <= ns_v; i++) {
        zs = sources[FA2(1, i, 2)] + sources[FA2(2, i, 2)] * I;
        zdis = (zs - zc) / rscale_v;
        ztemp = zs - zc;
        zdisc = conj(zdis);

        if (cabs(zdis) <= 1.0e-16) {
            for (idim = 1; idim <= nd_v; idim++) {
                mpole[MPIDX(idim, 4, 0, nd_v)] +=
                    creal(2.0 * c1[FA3(idim, 1, i, nd_v, 2)]);
                mpole[MPIDX(idim, 5, 0, nd_v)] +=
                    cimag(2.0 * c1[FA3(idim, 1, i, nd_v, 2)]);
                mpole[MPIDX(idim, 3, 1, nd_v)] +=
                    c1[FA3(idim, 2, i, nd_v, 2)] * rinv;
            }
        }

        if (cabs(zdis) > 1.0e-16) {
            ztempc = 1.0 / conj(ztemp);
            for (idim = 1; idim <= nd_v; idim++) {
                mpole[MPIDX(idim, 4, 0, nd_v)] +=
                    creal(2.0 * c1[FA3(idim, 1, i, nd_v, 2)]);
                mpole[MPIDX(idim, 5, 0, nd_v)] +=
                    cimag(2.0 * c1[FA3(idim, 1, i, nd_v, 2)]);
            }
            for (j = 1; j <= nterms_v; j++) {
                for (idim = 1; idim <= nd_v; idim++) {
                    zt1 = c1[FA3(idim, 2, i, nd_v, 2)] * ztempc;
                    /* mp(4,j) -= Re(2*c1(1,i))*zdis/j */
                    mpole[MPIDX(idim, 4, j, nd_v)] -=
                        creal(2.0 * c1[FA3(idim, 1, i, nd_v, 2)])
                        * zdis / j;
                    /* mp(5,j) -= Im(2*c1(1,i))*zdis/j */
                    mpole[MPIDX(idim, 5, j, nd_v)] -=
                        cimag(2.0 * c1[FA3(idim, 1, i, nd_v, 2)])
                        * zdis / j;
                    /* mp(2,j) -= zt1*zdisc*ztemp */
                    mpole[MPIDX(idim, 2, j, nd_v)] -=
                        zt1 * zdisc * ztemp;
                    /* mp(3,j) += zt1*zdisc */
                    mpole[MPIDX(idim, 3, j, nd_v)] +=
                        zt1 * zdisc;
                }
                zdis = zdis * ztemp * rinv;
                zdisc = zdisc / ztempc * rinv;
            }
        }
    }
}


/* ================================================================ */
/* bh2dformmpcd - form mpole expansion from charges + dipoles        */
/* ================================================================ */
void FNAME(bh2dformmpcd)(const fint *nd, const double *rscale,
                         const double *sources, const fint *ns,
                         const fcomplex *c1, const fcomplex *dip,
                         const double *center, const fint *nterms,
                         fcomplex *mpole)
{
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nterms_v = *nterms;
    double rscale_v = *rscale;
    double rinv = 1.0 / rscale_v;

    fint i, j, idim;
    fcomplex zc, zs, zdis, ztemp, zdisc, ztempc, ztempc2;
    fcomplex zt1, zt2, zt3;

    zc = center[0] + center[1] * I;

    for (i = 1; i <= ns_v; i++) {
        zs = sources[FA2(1, i, 2)] + sources[FA2(2, i, 2)] * I;
        zdis = (zs - zc) / rscale_v;
        ztemp = zs - zc;
        zdisc = conj(zdis);

        if (cabs(zdis) <= 1.0e-16) {
            for (idim = 1; idim <= nd_v; idim++) {
                mpole[MPIDX(idim, 4, 0, nd_v)] +=
                    creal(2.0 * c1[FA3(idim, 1, i, nd_v, 2)]);
                mpole[MPIDX(idim, 5, 0, nd_v)] +=
                    cimag(2.0 * c1[FA3(idim, 1, i, nd_v, 2)]);
                mpole[MPIDX(idim, 3, 1, nd_v)] +=
                    c1[FA3(idim, 2, i, nd_v, 2)] * rinv;
                mpole[MPIDX(idim, 1, 1, nd_v)] +=
                    dip[FA3(idim, 1, i, nd_v, 3)] * rinv;
                mpole[MPIDX(idim, 2, 1, nd_v)] +=
                    dip[FA3(idim, 3, i, nd_v, 3)] * rinv;
                mpole[MPIDX(idim, 3, 2, nd_v)] +=
                    dip[FA3(idim, 2, i, nd_v, 3)] * rinv * rinv;
            }
        }

        if (cabs(zdis) > 1.0e-16) {
            ztempc = 1.0 / conj(ztemp);
            ztempc2 = ztempc * ztempc;
            for (idim = 1; idim <= nd_v; idim++) {
                mpole[MPIDX(idim, 4, 0, nd_v)] +=
                    creal(2.0 * c1[FA3(idim, 1, i, nd_v, 2)]);
                mpole[MPIDX(idim, 5, 0, nd_v)] +=
                    cimag(2.0 * c1[FA3(idim, 1, i, nd_v, 2)]);
            }
            for (j = 1; j <= nterms_v; j++) {
                /* charge contribution */
                for (idim = 1; idim <= nd_v; idim++) {
                    zt1 = c1[FA3(idim, 2, i, nd_v, 2)] * ztempc;
                    mpole[MPIDX(idim, 4, j, nd_v)] -=
                        creal(2.0 * c1[FA3(idim, 1, i, nd_v, 2)])
                        * zdis / j;
                    mpole[MPIDX(idim, 5, j, nd_v)] -=
                        cimag(2.0 * c1[FA3(idim, 1, i, nd_v, 2)])
                        * zdis / j;
                    mpole[MPIDX(idim, 2, j, nd_v)] -=
                        zt1 * zdisc * ztemp;
                    mpole[MPIDX(idim, 3, j, nd_v)] +=
                        zt1 * zdisc;
                }
                /* dipole contribution */
                for (idim = 1; idim <= nd_v; idim++) {
                    zt2 = dip[FA3(idim, 2, i, nd_v, 3)] * ztempc2;
                    zt3 = dip[FA3(idim, 3, i, nd_v, 3)] * ztempc;
                    mpole[MPIDX(idim, 1, j, nd_v)] +=
                        dip[FA3(idim, 1, i, nd_v, 3)] * zdis / ztemp;
                    mpole[MPIDX(idim, 2, j, nd_v)] -=
                        zt2 * (j - 1) * zdisc * ztemp;
                    mpole[MPIDX(idim, 3, j, nd_v)] +=
                        zt2 * (j - 1) * zdisc;
                    mpole[MPIDX(idim, 2, j, nd_v)] +=
                        zt3 * zdisc;
                }
                zdis = zdis * ztemp * rinv;
                zdisc = zdisc / ztempc * rinv;
            }
        }
    }
}


/* ================================================================ */
/* bh2dformtac - form local expansion from charges                   */
/* ================================================================ */
void FNAME(bh2dformtac)(const fint *nd, const double *rscale,
                        const double *sources, const fint *ns,
                        const fcomplex *c1, const double *center,
                        const fint *nterms, fcomplex *mpole)
{
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nterms_v = *nterms;
    double rscale_v = *rscale;

    fint i, j, idim;
    fcomplex zc, zs, zdis, ztemp, zdisc, ztempc;

    zc = center[0] + center[1] * I;

    for (i = 1; i <= ns_v; i++) {
        zs = sources[FA2(1, i, 2)] + sources[FA2(2, i, 2)] * I;
        zdis = 1.0;
        ztemp = 1.0 / (zs - zc);
        zdisc = conj(zdis);
        ztempc = conj(ztemp);
        for (j = 0; j <= nterms_v; j++) {
            for (idim = 1; idim <= nd_v; idim++) {
                if (j == 0) {
                    mpole[MPIDX(idim, 4, j, nd_v)] +=
                        creal(2.0 * c1[FA3(idim, 1, i, nd_v, 2)])
                        * log(cabs(1.0 / ztemp));
                    mpole[MPIDX(idim, 5, j, nd_v)] +=
                        cimag(2.0 * c1[FA3(idim, 1, i, nd_v, 2)])
                        * log(cabs(1.0 / ztemp));
                } else {
                    mpole[MPIDX(idim, 4, j, nd_v)] -=
                        creal(2.0 * c1[FA3(idim, 1, i, nd_v, 2)])
                        * zdis / j;
                    mpole[MPIDX(idim, 5, j, nd_v)] -=
                        cimag(2.0 * c1[FA3(idim, 1, i, nd_v, 2)])
                        * zdis / j;
                }
                mpole[MPIDX(idim, 2, j, nd_v)] +=
                    c1[FA3(idim, 2, i, nd_v, 2)] * zdisc * ztempc / ztemp;
                mpole[MPIDX(idim, 3, j, nd_v)] -=
                    c1[FA3(idim, 2, i, nd_v, 2)] * zdisc * ztempc;
            }
            zdis = zdis * ztemp * rscale_v;
            zdisc = zdisc * ztempc * rscale_v;
        }
    }
}


/* ================================================================ */
/* bh2dformtad - form local expansion from dipoles                   */
/* ================================================================ */
void FNAME(bh2dformtad)(const fint *nd, const double *rscale,
                        const double *sources, const fint *ns,
                        const fcomplex *dip, const double *center,
                        const fint *nterms, fcomplex *mpole)
{
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nterms_v = *nterms;
    double rscale_v = *rscale;

    fint i, j, idim;
    fcomplex zc, zs, zdis, ztemp, zdisc, ztempc;

    zc = center[0] + center[1] * I;

    for (i = 1; i <= ns_v; i++) {
        zs = sources[FA2(1, i, 2)] + sources[FA2(2, i, 2)] * I;
        zdis = 1.0;
        ztemp = 1.0 / (zs - zc);
        zdisc = conj(zdis);
        ztempc = conj(ztemp);
        for (j = 0; j <= nterms_v; j++) {
            for (idim = 1; idim <= nd_v; idim++) {
                /* mp(1,j) -= dip(1,i)*zdis*ztemp */
                mpole[MPIDX(idim, 1, j, nd_v)] -=
                    dip[FA3(idim, 1, i, nd_v, 3)] * zdis * ztemp;
                /* mp(2,j) -= dip(3,i)*zdisc*ztempc */
                mpole[MPIDX(idim, 2, j, nd_v)] -=
                    dip[FA3(idim, 3, i, nd_v, 3)] * zdisc * ztempc;
                /* mp(2,j) -= dip(2,i)*(j+1)*zdisc*ztempc*ztempc/ztemp */
                mpole[MPIDX(idim, 2, j, nd_v)] -=
                    dip[FA3(idim, 2, i, nd_v, 3)] * (j + 1)
                    * zdisc * ztempc * ztempc / ztemp;
                /* mp(3,j) += dip(2,i)*(j+1)*zdisc*ztempc*ztempc */
                mpole[MPIDX(idim, 3, j, nd_v)] +=
                    dip[FA3(idim, 2, i, nd_v, 3)] * (j + 1)
                    * zdisc * ztempc * ztempc;
            }
            zdis = zdis * ztemp * rscale_v;
            zdisc = zdisc * ztempc * rscale_v;
        }
    }
}


/* ================================================================ */
/* bh2dformtacd - form local expansion from charges + dipoles        */
/* ================================================================ */
void FNAME(bh2dformtacd)(const fint *nd, const double *rscale,
                         const double *sources, const fint *ns,
                         const fcomplex *c1, const fcomplex *dip,
                         const double *center, const fint *nterms,
                         fcomplex *mpole)
{
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nterms_v = *nterms;
    double rscale_v = *rscale;

    fint i, j, idim;
    fcomplex zc, zs, zdis, ztemp, zdisc, ztempc;

    zc = center[0] + center[1] * I;

    for (i = 1; i <= ns_v; i++) {
        zs = sources[FA2(1, i, 2)] + sources[FA2(2, i, 2)] * I;
        zdis = 1.0;
        ztemp = 1.0 / (zs - zc);
        zdisc = conj(zdis);
        ztempc = conj(ztemp);
        for (j = 0; j <= nterms_v; j++) {
            for (idim = 1; idim <= nd_v; idim++) {
                if (j == 0) {
                    mpole[MPIDX(idim, 4, j, nd_v)] +=
                        creal(2.0 * c1[FA3(idim, 1, i, nd_v, 2)])
                        * log(cabs(1.0 / ztemp));
                    mpole[MPIDX(idim, 5, j, nd_v)] +=
                        cimag(2.0 * c1[FA3(idim, 1, i, nd_v, 2)])
                        * log(cabs(1.0 / ztemp));
                } else {
                    mpole[MPIDX(idim, 4, j, nd_v)] -=
                        creal(2.0 * c1[FA3(idim, 1, i, nd_v, 2)])
                        * zdis / j;
                    mpole[MPIDX(idim, 5, j, nd_v)] -=
                        cimag(2.0 * c1[FA3(idim, 1, i, nd_v, 2)])
                        * zdis / j;
                }
                mpole[MPIDX(idim, 2, j, nd_v)] +=
                    c1[FA3(idim, 2, i, nd_v, 2)] * zdisc * ztempc / ztemp;
                mpole[MPIDX(idim, 3, j, nd_v)] -=
                    c1[FA3(idim, 2, i, nd_v, 2)] * zdisc * ztempc;

                mpole[MPIDX(idim, 1, j, nd_v)] -=
                    dip[FA3(idim, 1, i, nd_v, 3)] * zdis * ztemp;
                mpole[MPIDX(idim, 2, j, nd_v)] -=
                    dip[FA3(idim, 3, i, nd_v, 3)] * zdisc * ztempc;
                mpole[MPIDX(idim, 2, j, nd_v)] -=
                    dip[FA3(idim, 2, i, nd_v, 3)] * (j + 1)
                    * zdisc * ztempc * ztempc / ztemp;
                mpole[MPIDX(idim, 3, j, nd_v)] +=
                    dip[FA3(idim, 2, i, nd_v, 3)] * (j + 1)
                    * zdisc * ztempc * ztempc;
            }
            zdis = zdis * ztemp * rscale_v;
            zdisc = zdisc * ztempc * rscale_v;
        }
    }
}


/* ================================================================ */
/* bh2dlocloc - local-to-local translation                           */
/* ================================================================ */
void FNAME(bh2dlocloc)(const fint *nd,
                       const double *rscale1, const double *c1,
                       const fcomplex *hexp, const fint *nterms1,
                       const double *rscale2, const double *c2,
                       fcomplex *jexp, const fint *nterms2,
                       const double *carray, const fint *ldc)
{
    fint nd_v = *nd;
    fint nterms1_v = *nterms1;
    fint nterms2_v = *nterms2;
    fint ldc_v = *ldc;
    double rscale1_v = *rscale1;
    double rscale2_v = *rscale2;
    double rinv1 = 1.0 / rscale1_v;

    fint i, j, idim, nmax;
    fcomplex zc1, zc2, zdis, zdisinv;
    fcomplex *zpow1, *zpow2;
    fcomplex *h1a, *h2a, *h3a, *h4a, *h5a;
    fcomplex *j1, *j2, *j3, *j4, *j5;

    zc1 = c1[0] + c1[1] * I;
    zc2 = c2[0] + c2[1] * I;
    zdis = zc2 - zc1;
    zdisinv = 1.0 / zdis;

    nmax = nterms1_v;
    if (nterms2_v > nmax) nmax = nterms2_v;

    zpow1 = (fcomplex *)malloc((nmax + 1) * sizeof(fcomplex));
    zpow2 = (fcomplex *)malloc((nmax + 1) * sizeof(fcomplex));
    h1a = (fcomplex *)malloc(nd_v * (nterms1_v + 1) * sizeof(fcomplex));
    h2a = (fcomplex *)malloc(nd_v * (nterms1_v + 1) * sizeof(fcomplex));
    h3a = (fcomplex *)malloc(nd_v * (nterms1_v + 1) * sizeof(fcomplex));
    h4a = (fcomplex *)malloc(nd_v * (nterms1_v + 1) * sizeof(fcomplex));
    h5a = (fcomplex *)malloc(nd_v * (nterms1_v + 1) * sizeof(fcomplex));
    j1 = (fcomplex *)malloc(nd_v * (nterms2_v + 1) * sizeof(fcomplex));
    j2 = (fcomplex *)malloc(nd_v * (nterms2_v + 1) * sizeof(fcomplex));
    j3 = (fcomplex *)malloc(nd_v * (nterms2_v + 1) * sizeof(fcomplex));
    j4 = (fcomplex *)malloc(nd_v * (nterms2_v + 1) * sizeof(fcomplex));
    j5 = (fcomplex *)malloc(nd_v * (nterms2_v + 1) * sizeof(fcomplex));

    zpow1[0] = 1.0;
    zpow2[0] = 1.0;

    for (i = 0; i <= nterms2_v; i++) {
        for (idim = 1; idim <= nd_v; idim++) {
            j1[HIDX(idim, i, nd_v)] = 0;
            j2[HIDX(idim, i, nd_v)] = 0;
            j3[HIDX(idim, i, nd_v)] = 0;
            j4[HIDX(idim, i, nd_v)] = 0;
            j5[HIDX(idim, i, nd_v)] = 0;
        }
    }

    for (i = 1; i <= nmax; i++) {
        zpow1[i] = zpow1[i - 1] * zdis * rinv1;
        zpow2[i] = zpow2[i - 1] * rscale2_v * zdisinv;
    }

    for (i = 0; i <= nterms1_v; i++) {
        for (idim = 1; idim <= nd_v; idim++) {
            h1a[HIDX(idim, i, nd_v)] =
                hexp[MPIDX(idim, 1, i, nd_v)] * zpow1[i];
            h2a[HIDX(idim, i, nd_v)] =
                (hexp[MPIDX(idim, 2, i, nd_v)]
                 + hexp[MPIDX(idim, 3, i, nd_v)] * zdis)
                * conj(zpow1[i]);
            h3a[HIDX(idim, i, nd_v)] =
                hexp[MPIDX(idim, 3, i, nd_v)] * conj(zpow1[i]);
            h4a[HIDX(idim, i, nd_v)] =
                hexp[MPIDX(idim, 4, i, nd_v)] * zpow1[i];
            h5a[HIDX(idim, i, nd_v)] =
                hexp[MPIDX(idim, 5, i, nd_v)] * zpow1[i];
        }
    }

    for (i = 0; i <= nterms2_v; i++) {
        for (j = i; j <= nterms1_v; j++) {
            for (idim = 1; idim <= nd_v; idim++) {
                j1[HIDX(idim, i, nd_v)] +=
                    h1a[HIDX(idim, j, nd_v)] * carray[CIDX(j, i, ldc_v)];
                j2[HIDX(idim, i, nd_v)] +=
                    h2a[HIDX(idim, j, nd_v)] * carray[CIDX(j, i, ldc_v)];
                j3[HIDX(idim, i, nd_v)] +=
                    h3a[HIDX(idim, j, nd_v)] * carray[CIDX(j, i, ldc_v)];
                j4[HIDX(idim, i, nd_v)] +=
                    h4a[HIDX(idim, j, nd_v)] * carray[CIDX(j, i, ldc_v)];
                j5[HIDX(idim, i, nd_v)] +=
                    h5a[HIDX(idim, j, nd_v)] * carray[CIDX(j, i, ldc_v)];
            }
        }
        for (idim = 1; idim <= nd_v; idim++) {
            j1[HIDX(idim, i, nd_v)] = j1[HIDX(idim, i, nd_v)] * zpow2[i];
            j2[HIDX(idim, i, nd_v)] = j2[HIDX(idim, i, nd_v)] * conj(zpow2[i]);
            j3[HIDX(idim, i, nd_v)] = j3[HIDX(idim, i, nd_v)] * conj(zpow2[i]);
            j4[HIDX(idim, i, nd_v)] = j4[HIDX(idim, i, nd_v)] * zpow2[i];
            j5[HIDX(idim, i, nd_v)] = j5[HIDX(idim, i, nd_v)] * zpow2[i];
            jexp[MPIDX(idim, 1, i, nd_v)] += j1[HIDX(idim, i, nd_v)];
            jexp[MPIDX(idim, 2, i, nd_v)] += j2[HIDX(idim, i, nd_v)];
            jexp[MPIDX(idim, 3, i, nd_v)] += j3[HIDX(idim, i, nd_v)];
            jexp[MPIDX(idim, 4, i, nd_v)] += j4[HIDX(idim, i, nd_v)];
            jexp[MPIDX(idim, 5, i, nd_v)] += j5[HIDX(idim, i, nd_v)];
        }
    }

    free(zpow1); free(zpow2);
    free(h1a); free(h2a); free(h3a); free(h4a); free(h5a);
    free(j1); free(j2); free(j3); free(j4); free(j5);
}


/* ================================================================ */
/* bh2dmpmp - mpole-to-mpole translation                             */
/* ================================================================ */
void FNAME(bh2dmpmp)(const fint *nd,
                     const double *rscale1, const double *c1,
                     const fcomplex *hexp, const fint *nterms1,
                     const double *rscale2, const double *c2,
                     fcomplex *jexp, const fint *nterms2,
                     const double *carray, const fint *ldc)
{
    fint nd_v = *nd;
    fint nterms1_v = *nterms1;
    fint nterms2_v = *nterms2;
    fint ldc_v = *ldc;
    double rscale1_v = *rscale1;
    double rscale2_v = *rscale2;
    double rinv2 = 1.0 / rscale2_v;

    fint i, j, idim, nmax, jmax;
    fcomplex zc1, zc2, zdis, zdisinv;
    fcomplex *zpow1, *zpow2;
    fcomplex *h1a, *h2a, *h3a, *h4a, *h5a;
    fcomplex *j1, *j2, *j3, *j4, *j5;

    nmax = nterms1_v;
    if (nterms2_v > nmax) nmax = nterms2_v;

    zc1 = c1[0] + c1[1] * I;
    zc2 = c2[0] + c2[1] * I;
    zdis = zc1 - zc2;
    zdisinv = 1.0 / zdis;

    zpow1 = (fcomplex *)malloc((nmax + 1) * sizeof(fcomplex));
    zpow2 = (fcomplex *)malloc((nmax + 1) * sizeof(fcomplex));
    h1a = (fcomplex *)malloc(nd_v * (nterms1_v + 1) * sizeof(fcomplex));
    h2a = (fcomplex *)malloc(nd_v * (nterms1_v + 1) * sizeof(fcomplex));
    h3a = (fcomplex *)malloc(nd_v * (nterms1_v + 1) * sizeof(fcomplex));
    h4a = (fcomplex *)malloc(nd_v * (nterms1_v + 1) * sizeof(fcomplex));
    h5a = (fcomplex *)malloc(nd_v * (nterms1_v + 1) * sizeof(fcomplex));
    j1 = (fcomplex *)malloc(nd_v * (nterms2_v + 1) * sizeof(fcomplex));
    j2 = (fcomplex *)malloc(nd_v * (nterms2_v + 1) * sizeof(fcomplex));
    j3 = (fcomplex *)malloc(nd_v * (nterms2_v + 1) * sizeof(fcomplex));
    j4 = (fcomplex *)malloc(nd_v * (nterms2_v + 1) * sizeof(fcomplex));
    j5 = (fcomplex *)malloc(nd_v * (nterms2_v + 1) * sizeof(fcomplex));

    zpow1[0] = 1.0;
    zpow2[0] = 1.0;
    for (i = 1; i <= nmax; i++) {
        zpow1[i] = zpow1[i - 1] * zdisinv * rscale1_v;
        zpow2[i] = zpow2[i - 1] * zdis * rinv2;
    }

    for (i = 0; i <= nterms2_v; i++) {
        for (idim = 1; idim <= nd_v; idim++) {
            j1[HIDX(idim, i, nd_v)] = 0.0;
            j2[HIDX(idim, i, nd_v)] = 0.0;
            j3[HIDX(idim, i, nd_v)] = 0.0;
            j4[HIDX(idim, i, nd_v)] = 0.0;
            j5[HIDX(idim, i, nd_v)] = 0.0;
        }
    }

    /* Handle log term in the expansion */
    for (idim = 1; idim <= nd_v; idim++) {
        jexp[MPIDX(idim, 4, 0, nd_v)] += hexp[MPIDX(idim, 4, 0, nd_v)];
        jexp[MPIDX(idim, 5, 0, nd_v)] += hexp[MPIDX(idim, 5, 0, nd_v)];
    }
    for (i = 1; i <= nterms2_v; i++) {
        for (idim = 1; idim <= nd_v; idim++) {
            j4[HIDX(idim, i, nd_v)] -=
                hexp[MPIDX(idim, 4, 0, nd_v)] / i;
            j5[HIDX(idim, i, nd_v)] -=
                hexp[MPIDX(idim, 5, 0, nd_v)] / i;
        }
    }

    jmax = nterms1_v;
    if (nterms2_v < jmax) jmax = nterms2_v;
    for (i = 1; i <= jmax; i++) {
        for (idim = 1; idim <= nd_v; idim++) {
            h1a[HIDX(idim, i, nd_v)] =
                hexp[MPIDX(idim, 1, i, nd_v)] * zpow1[i];
            h2a[HIDX(idim, i, nd_v)] =
                (hexp[MPIDX(idim, 2, i, nd_v)]
                 - hexp[MPIDX(idim, 3, i, nd_v)] * zdis)
                * conj(zpow1[i]);
            h3a[HIDX(idim, i, nd_v)] =
                hexp[MPIDX(idim, 3, i, nd_v)] * conj(zpow1[i]);
            h4a[HIDX(idim, i, nd_v)] =
                hexp[MPIDX(idim, 4, i, nd_v)] * zpow1[i];
            h5a[HIDX(idim, i, nd_v)] =
                hexp[MPIDX(idim, 5, i, nd_v)] * zpow1[i];
        }
    }

    for (i = 1; i <= nterms2_v; i++) {
        fint jhi = i;
        if (nterms1_v < jhi) jhi = nterms1_v;
        for (j = 1; j <= jhi; j++) {
            for (idim = 1; idim <= nd_v; idim++) {
                j1[HIDX(idim, i, nd_v)] +=
                    h1a[HIDX(idim, j, nd_v)]
                    * carray[CIDX(i - 1, j - 1, ldc_v)];
                j2[HIDX(idim, i, nd_v)] +=
                    h2a[HIDX(idim, j, nd_v)]
                    * carray[CIDX(i - 1, j - 1, ldc_v)];
                j3[HIDX(idim, i, nd_v)] +=
                    h3a[HIDX(idim, j, nd_v)]
                    * carray[CIDX(i - 1, j - 1, ldc_v)];
                j4[HIDX(idim, i, nd_v)] +=
                    h4a[HIDX(idim, j, nd_v)]
                    * carray[CIDX(i - 1, j - 1, ldc_v)];
                j5[HIDX(idim, i, nd_v)] +=
                    h5a[HIDX(idim, j, nd_v)]
                    * carray[CIDX(i - 1, j - 1, ldc_v)];
            }
        }
        for (idim = 1; idim <= nd_v; idim++) {
            j1[HIDX(idim, i, nd_v)] = j1[HIDX(idim, i, nd_v)] * zpow2[i];
            j2[HIDX(idim, i, nd_v)] = j2[HIDX(idim, i, nd_v)] * conj(zpow2[i]);
            j3[HIDX(idim, i, nd_v)] = j3[HIDX(idim, i, nd_v)] * conj(zpow2[i]);
            j4[HIDX(idim, i, nd_v)] = j4[HIDX(idim, i, nd_v)] * zpow2[i];
            j5[HIDX(idim, i, nd_v)] = j5[HIDX(idim, i, nd_v)] * zpow2[i];
            jexp[MPIDX(idim, 1, i, nd_v)] += j1[HIDX(idim, i, nd_v)];
            jexp[MPIDX(idim, 2, i, nd_v)] += j2[HIDX(idim, i, nd_v)];
            jexp[MPIDX(idim, 3, i, nd_v)] += j3[HIDX(idim, i, nd_v)];
            jexp[MPIDX(idim, 4, i, nd_v)] += j4[HIDX(idim, i, nd_v)];
            jexp[MPIDX(idim, 5, i, nd_v)] += j5[HIDX(idim, i, nd_v)];
        }
    }

    free(zpow1); free(zpow2);
    free(h1a); free(h2a); free(h3a); free(h4a); free(h5a);
    free(j1); free(j2); free(j3); free(j4); free(j5);
}


/* ================================================================ */
/* bh2dmploc - mpole-to-local translation                            */
/* ================================================================ */
void FNAME(bh2dmploc)(const fint *nd,
                      const double *rscale1, const double *c1,
                      const fcomplex *hexp, const fint *nterms1,
                      const double *rscale2, const double *c2,
                      fcomplex *jexp, const fint *nterms2,
                      const double *carray, const fint *ldc)
{
    fint nd_v = *nd;
    fint nterms1_v = *nterms1;
    fint nterms2_v = *nterms2;
    fint ldc_v = *ldc;
    double rscale1_v = *rscale1;
    double rscale2_v = *rscale2;

    fint i, j, idim, nmax;
    fcomplex zc1, zc2, zdis, zdis1;
    fcomplex *zpow1, *zpow2;
    fcomplex *h1a, *h2a, *h3a, *h4a, *h5a;
    fcomplex *j1, *j2, *j3, *j4, *j5;

    nmax = nterms1_v + nterms2_v;

    zc1 = c1[0] + c1[1] * I;
    zc2 = c2[0] + c2[1] * I;
    zdis = 1.0 / (zc1 - zc2);
    zdis1 = (zc1 - zc2);

    zpow1 = (fcomplex *)malloc((nmax + 1) * sizeof(fcomplex));
    zpow2 = (fcomplex *)malloc((nmax + 1) * sizeof(fcomplex));
    h1a = (fcomplex *)malloc(nd_v * (nterms1_v + 1) * sizeof(fcomplex));
    h2a = (fcomplex *)malloc(nd_v * (nterms1_v + 1) * sizeof(fcomplex));
    h3a = (fcomplex *)malloc(nd_v * (nterms1_v + 1) * sizeof(fcomplex));
    h4a = (fcomplex *)malloc(nd_v * (nterms1_v + 1) * sizeof(fcomplex));
    h5a = (fcomplex *)malloc(nd_v * (nterms1_v + 1) * sizeof(fcomplex));
    j1 = (fcomplex *)malloc(nd_v * (nterms2_v + 1) * sizeof(fcomplex));
    j2 = (fcomplex *)malloc(nd_v * (nterms2_v + 1) * sizeof(fcomplex));
    j3 = (fcomplex *)malloc(nd_v * (nterms2_v + 1) * sizeof(fcomplex));
    j4 = (fcomplex *)malloc(nd_v * (nterms2_v + 1) * sizeof(fcomplex));
    j5 = (fcomplex *)malloc(nd_v * (nterms2_v + 1) * sizeof(fcomplex));

    zpow1[0] = 1.0;
    zpow2[0] = 1.0;
    for (i = 1; i <= nmax; i++) {
        zpow1[i] = -zpow1[i - 1] * zdis * rscale1_v;
        zpow2[i] = zpow2[i - 1] * zdis * rscale2_v;
    }

    for (i = 0; i <= nterms2_v; i++) {
        for (idim = 1; idim <= nd_v; idim++) {
            j1[HIDX(idim, i, nd_v)] = 0.0;
            j2[HIDX(idim, i, nd_v)] = 0.0;
            j3[HIDX(idim, i, nd_v)] = 0.0;
            j4[HIDX(idim, i, nd_v)] = 0.0;
            j5[HIDX(idim, i, nd_v)] = 0.0;
        }
    }

    /* Handle the log term in the expansion */
    for (idim = 1; idim <= nd_v; idim++) {
        j4[HIDX(idim, 0, nd_v)] =
            hexp[MPIDX(idim, 4, 0, nd_v)] * log(cabs(zdis1));
        j5[HIDX(idim, 0, nd_v)] =
            hexp[MPIDX(idim, 5, 0, nd_v)] * log(cabs(zdis1));
    }
    for (i = 1; i <= nterms2_v; i++) {
        for (idim = 1; idim <= nd_v; idim++) {
            j4[HIDX(idim, i, nd_v)] -=
                hexp[MPIDX(idim, 4, 0, nd_v)] / i;
            j5[HIDX(idim, i, nd_v)] -=
                hexp[MPIDX(idim, 5, 0, nd_v)] / i;
        }
    }

    for (i = 1; i <= nterms1_v; i++) {
        for (idim = 1; idim <= nd_v; idim++) {
            h1a[HIDX(idim, i, nd_v)] =
                hexp[MPIDX(idim, 1, i, nd_v)] * zpow1[i];
            h2a[HIDX(idim, i, nd_v)] =
                (hexp[MPIDX(idim, 2, i, nd_v)]
                 - hexp[MPIDX(idim, 3, i, nd_v)] * zdis1)
                * conj(zpow1[i]);
            h3a[HIDX(idim, i, nd_v)] =
                hexp[MPIDX(idim, 3, i, nd_v)] * conj(zpow1[i]);
            h4a[HIDX(idim, i, nd_v)] =
                hexp[MPIDX(idim, 4, i, nd_v)] * zpow1[i];
            h5a[HIDX(idim, i, nd_v)] =
                hexp[MPIDX(idim, 5, i, nd_v)] * zpow1[i];
        }
    }

    for (i = 0; i <= nterms2_v; i++) {
        for (j = 1; j <= nterms1_v; j++) {
            for (idim = 1; idim <= nd_v; idim++) {
                j1[HIDX(idim, i, nd_v)] +=
                    h1a[HIDX(idim, j, nd_v)]
                    * carray[CIDX(i + j - 1, i, ldc_v)];
                j2[HIDX(idim, i, nd_v)] +=
                    h2a[HIDX(idim, j, nd_v)]
                    * carray[CIDX(i + j - 1, i, ldc_v)];
                j3[HIDX(idim, i, nd_v)] +=
                    h3a[HIDX(idim, j, nd_v)]
                    * carray[CIDX(i + j - 1, i, ldc_v)];
                j4[HIDX(idim, i, nd_v)] +=
                    h4a[HIDX(idim, j, nd_v)]
                    * carray[CIDX(i + j - 1, i, ldc_v)];
                j5[HIDX(idim, i, nd_v)] +=
                    h5a[HIDX(idim, j, nd_v)]
                    * carray[CIDX(i + j - 1, i, ldc_v)];
            }
        }
        for (idim = 1; idim <= nd_v; idim++) {
            j1[HIDX(idim, i, nd_v)] = j1[HIDX(idim, i, nd_v)] * zpow2[i];
            j2[HIDX(idim, i, nd_v)] = j2[HIDX(idim, i, nd_v)] * conj(zpow2[i]);
            j3[HIDX(idim, i, nd_v)] = j3[HIDX(idim, i, nd_v)] * conj(zpow2[i]);
            j4[HIDX(idim, i, nd_v)] = j4[HIDX(idim, i, nd_v)] * zpow2[i];
            j5[HIDX(idim, i, nd_v)] = j5[HIDX(idim, i, nd_v)] * zpow2[i];
            jexp[MPIDX(idim, 1, i, nd_v)] += j1[HIDX(idim, i, nd_v)];
            jexp[MPIDX(idim, 2, i, nd_v)] += j2[HIDX(idim, i, nd_v)];
            jexp[MPIDX(idim, 3, i, nd_v)] += j3[HIDX(idim, i, nd_v)];
            jexp[MPIDX(idim, 4, i, nd_v)] += j4[HIDX(idim, i, nd_v)];
            jexp[MPIDX(idim, 5, i, nd_v)] += j5[HIDX(idim, i, nd_v)];
        }
    }

    free(zpow1); free(zpow2);
    free(h1a); free(h2a); free(h3a); free(h4a); free(h5a);
    free(j1); free(j2); free(j3); free(j4); free(j5);
}


/* ================================================================ */
/* bh2dmpzero - zero out a coefficient buffer                        */
/* ================================================================ */
void FNAME(bh2dmpzero)(const fint *nd, fcomplex *mpole, const fint *nterms)
{
    fint nd_v = *nd;
    fint nterms_v = *nterms;
    fint i, j, idim;

    for (i = 0; i <= nterms_v; i++) {
        for (j = 1; j <= 5; j++) {
            for (idim = 1; idim <= nd_v; idim++) {
                mpole[MPIDX(idim, j, i, nd_v)] = 0;
            }
        }
    }
}
