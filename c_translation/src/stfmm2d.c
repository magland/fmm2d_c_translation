/*
 * stfmm2d.c - C translation of src/stokes/stfmm2d.f
 *
 * 2D Stokes FMM driver. Thin wrapper around bhfmm2d:
 *   1. Allocate complex *16 (nd, 2, nsource) charge and (nd, 3, nsource)
 *      dip work buffers, plus complex (nd, *) potl/pottargl and
 *      (nd, 3, *) gradl/gradtargl.
 *   2. Encode the Stokeslet stoklet(nd,2,*) into the biharmonic
 *      complex charge(nd,2,*) array, and the stresslet
 *      strslet(nd,2,*)/strsvec(nd,2,*) into the biharmonic
 *      complex dip(nd,3,*) array.
 *   3. Call bhfmm2d.
 *   4. Decode the biharmonic potl/gradl outputs back into the Stokes
 *      pot (velocity, real (nd,2,*)), pre (pressure, real (nd,*))
 *      and grad (velocity gradient, real (nd,2,2,*)).
 *
 * Translated 1:1 from the Fortran reference. The Fortran source's
 * grad output (ifppreg/ifppregtarg .ge. 3) is implemented but the
 * Fortran source comments note "GRADIENT NOT IMPLEMENTED" as a
 * caveat — we still translate the code path for fidelity.
 *
 * Cross-file call to bhfmm2d is via the bare Fortran symbol name so
 * the diff test isolates this file from bhfmm2d's own implementation.
 */

#include "stfmm2d.h"

extern void bhfmm2d_(const fint *nd, const double *eps,
                     const fint *ns, const double *sources,
                     const fint *ifcharge, const fcomplex *charge,
                     const fint *ifdipole, const fcomplex *dip,
                     fint *iper, const fint *ifpgh,
                     fcomplex *pot, fcomplex *grad, fcomplex *hess,
                     const fint *nt, const double *targ,
                     const fint *ifpghtarg,
                     fcomplex *pottarg, fcomplex *gradtarg,
                     fcomplex *hesstarg, fint *ier);


void FNAME(stfmm2d)(const fint *nd, const double *eps,
                    const fint *nsource, const double *source,
                    const fint *ifstoklet, const double *stoklet,
                    const fint *ifstrslet, const double *strslet,
                    const double *strsvec,
                    const fint *ifppreg, double *pot, double *pre,
                    double *grad,
                    const fint *ntarg, const double *targ,
                    const fint *ifppregtarg, double *pottarg,
                    double *pretarg, double *gradtarg, fint *ier)
{
    fint nd_v = *nd;
    fint nsource_v = *nsource;
    fint ntarg_v = *ntarg;
    fint ifstoklet_v = *ifstoklet;
    fint ifstrslet_v = *ifstrslet;
    fint ifppreg_v = *ifppreg;
    fint ifppregtarg_v = *ifppregtarg;

    fint ifchargel = 0, ifdipolel = 0;
    fint ifpghl = 0, ifpghtargl = 0;
    fint iper;

    /* hesstmp(10): single complex *16 placeholder array shared between
     * the hess and hesstarg arguments to bhfmm2d (which is unused). */
    fcomplex hesstmp[10];

    fcomplex *charge = NULL, *dip = NULL;
    fcomplex *potl = NULL, *gradl = NULL;
    fcomplex *pottargl = NULL, *gradtargl = NULL;
    fcomplex *zsum = NULL;
    fcomplex zd1, zd2;
    fcomplex tmp;
    fint i, j;

    if (ifstoklet_v == 1) ifchargel = 1;
    if (ifstrslet_v == 1) ifdipolel = 1;

    /* allocate work arrays */
    charge = (fcomplex *)malloc(nd_v * 2 * nsource_v * sizeof(fcomplex));
    dip    = (fcomplex *)malloc(nd_v * 3 * nsource_v * sizeof(fcomplex));
    potl   = (fcomplex *)malloc(nd_v * nsource_v * sizeof(fcomplex));
    pottargl = (fcomplex *)malloc(nd_v * ntarg_v * sizeof(fcomplex));
    gradl  = (fcomplex *)malloc(nd_v * 3 * nsource_v * sizeof(fcomplex));
    gradtargl = (fcomplex *)malloc(nd_v * 3 * ntarg_v * sizeof(fcomplex));
    zsum = (fcomplex *)malloc(nd_v * sizeof(fcomplex));

    /* Initialize hesstmp (10 complex slots, used as a placeholder). */
    for (i = 0; i < 10; i++) hesstmp[i] = 0;

    /* zsum = 0 */
    for (j = 1; j <= nd_v; j++) {
        zsum[j - 1] = 0;
    }

    /* set-up appropriate vector charge and dipole arrays */
    for (i = 1; i <= nsource_v; i++) {
        for (j = 1; j <= nd_v; j++) {
            charge[FA3(j, 1, i, nd_v, 2)] = 0;
            charge[FA3(j, 2, i, nd_v, 2)] = 0;
            dip[FA3(j, 1, i, nd_v, 3)] = 0;
            dip[FA3(j, 2, i, nd_v, 3)] = 0;
            dip[FA3(j, 3, i, nd_v, 3)] = 0;
        }

        if (ifstoklet_v == 1) {
            for (j = 1; j <= nd_v; j++) {
                /* charge(j,1,i) = (-i*stoklet(j,1,i) + stoklet(j,2,i))/4 */
                charge[FA3(j, 1, i, nd_v, 2)] =
                    (-I * stoklet[FA3(j, 1, i, nd_v, 2)]
                     + stoklet[FA3(j, 2, i, nd_v, 2)]) / 4.0;
                charge[FA3(j, 2, i, nd_v, 2)] =
                    conj(charge[FA3(j, 1, i, nd_v, 2)]);
                zsum[j - 1] -= charge[FA3(j, 1, i, nd_v, 2)];
            }
        }

        if (ifstrslet_v == 1) {
            for (j = 1; j <= nd_v; j++) {
                /* zd1 = -(-strslet(j,2,i) + i*strslet(j,1,i))/2 */
                zd1 = -(-strslet[FA3(j, 2, i, nd_v, 2)]
                        + I * strslet[FA3(j, 1, i, nd_v, 2)]) / 2.0;
                /* zd2 =  (-strsvec(j,2,i) + i*strsvec(j,1,i)) */
                zd2 = (-strsvec[FA3(j, 2, i, nd_v, 2)]
                       + I * strsvec[FA3(j, 1, i, nd_v, 2)]);
                /* dip(j,1,i) = -i*zd1*zd2 */
                dip[FA3(j, 1, i, nd_v, 3)] = -I * zd1 * zd2;
                /* dip(j,2,i) = -conj(dip(j,1,i)) */
                dip[FA3(j, 2, i, nd_v, 3)] =
                    -conj(dip[FA3(j, 1, i, nd_v, 3)]);
                /* dip(j,3,i) = i*(zd1*conj(zd2) + conj(zd1)*zd2) */
                dip[FA3(j, 3, i, nd_v, 3)] =
                    I * (zd1 * conj(zd2) + conj(zd1) * zd2);
            }
        }
    }

    /* Map (ifppreg, ifppregtarg) -> (ifpghl, ifpghtargl). bhfmm2d only
     * supports up to 2 (no hessian) so we cap at 2 for both. */
    if (ifppreg_v == 1) ifpghl = 1;
    if (ifppreg_v >= 2) ifpghl = 2;
    if (ifppregtarg_v == 1) ifpghtargl = 1;
    if (ifppregtarg_v >= 2) ifpghtargl = 2;

    /* call biharmonic FMM */
    iper = 0;
    *ier = 0;

    bhfmm2d_(&nd_v, eps, &nsource_v, source, &ifchargel, charge,
             &ifdipolel, dip, &iper, &ifpghl, potl, gradl, hesstmp,
             &ntarg_v, targ, &ifpghtargl, pottargl, gradtargl,
             hesstmp, ier);

    /* unpack: pot (velocity) at sources */
    if (ifppreg_v >= 1) {
        for (i = 1; i <= nsource_v; i++) {
            for (j = 1; j <= nd_v; j++) {
                /* tmp = potl(j,i) + zsum(j) + charge(j,1,i) */
                tmp = potl[FA2(j, i, nd_v)];
                tmp += zsum[j - 1];
                tmp += charge[FA3(j, 1, i, nd_v, 2)];
                /* pot(j,1,i) =  imag(tmp) */
                pot[FA3(j, 1, i, nd_v, 2)] = cimag(tmp);
                /* pot(j,2,i) = -real(tmp) */
                pot[FA3(j, 2, i, nd_v, 2)] = -creal(tmp);
            }
        }
    }

    /* unpack: pre (pressure) at sources */
    if (ifppreg_v >= 2) {
        for (i = 1; i <= nsource_v; i++) {
            for (j = 1; j <= nd_v; j++) {
                pre[FA2(j, i, nd_v)] =
                    -4.0 * cimag(gradl[FA3(j, 1, i, nd_v, 3)]);
            }
        }
    }

    /* unpack: grad (velocity gradient) at sources */
    if (ifppreg_v >= 3) {
        for (i = 1; i <= nsource_v; i++) {
            for (j = 1; j <= nd_v; j++) {
                /* grad(j,1,1,i) =  imag(gradl(j,3,i)) */
                grad[FA4(j, 1, 1, i, nd_v, 2, 2)] =
                    cimag(gradl[FA3(j, 3, i, nd_v, 3)]);
                /* grad(j,2,2,i) = -imag(gradl(j,3,i)) */
                grad[FA4(j, 2, 2, i, nd_v, 2, 2)] =
                    -cimag(gradl[FA3(j, 3, i, nd_v, 3)]);
                /* grad(j,2,1,i) =  real(2*gradl(j,1,i)-gradl(j,3,i)) */
                grad[FA4(j, 2, 1, i, nd_v, 2, 2)] =
                    creal(2.0 * gradl[FA3(j, 1, i, nd_v, 3)]
                          - gradl[FA3(j, 3, i, nd_v, 3)]);
                /* grad(j,1,2,i) = -real(2*gradl(j,1,i)+gradl(j,3,i)) */
                grad[FA4(j, 1, 2, i, nd_v, 2, 2)] =
                    -creal(2.0 * gradl[FA3(j, 1, i, nd_v, 3)]
                           + gradl[FA3(j, 3, i, nd_v, 3)]);
            }
        }
    }

    /* unpack: pot (velocity) at targets */
    if (ifppregtarg_v >= 1) {
        for (i = 1; i <= ntarg_v; i++) {
            for (j = 1; j <= nd_v; j++) {
                /* tmp = pottargl(j,i) + zsum(j)  (no charge term here) */
                tmp = pottargl[FA2(j, i, nd_v)];
                tmp += zsum[j - 1];
                pottarg[FA3(j, 1, i, nd_v, 2)] = cimag(tmp);
                pottarg[FA3(j, 2, i, nd_v, 2)] = -creal(tmp);
            }
        }
    }

    /* unpack: pre (pressure) at targets */
    if (ifppregtarg_v >= 2) {
        for (i = 1; i <= ntarg_v; i++) {
            for (j = 1; j <= nd_v; j++) {
                pretarg[FA2(j, i, nd_v)] =
                    -4.0 * cimag(gradtargl[FA3(j, 1, i, nd_v, 3)]);
            }
        }
    }

    /* unpack: grad (velocity gradient) at targets */
    if (ifppregtarg_v >= 3) {
        for (i = 1; i <= ntarg_v; i++) {
            for (j = 1; j <= nd_v; j++) {
                gradtarg[FA4(j, 1, 1, i, nd_v, 2, 2)] =
                    cimag(gradtargl[FA3(j, 3, i, nd_v, 3)]);
                gradtarg[FA4(j, 2, 2, i, nd_v, 2, 2)] =
                    -cimag(gradtargl[FA3(j, 3, i, nd_v, 3)]);
                gradtarg[FA4(j, 2, 1, i, nd_v, 2, 2)] =
                    creal(2.0 * gradtargl[FA3(j, 1, i, nd_v, 3)]
                          - gradtargl[FA3(j, 3, i, nd_v, 3)]);
                gradtarg[FA4(j, 1, 2, i, nd_v, 2, 2)] =
                    -creal(2.0 * gradtargl[FA3(j, 1, i, nd_v, 3)]
                           + gradtargl[FA3(j, 3, i, nd_v, 3)]);
            }
        }
    }

    free(charge);
    free(dip);
    free(potl);
    free(pottargl);
    free(gradl);
    free(gradtargl);
    free(zsum);
}
