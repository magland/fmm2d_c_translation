/*
 * lfmm2d_ndiv.c - C translation of src/laplace/lfmm2d_ndiv.f
 *
 * Top-level 2D Laplace FMM driver with COMPLEX-valued charges and
 * dipole strengths and externally supplied subdivision criterion. Thin
 * wrapper around cfmm2d_ndiv with nd doubled:
 *   1. Allocate complex *16 buffers charge1, dipstr1, pot1, grad1,
 *      hess1, pottarg1, gradtarg1, hesstarg1 of shape (2, nd, *).
 *   2. Convert complex inputs into the (re, im) pair packing:
 *        charge1(1,j,i) = dble(charge(j,i))      (zero imag)
 *        charge1(2,j,i) = imag(charge(j,i))      (zero imag)
 *        ztmp = -(dipvec(j,1,i) + I*dipvec(j,2,i))
 *        dipstr1(1,j,i) = dble(dipstr(j,i)) * ztmp
 *        dipstr1(2,j,i) = imag(dipstr(j,i)) * ztmp
 *   3. Call cfmm2d_ndiv with nd2 = 2*nd, viewing the (2, nd, *)
 *      buffers as (2*nd, *).
 *   4. Recombine the (2, nd, *) complex outputs back to the original
 *      complex shapes:
 *        pot(j,i)     = Re(pot1(1,j,i))    + i * Re(pot1(2,j,i))
 *        grad(j,1,i)  = Re(grad1(1,j,i))   + i * Re(grad1(2,j,i))
 *        grad(j,2,i)  = -Im(grad1(1,j,i))  - i * Im(grad1(2,j,i))
 *        hess(j,1,i)  = Re(hess1(1,j,i))   + i * Re(hess1(2,j,i))
 *        hess(j,2,i)  = -Im(hess1(1,j,i))  - i * Im(hess1(2,j,i))
 *        hess(j,3,i)  = -hess(j,1,i)
 *      (analogous for pottarg / gradtarg / hesstarg)
 *
 * Translated 1:1 from the Fortran reference (same control flow, same
 * allocations, same operation order). OpenMP, error/print logging, and
 * timing are stripped.
 *
 * Cross-file call to cfmm2d_ndiv is via the bare Fortran symbol name so
 * the diff test isolates this file from cfmm2d_ndiv's own
 * implementation.
 */

#include "lfmm2d_ndiv.h"

/* Cross-file call: bare Fortran symbol name. Signature mirrors
 * cfmm2d_ndiv.h (with FNAME stripped). */
extern void cfmm2d_ndiv_(const fint *nd, const double *eps,
                         const fint *ns, const double *sources,
                         const fint *ifcharge, const fcomplex *charge,
                         const fint *ifdipole, const fcomplex *dipstr,
                         fint *iper, const fint *ifpgh,
                         fcomplex *pot, fcomplex *grad, fcomplex *hess,
                         const fint *nt, const double *targ,
                         const fint *ifpghtarg,
                         fcomplex *pottarg, fcomplex *gradtarg,
                         fcomplex *hesstarg,
                         const fint *ndiv, const fint *idivflag,
                         const fint *ifnear, double *timeinfo, fint *ier);


/* ---------------------------------------------------------------- */
/* lfmm2d_ndiv - complex-valued top-level user-callable driver      */
/* ---------------------------------------------------------------- */
void FNAME(lfmm2d_ndiv)(const fint *nd, const double *eps,
                        const fint *ns, const double *sources,
                        const fint *ifcharge, const fcomplex *charge,
                        const fint *ifdipole, const fcomplex *dipstr,
                        const double *dipvec,
                        fint *iper, const fint *ifpgh,
                        fcomplex *pot, fcomplex *grad, fcomplex *hess,
                        const fint *nt, const double *targ,
                        const fint *ifpghtarg,
                        fcomplex *pottarg, fcomplex *gradtarg,
                        fcomplex *hesstarg,
                        const fint *ndiv, const fint *idivflag,
                        const fint *ifnear, double *timeinfo, fint *ier)
{
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    fint ifcharge_v = *ifcharge;
    fint ifdipole_v = *ifdipole;
    fint ifpgh_v = *ifpgh;
    fint ifpghtarg_v = *ifpghtarg;

    fcomplex *charge1 = NULL;
    fcomplex *dipstr1 = NULL;
    fcomplex *pot1 = NULL;
    fcomplex *grad1 = NULL;
    fcomplex *hess1 = NULL;
    fcomplex *pottarg1 = NULL;
    fcomplex *gradtarg1 = NULL;
    fcomplex *hesstarg1 = NULL;

    fint i, j;
    fcomplex ztmp;
    fint nd2;

    /* allocate variables for cfmm call. The Fortran allocates
     * complex *16 (2, nd, ns) etc. */

    if (ifcharge_v == 1) {
        charge1 = (fcomplex *)malloc(2 * nd_v * ns_v * sizeof(fcomplex));
    } else {
        charge1 = (fcomplex *)malloc(2 * nd_v * 1 * sizeof(fcomplex));
    }

    if (ifdipole_v == 1) {
        dipstr1 = (fcomplex *)malloc(2 * nd_v * ns_v * sizeof(fcomplex));
    } else {
        dipstr1 = (fcomplex *)malloc(2 * nd_v * 1 * sizeof(fcomplex));
    }

    if (ifpgh_v == 1) {
        pot1  = (fcomplex *)malloc(2 * nd_v * ns_v * sizeof(fcomplex));
        grad1 = (fcomplex *)malloc(2 * nd_v * 1    * sizeof(fcomplex));
        hess1 = (fcomplex *)malloc(2 * nd_v * 1    * sizeof(fcomplex));
    }
    if (ifpgh_v == 2) {
        pot1  = (fcomplex *)malloc(2 * nd_v * ns_v * sizeof(fcomplex));
        grad1 = (fcomplex *)malloc(2 * nd_v * ns_v * sizeof(fcomplex));
        hess1 = (fcomplex *)malloc(2 * nd_v * 1    * sizeof(fcomplex));
    }
    if (ifpgh_v == 3) {
        pot1  = (fcomplex *)malloc(2 * nd_v * ns_v * sizeof(fcomplex));
        grad1 = (fcomplex *)malloc(2 * nd_v * ns_v * sizeof(fcomplex));
        hess1 = (fcomplex *)malloc(2 * nd_v * ns_v * sizeof(fcomplex));
    }

    if (ifpghtarg_v == 1) {
        pottarg1  = (fcomplex *)malloc(2 * nd_v * nt_v * sizeof(fcomplex));
        gradtarg1 = (fcomplex *)malloc(2 * nd_v * 1    * sizeof(fcomplex));
        hesstarg1 = (fcomplex *)malloc(2 * nd_v * 1    * sizeof(fcomplex));
    }
    if (ifpghtarg_v == 2) {
        pottarg1  = (fcomplex *)malloc(2 * nd_v * nt_v * sizeof(fcomplex));
        gradtarg1 = (fcomplex *)malloc(2 * nd_v * nt_v * sizeof(fcomplex));
        hesstarg1 = (fcomplex *)malloc(2 * nd_v * 1    * sizeof(fcomplex));
    }
    if (ifpghtarg_v == 3) {
        pottarg1  = (fcomplex *)malloc(2 * nd_v * nt_v * sizeof(fcomplex));
        gradtarg1 = (fcomplex *)malloc(2 * nd_v * nt_v * sizeof(fcomplex));
        hesstarg1 = (fcomplex *)malloc(2 * nd_v * nt_v * sizeof(fcomplex));
    }

    /* real and complex parts must be split */

    if (ifcharge_v == 1) {
        for (i = 1; i <= ns_v; i++) {
            for (j = 1; j <= nd_v; j++) {
                charge1[FA3(1, j, i, 2, nd_v)] = creal(charge[FA2(j, i, nd_v)]);
                charge1[FA3(2, j, i, 2, nd_v)] = cimag(charge[FA2(j, i, nd_v)]);
            }
        }
    }

    if (ifdipole_v == 1) {
        for (i = 1; i <= ns_v; i++) {
            for (j = 1; j <= nd_v; j++) {
                ztmp = -(dipvec[FA3(j, 1, i, nd_v, 2)]
                         + dipvec[FA3(j, 2, i, nd_v, 2)] * I);
                dipstr1[FA3(1, j, i, 2, nd_v)] =
                    creal(dipstr[FA2(j, i, nd_v)]) * ztmp;
                dipstr1[FA3(2, j, i, 2, nd_v)] =
                    cimag(dipstr[FA2(j, i, nd_v)]) * ztmp;
            }
        }
    }

    /* cfmm does the work */

    nd2 = 2 * nd_v;
    cfmm2d_ndiv_(&nd2, eps, ns, sources, ifcharge, charge1,
                 ifdipole, dipstr1, iper, ifpgh, pot1, grad1, hess1,
                 nt, targ, ifpghtarg, pottarg1, gradtarg1,
                 hesstarg1, ndiv, idivflag, ifnear, timeinfo, ier);

    /* unpack the d/dz, d^2/dz^2 as grad/hess and combine real/imag */

    if (ifpgh_v == 1 || ifpgh_v == 2 || ifpgh_v == 3) {
        for (i = 1; i <= ns_v; i++) {
            for (j = 1; j <= nd_v; j++) {
                pot[FA2(j, i, nd_v)] =
                    creal(pot1[FA3(1, j, i, 2, nd_v)])
                    + I * creal(pot1[FA3(2, j, i, 2, nd_v)]);
            }
        }
    }
    if (ifpgh_v == 2 || ifpgh_v == 3) {
        for (i = 1; i <= ns_v; i++) {
            for (j = 1; j <= nd_v; j++) {
                grad[FA3(j, 1, i, nd_v, 2)] =
                    creal(grad1[FA3(1, j, i, 2, nd_v)])
                    + I * creal(grad1[FA3(2, j, i, 2, nd_v)]);
                grad[FA3(j, 2, i, nd_v, 2)] =
                    -cimag(grad1[FA3(1, j, i, 2, nd_v)])
                    - I * cimag(grad1[FA3(2, j, i, 2, nd_v)]);
            }
        }
    }
    if (ifpgh_v == 3) {
        for (i = 1; i <= ns_v; i++) {
            for (j = 1; j <= nd_v; j++) {
                hess[FA3(j, 1, i, nd_v, 3)] =
                    creal(hess1[FA3(1, j, i, 2, nd_v)])
                    + I * creal(hess1[FA3(2, j, i, 2, nd_v)]);
                hess[FA3(j, 2, i, nd_v, 3)] =
                    -cimag(hess1[FA3(1, j, i, 2, nd_v)])
                    - I * cimag(hess1[FA3(2, j, i, 2, nd_v)]);
                hess[FA3(j, 3, i, nd_v, 3)] =
                    -hess[FA3(j, 1, i, nd_v, 3)];
            }
        }
    }

    if (ifpghtarg_v == 1 || ifpghtarg_v == 2 || ifpghtarg_v == 3) {
        for (i = 1; i <= nt_v; i++) {
            for (j = 1; j <= nd_v; j++) {
                pottarg[FA2(j, i, nd_v)] =
                    creal(pottarg1[FA3(1, j, i, 2, nd_v)])
                    + I * creal(pottarg1[FA3(2, j, i, 2, nd_v)]);
            }
        }
    }
    if (ifpghtarg_v == 2 || ifpghtarg_v == 3) {
        for (i = 1; i <= nt_v; i++) {
            for (j = 1; j <= nd_v; j++) {
                gradtarg[FA3(j, 1, i, nd_v, 2)] =
                    creal(gradtarg1[FA3(1, j, i, 2, nd_v)])
                    + I * creal(gradtarg1[FA3(2, j, i, 2, nd_v)]);
                gradtarg[FA3(j, 2, i, nd_v, 2)] =
                    -cimag(gradtarg1[FA3(1, j, i, 2, nd_v)])
                    - I * cimag(gradtarg1[FA3(2, j, i, 2, nd_v)]);
            }
        }
    }
    if (ifpghtarg_v == 3) {
        for (i = 1; i <= nt_v; i++) {
            for (j = 1; j <= nd_v; j++) {
                hesstarg[FA3(j, 1, i, nd_v, 3)] =
                    creal(hesstarg1[FA3(1, j, i, 2, nd_v)])
                    + I * creal(hesstarg1[FA3(2, j, i, 2, nd_v)]);
                hesstarg[FA3(j, 2, i, nd_v, 3)] =
                    -cimag(hesstarg1[FA3(1, j, i, 2, nd_v)])
                    - I * cimag(hesstarg1[FA3(2, j, i, 2, nd_v)]);
                hesstarg[FA3(j, 3, i, nd_v, 3)] =
                    -hesstarg[FA3(j, 1, i, nd_v, 3)];
            }
        }
    }

    free(charge1);
    free(dipstr1);
    free(pot1);
    free(grad1);
    free(hess1);
    free(pottarg1);
    free(gradtarg1);
    free(hesstarg1);
}
