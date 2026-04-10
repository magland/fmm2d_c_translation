/*
 * h2dterms.c - C translation of src/helmholtz/h2dterms.f
 *
 * Only h2dterms is translated (the only routine reachable from hfmm2d).
 * h2dterms2, h2dterms_far, h2dterms_list2, h2dterms_list2e,
 * h2dterms_eval are not translated.
 *
 * Strict 1:1 translation. Uses implicit real*8 (a-h,o-z).
 */

#include "h2dterms.h"

/* Cross-file dependencies */
extern void h2dall_(const fint *nterms, const fcomplex *z,
                    const double *rscale, fcomplex *hvec,
                    const fint *ifder, fcomplex *hder);

extern void jbessel2d_(const fint *nterms, const fcomplex *z,
                       const double *rscale, fcomplex *fjs,
                       const fint *ifder, fcomplex *fjder);

void FNAME(h2dterms)(const double *bsize_p, const fcomplex *zk_p,
                     const double *eps_p, fint *nterms_p, fint *ier_p)
{
    double bsize_v = *bsize_p;
    fcomplex zk_v = *zk_p;
    double eps_v = *eps_p;

    double pi2;
    fcomplex z1, z2;
    fint ntmax, ifder, j;
    double rscale;
    double xtemp0, xtemp1, xtemp2, xtemp;
    fcomplex fjder_dum[2], fhder_dum[2];

    fcomplex *jfun, *hfun;

    *ier_p = 0;
    pi2 = 8 * atan(1.0);

    z1 = (zk_v * bsize_v) * 1.5;

    ntmax = 50000;
    jfun = (fcomplex *)malloc((ntmax + 101) * sizeof(fcomplex));
    hfun = (fcomplex *)malloc((ntmax + 101) * sizeof(fcomplex));
    ifder = 0;
    rscale = 1.0;
    if (cabs(zk_v * bsize_v) < pi2) rscale = cabs(zk_v * bsize_v);
    h2dall_(&ntmax, &z1, &rscale, hfun, &ifder, fhder_dum);

    z2 = (zk_v * bsize_v) * sqrt(2.0) / 2.0;

    jbessel2d_(&ntmax, &z2, &rscale, jfun, &ifder, fjder_dum);

    xtemp1 = cabs(jfun[0] * hfun[0]);
    xtemp2 = cabs(jfun[1] * hfun[1]);
    xtemp0 = xtemp1 + xtemp2;
    *nterms_p = 1;
    for (j = 2; j <= ntmax; j++) {
        xtemp1 = cabs(jfun[j] * hfun[j]);
        xtemp2 = cabs(jfun[j - 1] * hfun[j - 1]);
        xtemp = (xtemp1 + xtemp2) * cabs(hfun[0]);
        if (xtemp < eps_v * xtemp0) {
            *nterms_p = j + 1;
            free(jfun);
            free(hfun);
            return;
        }
    }

    /* computational box is too big */
    *ier_p = 13;
    *nterms_p = 10001;

    free(jfun);
    free(hfun);
    return;
}
