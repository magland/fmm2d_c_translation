/*
 * stokkernels2d.c - C translation of src/stokes/stokkernels2d.f
 *
 * Direct evaluation kernels for the 2D Stokes FMM.
 * Each routine INCREMENTS its output arrays.
 */

#include "stokkernels2d.h"

/*
 * st2ddirectstokg: Stokeslet -> velocity + pressure + gradient
 * Threshold: r2 .lt. threshsq
 */
void FNAME(st2ddirectstokg)(const fint *nd, const double *sources,
                            const double *stoklet, const fint *ns,
                            const double *targ, const fint *nt,
                            double *pot, double *pre, double *grad,
                            const double *thresh)
{
    fint i, j, idim;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    double zdiff1, zdiff2, r2, r4, rtmp, pl, pl2;
    double threshsq;

    threshsq = (*thresh) * (*thresh);

    for (i = 1; i <= nt_v; i++) {
        for (j = 1; j <= ns_v; j++) {
            zdiff1 = targ[FA2(1, i, 2)] - sources[FA2(1, j, 2)];
            zdiff2 = targ[FA2(2, i, 2)] - sources[FA2(2, j, 2)];

            r2 = zdiff1 * zdiff1 + zdiff2 * zdiff2;
            if (r2 < threshsq) goto skip;

            rtmp = zdiff1 * zdiff1 - zdiff2 * zdiff2;
            r4 = r2 * r2;

            for (idim = 1; idim <= nd_v; idim++) {

                pot[FA3(idim, 1, i, nd_v, 2)] -=
                    stoklet[FA3(idim, 1, j, nd_v, 2)] * log(r2) / 4.0;
                pot[FA3(idim, 2, i, nd_v, 2)] -=
                    stoklet[FA3(idim, 2, j, nd_v, 2)] * log(r2) / 4.0;

                pl = zdiff1 * stoklet[FA3(idim, 1, j, nd_v, 2)]
                   + zdiff2 * stoklet[FA3(idim, 2, j, nd_v, 2)];

                pl2 = zdiff1 * stoklet[FA3(idim, 2, j, nd_v, 2)]
                    + zdiff2 * stoklet[FA3(idim, 1, j, nd_v, 2)];

                pot[FA3(idim, 1, i, nd_v, 2)] += zdiff1 * pl / r2 / 2.0;
                pot[FA3(idim, 2, i, nd_v, 2)] += zdiff2 * pl / r2 / 2.0;

                grad[FA4(idim, 1, 1, i, nd_v, 2, 2)] -= rtmp * pl / r4 / 2.0;
                grad[FA4(idim, 2, 2, i, nd_v, 2, 2)] += rtmp * pl / r4 / 2.0;
                grad[FA4(idim, 2, 1, i, nd_v, 2, 2)] +=
                    (pl2 * rtmp
                   - stoklet[FA3(idim, 1, j, nd_v, 2)] * 4.0 * zdiff2
                     * (zdiff1 * zdiff1)) / 2.0 / r4;

                grad[FA4(idim, 1, 2, i, nd_v, 2, 2)] -=
                    (pl2 * rtmp
                   + stoklet[FA3(idim, 2, j, nd_v, 2)] * 4.0 * zdiff1
                     * (zdiff2 * zdiff2)) / 2.0 / r4;

                pre[FA2(idim, i, nd_v)] += pl / r2;
            }
        skip:;
        }
    }
}

/*
 * st2ddirectstokstrsg: Stokeslet + stresslet -> velocity + pressure + gradient
 * Threshold: r2 .lt. threshsq
 */
void FNAME(st2ddirectstokstrsg)(const fint *nd, const double *sources,
                                const fint *ifstoklet, const double *stoklet,
                                const fint *istress, const double *strslet,
                                const double *strsvec,
                                const fint *ns, const double *targ, const fint *nt,
                                double *pot, double *pre, double *grad,
                                const double *thresh)
{
    fint i, j, idim;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    double zdiff1, zdiff2, r2, r4, r6, rtmp, pl, pl2, pv;
    double threshsq;

    if (*ifstoklet == 1) {
        FNAME(st2ddirectstokg)(nd, sources, stoklet, ns,
                               targ, nt, pot, pre, grad, thresh);
    }

    threshsq = (*thresh) * (*thresh);
    if (*istress != 1) goto label_100;

    /* stresslet contribution */
    for (i = 1; i <= nt_v; i++) {
        for (j = 1; j <= ns_v; j++) {
            zdiff1 = targ[FA2(1, i, 2)] - sources[FA2(1, j, 2)];
            zdiff2 = targ[FA2(2, i, 2)] - sources[FA2(2, j, 2)];

            r2 = zdiff1 * zdiff1 + zdiff2 * zdiff2;
            if (r2 < threshsq) goto skip;

            rtmp = zdiff1 * zdiff1 - zdiff2 * zdiff2;
            r4 = r2 * r2;
            r6 = r4 * r2;

            for (idim = 1; idim <= nd_v; idim++) {

                pl = zdiff1 * strslet[FA3(idim, 1, j, nd_v, 2)]
                   + zdiff2 * strslet[FA3(idim, 2, j, nd_v, 2)];

                pl2 = zdiff1 * strsvec[FA3(idim, 1, j, nd_v, 2)]
                    + zdiff2 * strsvec[FA3(idim, 2, j, nd_v, 2)];

                pv = strslet[FA3(idim, 1, j, nd_v, 2)]
                       * strsvec[FA3(idim, 1, j, nd_v, 2)]
                   + strslet[FA3(idim, 2, j, nd_v, 2)]
                       * strsvec[FA3(idim, 2, j, nd_v, 2)];

                pot[FA3(idim, 1, i, nd_v, 2)] -= 2 * zdiff1 * pl * pl2 / r4;
                pot[FA3(idim, 2, i, nd_v, 2)] -= 2 * zdiff2 * pl * pl2 / r4;

                /* grad(1,1) = grad(1,1) - A - B - C + D  (split) */
                grad[FA4(idim, 1, 1, i, nd_v, 2, 2)] -= 2 * pl * pl2 / r4;
                grad[FA4(idim, 1, 1, i, nd_v, 2, 2)] -=
                    2 * zdiff1 * strslet[FA3(idim, 1, j, nd_v, 2)] * pl2 / r4;
                grad[FA4(idim, 1, 1, i, nd_v, 2, 2)] -=
                    2 * zdiff1 * strsvec[FA3(idim, 1, j, nd_v, 2)] * pl / r4;
                grad[FA4(idim, 1, 1, i, nd_v, 2, 2)] +=
                    8 * zdiff1 * zdiff1 * pl * pl2 / r6;

                /* grad(2,1) = ... - A - B + C */
                grad[FA4(idim, 2, 1, i, nd_v, 2, 2)] -=
                    2 * zdiff1 * strslet[FA3(idim, 2, j, nd_v, 2)] * pl2 / r4;
                grad[FA4(idim, 2, 1, i, nd_v, 2, 2)] -=
                    2 * zdiff1 * strsvec[FA3(idim, 2, j, nd_v, 2)] * pl / r4;
                grad[FA4(idim, 2, 1, i, nd_v, 2, 2)] +=
                    8 * zdiff1 * zdiff2 * pl * pl2 / r6;

                /* grad(1,2) = ... - A - B + C */
                grad[FA4(idim, 1, 2, i, nd_v, 2, 2)] -=
                    2 * zdiff2 * strslet[FA3(idim, 1, j, nd_v, 2)] * pl2 / r4;
                grad[FA4(idim, 1, 2, i, nd_v, 2, 2)] -=
                    2 * zdiff2 * strsvec[FA3(idim, 1, j, nd_v, 2)] * pl / r4;
                grad[FA4(idim, 1, 2, i, nd_v, 2, 2)] +=
                    8 * zdiff1 * zdiff2 * pl * pl2 / r6;

                /* grad(2,2) = grad(2,2) - A - B - C + D */
                grad[FA4(idim, 2, 2, i, nd_v, 2, 2)] -= 2 * pl * pl2 / r4;
                grad[FA4(idim, 2, 2, i, nd_v, 2, 2)] -=
                    2 * zdiff2 * strslet[FA3(idim, 2, j, nd_v, 2)] * pl2 / r4;
                grad[FA4(idim, 2, 2, i, nd_v, 2, 2)] -=
                    2 * zdiff2 * strsvec[FA3(idim, 2, j, nd_v, 2)] * pl / r4;
                grad[FA4(idim, 2, 2, i, nd_v, 2, 2)] +=
                    8 * zdiff2 * zdiff2 * pl * pl2 / r6;

                /* pre = pre - 4*pl*pl2/r4 + 2*pv/r2 */
                pre[FA2(idim, i, nd_v)] -= 4 * pl * pl2 / r4;
                pre[FA2(idim, i, nd_v)] += 2 * pv / r2;
            }
        skip:;
        }
    }
label_100:;
}
