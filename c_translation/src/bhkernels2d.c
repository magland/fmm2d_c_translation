/*
 * bhkernels2d.c - C translation of src/biharmonic/bhkernels2d.f
 *
 * Direct evaluation kernels for the 2D biharmonic FMM. Each routine
 * INCREMENTS its output arrays (vel, grad) at every target due to
 * charges and/or dipoles at every source.
 *
 * Threshold check: cabs(zdis) .le. thresh (absolute value of the
 * complex distance compared against thresh, using .le.). This is
 * preserved exactly from the Fortran source.
 *
 * Multi-term assignment statements in the Fortran source such as
 *
 *     vel(idim,j) = vel(idim,j) + A + B
 *
 * are evaluated left-to-right as ((vel + A) + B). Writing this as the
 * C compound `vel += A + B` instead parses as `vel + (A + B)`, giving
 * a different rounding order at -O0. To keep bit-for-bit equality with
 * the Fortran reference every multi-term update is split into separate
 * += statements in the same left-to-right order as the Fortran source.
 */

#include "bhkernels2d.h"

/*
 * bh2d_directcp: charges -> velocity
 *
 * vel += 2*charges(1)*log|zdis| + charges(2)*conj(1/zdis)*zdis
 */
void FNAME(bh2d_directcp)(const fint *nd, const double *sources,
                          const fint *ns, const fcomplex *charges,
                          const double *targ, const fint *nt,
                          fcomplex *vel, const double *thresh)
{
    fint i, j, idim;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    fcomplex zs, zt, zdis, zdis1;
    double logabs;
    /* unused locals preserved from the Fortran declarations */
    fcomplex zdis2, eye;
    (void)zdis2;
    (void)eye;

    for (j = 1; j <= nt_v; j++) {
        zt = targ[FA2(1, j, 2)] + targ[FA2(2, j, 2)] * I;
        for (i = 1; i <= ns_v; i++) {
            zs = sources[FA2(1, i, 2)] + sources[FA2(2, i, 2)] * I;
            zdis = zt - zs;
            if (cabs(zdis) <= *thresh) goto skip;
            zdis1 = 1.0 / zdis;
            logabs = log(cabs(zdis));
            for (idim = 1; idim <= nd_v; idim++) {
                /* Fortran:
                 *   vel(idim,j) = vel(idim,j)
                 *       + 2*charges(idim,1,i)*log(cdabs(zdis))
                 *       + charges(idim,2,i)*dconjg(zdis1)*zdis
                 * Left-to-right: split into two += statements. */
                vel[FA2(idim, j, nd_v)] +=
                    2.0 * charges[FA3(idim, 1, i, nd_v, 2)] * logabs;
                vel[FA2(idim, j, nd_v)] +=
                    charges[FA3(idim, 2, i, nd_v, 2)] * conj(zdis1) * zdis;
            }
        skip:;
        }
    }
}

/*
 * bh2d_directcg: charges -> velocity + gradient
 *
 * vel += 2*charges(1)*log|zdis| + charges(2)*conj(1/zdis)*zdis
 * grad(1) += charges(1) * (1/zdis)
 * grad(2) += charges(2) * conj(1/zdis)
 * grad(3) += charges(1) * conj(1/zdis)
 * grad(3) -= charges(2) * conj((1/zdis)^2) * zdis
 */
void FNAME(bh2d_directcg)(const fint *nd, const double *sources,
                          const fint *ns, const fcomplex *charges,
                          const double *targ, const fint *nt,
                          fcomplex *vel, fcomplex *grad,
                          const double *thresh)
{
    fint i, j, idim;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    fcomplex zs, zt, zdis, zdis1, zdis2;
    double logabs;
    fcomplex eye;
    (void)eye;

    for (j = 1; j <= nt_v; j++) {
        zt = targ[FA2(1, j, 2)] + targ[FA2(2, j, 2)] * I;
        for (i = 1; i <= ns_v; i++) {
            zs = sources[FA2(1, i, 2)] + sources[FA2(2, i, 2)] * I;
            zdis = zt - zs;
            if (cabs(zdis) <= *thresh) goto skip;
            zdis1 = 1.0 / zdis;
            zdis2 = zdis1 * zdis1;
            logabs = log(cabs(zdis));
            for (idim = 1; idim <= nd_v; idim++) {
                /* vel(idim,j) = vel(idim,j)
                 *     + 2*charges(idim,1,i)*log(cdabs(zdis))
                 *     + charges(idim,2,i)*dconjg(zdis1)*zdis  */
                vel[FA2(idim, j, nd_v)] +=
                    2.0 * charges[FA3(idim, 1, i, nd_v, 2)] * logabs;
                vel[FA2(idim, j, nd_v)] +=
                    charges[FA3(idim, 2, i, nd_v, 2)] * conj(zdis1) * zdis;

                /* grad(idim,1,j) = grad(idim,1,j)
                 *     + charges(idim,1,i)*zdis1 */
                grad[FA3(idim, 1, j, nd_v, 3)] +=
                    charges[FA3(idim, 1, i, nd_v, 2)] * zdis1;

                /* grad(idim,2,j) = grad(idim,2,j)
                 *     + charges(idim,2,i)*dconjg(zdis1) */
                grad[FA3(idim, 2, j, nd_v, 3)] +=
                    charges[FA3(idim, 2, i, nd_v, 2)] * conj(zdis1);

                /* grad(idim,3,j) = grad(idim,3,j)
                 *     + charges(idim,1,i)*dconjg(zdis1)
                 * grad(idim,3,j) = grad(idim,3,j)
                 *     - charges(idim,2,i)*dconjg(zdis2)*zdis
                 * (two separate Fortran statements). */
                grad[FA3(idim, 3, j, nd_v, 3)] +=
                    charges[FA3(idim, 1, i, nd_v, 2)] * conj(zdis1);
                grad[FA3(idim, 3, j, nd_v, 3)] -=
                    charges[FA3(idim, 2, i, nd_v, 2)] * conj(zdis2) * zdis;
            }
        skip:;
        }
    }
}

/*
 * bh2d_directdp: dipoles -> velocity
 *
 * vel += dippar(1)*(1/zdis) + dippar(3)*conj(1/zdis)
 * vel += dippar(2)*conj((1/zdis)^2)*zdis
 */
void FNAME(bh2d_directdp)(const fint *nd, const double *sources,
                          const fint *ns, const fcomplex *dippar,
                          const double *targ, const fint *nt,
                          fcomplex *vel, const double *thresh)
{
    fint i, j, idim;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    fcomplex zs, zt, zdis, zdis1, zdis2;
    fcomplex eye;
    (void)eye;

    for (j = 1; j <= nt_v; j++) {
        zt = targ[FA2(1, j, 2)] + targ[FA2(2, j, 2)] * I;
        for (i = 1; i <= ns_v; i++) {
            zs = sources[FA2(1, i, 2)] + sources[FA2(2, i, 2)] * I;
            zdis = zt - zs;
            if (cabs(zdis) <= *thresh) goto skip;
            zdis1 = 1.0 / zdis;
            zdis2 = zdis1 * zdis1;
            for (idim = 1; idim <= nd_v; idim++) {
                /* vel(idim,j) = vel(idim,j)
                 *     + dippar(idim,1,i)*zdis1
                 *     + dippar(idim,3,i)*dconjg(zdis1) */
                vel[FA2(idim, j, nd_v)] +=
                    dippar[FA3(idim, 1, i, nd_v, 3)] * zdis1;
                vel[FA2(idim, j, nd_v)] +=
                    dippar[FA3(idim, 3, i, nd_v, 3)] * conj(zdis1);

                /* vel(idim,j) = vel(idim,j)
                 *     + dippar(idim,2,i)*dconjg(zdis2)*zdis */
                vel[FA2(idim, j, nd_v)] +=
                    dippar[FA3(idim, 2, i, nd_v, 3)] * conj(zdis2) * zdis;
            }
        skip:;
        }
    }
}

/*
 * bh2d_directdg: dipoles -> velocity + gradient
 *
 * vel += dippar(1)*(1/zdis) + dippar(3)*conj(1/zdis)
 * vel += dippar(2)*conj((1/zdis)^2)*zdis
 * grad(1) -= dippar(1) * (1/zdis)^2
 * grad(2) += dippar(2) * conj((1/zdis)^2)
 * grad(3) -= dippar(3) * conj((1/zdis)^2)
 * grad(3) -= 2 * dippar(2) * conj((1/zdis)^2 * (1/zdis)) * zdis
 */
void FNAME(bh2d_directdg)(const fint *nd, const double *sources,
                          const fint *ns, const fcomplex *dippar,
                          const double *targ, const fint *nt,
                          fcomplex *vel, fcomplex *grad,
                          const double *thresh)
{
    fint i, j, idim;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    fcomplex zs, zt, zdis, zdis1, zdis2;
    fcomplex eye;
    (void)eye;

    for (j = 1; j <= nt_v; j++) {
        zt = targ[FA2(1, j, 2)] + targ[FA2(2, j, 2)] * I;
        for (i = 1; i <= ns_v; i++) {
            zs = sources[FA2(1, i, 2)] + sources[FA2(2, i, 2)] * I;
            zdis = zt - zs;
            if (cabs(zdis) <= *thresh) goto skip;
            zdis1 = 1.0 / zdis;
            zdis2 = zdis1 * zdis1;
            for (idim = 1; idim <= nd_v; idim++) {
                /* vel(idim,j) = vel(idim,j)
                 *     + dippar(idim,1,i)*zdis1
                 *     + dippar(idim,3,i)*dconjg(zdis1) */
                vel[FA2(idim, j, nd_v)] +=
                    dippar[FA3(idim, 1, i, nd_v, 3)] * zdis1;
                vel[FA2(idim, j, nd_v)] +=
                    dippar[FA3(idim, 3, i, nd_v, 3)] * conj(zdis1);

                /* vel(idim,j) = vel(idim,j)
                 *     + dippar(idim,2,i)*dconjg(zdis2)*zdis */
                vel[FA2(idim, j, nd_v)] +=
                    dippar[FA3(idim, 2, i, nd_v, 3)] * conj(zdis2) * zdis;

                /* grad(idim,1,j) = grad(idim,1,j)
                 *     - dippar(idim,1,i)*(zdis2) */
                grad[FA3(idim, 1, j, nd_v, 3)] -=
                    dippar[FA3(idim, 1, i, nd_v, 3)] * zdis2;

                /* grad(idim,2,j) = grad(idim,2,j)
                 *     + dippar(idim,2,i)*dconjg(zdis2) */
                grad[FA3(idim, 2, j, nd_v, 3)] +=
                    dippar[FA3(idim, 2, i, nd_v, 3)] * conj(zdis2);

                /* grad(idim,3,j) = grad(idim,3,j)
                 *     - dippar(idim,3,i)*dconjg(zdis2)
                 * grad(idim,3,j) = grad(idim,3,j)
                 *     - 2*dippar(idim,2,i)*dconjg(zdis2*zdis1)*zdis
                 * (two separate Fortran statements). */
                grad[FA3(idim, 3, j, nd_v, 3)] -=
                    dippar[FA3(idim, 3, i, nd_v, 3)] * conj(zdis2);
                grad[FA3(idim, 3, j, nd_v, 3)] -=
                    2.0 * dippar[FA3(idim, 2, i, nd_v, 3)]
                        * conj(zdis2 * zdis1) * zdis;
            }
        skip:;
        }
    }
}

/*
 * bh2d_directcdp: charges + dipoles -> velocity
 */
void FNAME(bh2d_directcdp)(const fint *nd, const double *sources,
                           const fint *ns, const fcomplex *charges,
                           const fcomplex *dippar,
                           const double *targ, const fint *nt,
                           fcomplex *vel, const double *thresh)
{
    fint i, j, idim;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    fcomplex zs, zt, zdis, zdis1, zdis2;
    double logabs;
    fcomplex eye;
    (void)eye;

    for (j = 1; j <= nt_v; j++) {
        zt = targ[FA2(1, j, 2)] + targ[FA2(2, j, 2)] * I;
        for (i = 1; i <= ns_v; i++) {
            zs = sources[FA2(1, i, 2)] + sources[FA2(2, i, 2)] * I;
            zdis = zt - zs;
            if (cabs(zdis) <= *thresh) goto skip;
            zdis1 = 1.0 / zdis;
            zdis2 = zdis1 * zdis1;
            logabs = log(cabs(zdis));
            for (idim = 1; idim <= nd_v; idim++) {
                /* vel(idim,j) = vel(idim,j)
                 *     + 2*charges(idim,1,i)*log(cdabs(zdis))
                 *     + charges(idim,2,i)*dconjg(zdis1)*zdis */
                vel[FA2(idim, j, nd_v)] +=
                    2.0 * charges[FA3(idim, 1, i, nd_v, 2)] * logabs;
                vel[FA2(idim, j, nd_v)] +=
                    charges[FA3(idim, 2, i, nd_v, 2)] * conj(zdis1) * zdis;

                /* vel(idim,j) = vel(idim,j)
                 *     + dippar(idim,1,i)*zdis1
                 *     + dippar(idim,3,i)*dconjg(zdis1) */
                vel[FA2(idim, j, nd_v)] +=
                    dippar[FA3(idim, 1, i, nd_v, 3)] * zdis1;
                vel[FA2(idim, j, nd_v)] +=
                    dippar[FA3(idim, 3, i, nd_v, 3)] * conj(zdis1);

                /* vel(idim,j) = vel(idim,j)
                 *     + dippar(idim,2,i)*dconjg(zdis2)*zdis */
                vel[FA2(idim, j, nd_v)] +=
                    dippar[FA3(idim, 2, i, nd_v, 3)] * conj(zdis2) * zdis;
            }
        skip:;
        }
    }
}

/*
 * bh2d_directcdg: charges + dipoles -> velocity + gradient
 */
void FNAME(bh2d_directcdg)(const fint *nd, const double *sources,
                           const fint *ns, const fcomplex *charges,
                           const fcomplex *dippar,
                           const double *targ, const fint *nt,
                           fcomplex *vel, fcomplex *grad,
                           const double *thresh)
{
    fint i, j, idim;
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    fcomplex zs, zt, zdis, zdis1, zdis2;
    double logabs;
    fcomplex eye;
    (void)eye;

    for (j = 1; j <= nt_v; j++) {
        zt = targ[FA2(1, j, 2)] + targ[FA2(2, j, 2)] * I;
        for (i = 1; i <= ns_v; i++) {
            zs = sources[FA2(1, i, 2)] + sources[FA2(2, i, 2)] * I;
            zdis = zt - zs;
            if (cabs(zdis) <= *thresh) goto skip;
            zdis1 = 1.0 / zdis;
            zdis2 = zdis1 * zdis1;
            logabs = log(cabs(zdis));
            for (idim = 1; idim <= nd_v; idim++) {
                /* vel(idim,j) = vel(idim,j)
                 *     + 2*charges(idim,1,i)*log(cdabs(zdis))
                 *     + charges(idim,2,i)*dconjg(zdis1)*zdis */
                vel[FA2(idim, j, nd_v)] +=
                    2.0 * charges[FA3(idim, 1, i, nd_v, 2)] * logabs;
                vel[FA2(idim, j, nd_v)] +=
                    charges[FA3(idim, 2, i, nd_v, 2)] * conj(zdis1) * zdis;

                /* vel(idim,j) = vel(idim,j)
                 *     + dippar(idim,1,i)*zdis1
                 *     + dippar(idim,3,i)*dconjg(zdis1) */
                vel[FA2(idim, j, nd_v)] +=
                    dippar[FA3(idim, 1, i, nd_v, 3)] * zdis1;
                vel[FA2(idim, j, nd_v)] +=
                    dippar[FA3(idim, 3, i, nd_v, 3)] * conj(zdis1);

                /* vel(idim,j) = vel(idim,j)
                 *     + dippar(idim,2,i)*dconjg(zdis2)*zdis */
                vel[FA2(idim, j, nd_v)] +=
                    dippar[FA3(idim, 2, i, nd_v, 3)] * conj(zdis2) * zdis;

                /* grad(idim,1,j) = grad(idim,1,j)
                 *     + charges(idim,1,i)*zdis1
                 * grad(idim,1,j) = grad(idim,1,j)
                 *     - dippar(idim,1,i)*(zdis2)   */
                grad[FA3(idim, 1, j, nd_v, 3)] +=
                    charges[FA3(idim, 1, i, nd_v, 2)] * zdis1;
                grad[FA3(idim, 1, j, nd_v, 3)] -=
                    dippar[FA3(idim, 1, i, nd_v, 3)] * zdis2;

                /* grad(idim,2,j) = grad(idim,2,j)
                 *     + charges(idim,2,i)*dconjg(zdis1)
                 * grad(idim,2,j) = grad(idim,2,j)
                 *     + dippar(idim,2,i)*dconjg(zdis2) */
                grad[FA3(idim, 2, j, nd_v, 3)] +=
                    charges[FA3(idim, 2, i, nd_v, 2)] * conj(zdis1);
                grad[FA3(idim, 2, j, nd_v, 3)] +=
                    dippar[FA3(idim, 2, i, nd_v, 3)] * conj(zdis2);

                /* grad(idim,3,j) = grad(idim,3,j)
                 *     + charges(idim,1,i)*dconjg(zdis1)
                 * grad(idim,3,j) = grad(idim,3,j)
                 *     - charges(idim,2,i)*dconjg(zdis2)*zdis
                 * grad(idim,3,j) = grad(idim,3,j)
                 *     - dippar(idim,3,i)*dconjg(zdis2)
                 * grad(idim,3,j) = grad(idim,3,j)
                 *     - 2*dippar(idim,2,i)*dconjg(zdis2*zdis1)*zdis */
                grad[FA3(idim, 3, j, nd_v, 3)] +=
                    charges[FA3(idim, 1, i, nd_v, 2)] * conj(zdis1);
                grad[FA3(idim, 3, j, nd_v, 3)] -=
                    charges[FA3(idim, 2, i, nd_v, 2)] * conj(zdis2) * zdis;
                grad[FA3(idim, 3, j, nd_v, 3)] -=
                    dippar[FA3(idim, 3, i, nd_v, 3)] * conj(zdis2);
                grad[FA3(idim, 3, j, nd_v, 3)] -=
                    2.0 * dippar[FA3(idim, 2, i, nd_v, 3)]
                        * conj(zdis2 * zdis1) * zdis;
            }
        skip:;
        }
    }
}
