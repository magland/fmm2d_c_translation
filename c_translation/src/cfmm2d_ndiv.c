/*
 * cfmm2d_ndiv.c - C translation of src/laplace/cfmm2d_ndiv.f
 *
 * Top-level 2D Laplace (Cauchy form) FMM driver with externally
 * supplied subdivision criterion. Same body as cfmm2d, except:
 *   - it does NOT call lndiv2d; ndiv and idivflag are INPUT parameters
 *   - ifnear, timeinfo, ier are pass-through parameters forwarded to
 *     cfmm2dmain
 *
 * Translated 1:1 from the Fortran reference (same control flow, same
 * allocations, same operation order). OpenMP, error/print logging, and
 * timing (cpu_time / omp_get_wtime) are stripped. time1/time2 are left
 * at zero so the timeinfo output from cfmm2dmain becomes deterministic
 * zero values. The diff test driver does not compare timeinfo (wall-
 * clock noise by design).
 *
 * Cross-file dependencies (pts_tree_*, dreorder*, l2dterms,
 * cfmm2dmain, l2dmpalloc) are issued via bare Fortran symbol names so
 * the diff test isolates this file. cfmm2dmain_ and l2dmpalloc_ are
 * resolved from cfmm2d.c's drop-in build or libfmm2d.a's reference
 * build.
 */

#include "cfmm2d_ndiv.h"

/* Cross-file calls: bare Fortran symbol names. */

/* pts_tree2d. */
extern void pts_tree_mem_(const double *src, const fint *ns,
                          const double *targ, const fint *nt,
                          const fint *idivflag, const fint *ndiv,
                          const fint *nlmin, const fint *nlmax,
                          const fint *ifunif, const fint *iper,
                          fint *nlevels, fint *nboxes, fint *ltree);

extern void pts_tree_build_(const double *src, const fint *ns,
                            const double *targ, const fint *nt,
                            const fint *idivflag, const fint *ndiv,
                            const fint *nlmin, const fint *nlmax,
                            const fint *ifunif, const fint *iper,
                            fint *nlevels, const fint *nboxes,
                            const fint *ltree, fint *itree, fint *iptr,
                            double *centers, double *boxsize);

extern void pts_tree_sort_(const fint *n, const double *xys,
                           const fint *itree, const fint *ltree,
                           const fint *nboxes, const fint *nlevels,
                           const fint *iptr, const double *centers,
                           fint *ixy, fint *ixyse);

/* fmmcommon2d. */
extern void dreorderf_(const fint *ndim, const fint *n, const double *arr,
                       double *arrsort, const fint *iarr);
extern void dreorderi_(const fint *ndim, const fint *n, const double *arr,
                       double *arrsort, const fint *iarr);

/* l2dterms. */
extern void l2dterms_(const double *eps, fint *nterms, fint *ier);

/* cfmm2dmain and l2dmpalloc - defined in cfmm2d.c (same-file FNAME
 * dispatch there). Here they are cross-file, so use bare names. */
extern void cfmm2dmain_(const fint *nd, const double *eps,
                        const fint *nsource, const double *sourcesort,
                        const fint *ifcharge, const fcomplex *chargesort,
                        const fint *ifdipole, const fcomplex *dipstrsort,
                        const fint *ntarget, const double *targetsort,
                        const fint *nexpc, const double *expcsort,
                        const fint *iaddr, double *rmlexp,
                        fcomplex *mptemp, const fint *lmptmp,
                        const fint *itree, const fint *ltree,
                        const fint *iptr, const fint *ndiv,
                        const fint *nlevels, const fint *nboxes,
                        const fint *iper, const double *boxsize,
                        const double *rscales, const double *centers,
                        const fint *laddr,
                        const fint *isrcse, const fint *itargse,
                        const fint *iexpcse, const fint *nterms,
                        const fint *ntj,
                        const fint *ifpgh, fcomplex *pot,
                        fcomplex *grad, fcomplex *hess,
                        const fint *ifpghtarg, fcomplex *pottarg,
                        fcomplex *gradtarg, fcomplex *hesstarg,
                        fcomplex *jsort, double *scjsort,
                        const fint *ifnear, double *timeinfo, fint *ier);

extern void l2dmpalloc_(const fint *nd, const fint *laddr, fint *iaddr,
                        const fint *nlevels, fint *lmptot,
                        const fint *nterms);


/* ---------------------------------------------------------------- */
/* cfmm2d_ndiv - top-level user-callable driver with external ndiv  */
/* ---------------------------------------------------------------- */
void FNAME(cfmm2d_ndiv)(const fint *nd, const double *eps,
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
                        const fint *ifnear, double *timeinfo, fint *ier)
{
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    fint ifcharge_v = *ifcharge;
    fint ifdipole_v = *ifdipole;
    fint ifpgh_v = *ifpgh;
    fint ifpghtarg_v = *ifpghtarg;

    fint *itree = NULL;
    fint iptr[8];
    fint nlmin, nlmax, ifunif;
    double *tcenters = NULL;
    double *boxsize = NULL;
    fint nexpc, ntj;
    double expc[2];
    double scj;
    fcomplex jexps[100];
    fint nlevels, nboxes;
    fint ltree;

    fint *isrc = NULL, *isrcse = NULL;
    fint *itarg = NULL, *itargse = NULL, *iexpcse = NULL;
    double *sourcesort = NULL;
    double *targsort = NULL;
    fcomplex *chargesort = NULL, *dipstrsort = NULL;
    fcomplex *potsort = NULL, *gradsort = NULL, *hesssort = NULL;
    fcomplex *pottargsort = NULL, *gradtargsort = NULL, *hesstargsort = NULL;

    fint lmptot;
    double *rscales = NULL;
    fint *nterms = NULL;
    fint *iaddr = NULL;
    double *rmlexp = NULL;
    fcomplex *mptemp = NULL;

    fint i, lmptmp, nmax, idim;

    /* Suppress unused. */
    (void)eps;

    nexpc = 0;
    nlevels = 0;
    nboxes = 0;
    ntj = 0;
    scj = 0.0;
    expc[0] = 0.0;
    expc[1] = 0.0;
    for (i = 0; i < 100; i++) jexps[i] = 0.0;

    /* No lndiv2d call: ndiv and idivflag are inputs. */

    ltree = 0;
    nlmin = 0;
    nlmax = 51;
    ifunif = 0;
    *iper = 0;

    /* initialize timeinfo */
    for (i = 0; i < 8; i++) {
        timeinfo[i] = 0.0;
    }

    /* tree memory sizing */
    pts_tree_mem_(sources, ns, targ, nt, idivflag, ndiv, &nlmin, &nlmax,
                  &ifunif, iper, &nlevels, &nboxes, &ltree);

    itree = (fint *)malloc(ltree * sizeof(fint));
    boxsize = (double *)malloc((nlevels + 1) * sizeof(double));
    tcenters = (double *)malloc(2 * nboxes * sizeof(double));

    /* build the tree */
    pts_tree_build_(sources, ns, targ, nt, idivflag, ndiv, &nlmin, &nlmax,
                    &ifunif, iper, &nlevels, &nboxes, &ltree, itree, iptr,
                    tcenters, boxsize);

    isrc = (fint *)malloc(ns_v * sizeof(fint));
    isrcse = (fint *)malloc(2 * nboxes * sizeof(fint));
    itarg = (fint *)malloc(nt_v * sizeof(fint));
    itargse = (fint *)malloc(2 * nboxes * sizeof(fint));
    iexpcse = (fint *)malloc(2 * nboxes * sizeof(fint));

    for (i = 1; i <= nboxes; i++) {
        iexpcse[FA2(1, i, 2)] = 1;
        iexpcse[FA2(2, i, 2)] = 0;
    }

    pts_tree_sort_(ns, sources, itree, &ltree, &nboxes, &nlevels, iptr,
                   tcenters, isrc, isrcse);

    pts_tree_sort_(nt, targ, itree, &ltree, &nboxes, &nlevels, iptr,
                   tcenters, itarg, itargse);

    sourcesort = (double *)malloc(2 * ns_v * sizeof(double));
    targsort = (double *)malloc(2 * nt_v * sizeof(double));

    if (ifcharge_v == 1 && ifdipole_v == 0) {
        chargesort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        dipstrsort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
    }
    if (ifcharge_v == 0 && ifdipole_v == 1) {
        chargesort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
        dipstrsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
    }
    if (ifcharge_v == 1 && ifdipole_v == 1) {
        chargesort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        dipstrsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
    }
    /* If both ifcharge=0 and ifdipole=0, allocate placeholders so the
     * pointers passed to cfmm2dmain are valid. */
    if (chargesort == NULL) {
        chargesort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
    }
    if (dipstrsort == NULL) {
        dipstrsort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
    }

    if (ifpgh_v == 1) {
        potsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        gradsort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
        hesssort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
    } else if (ifpgh_v == 2) {
        potsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        gradsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        hesssort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
    } else if (ifpgh_v == 3) {
        potsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        gradsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        hesssort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
    } else {
        potsort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
        gradsort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
        hesssort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
    }

    if (ifpghtarg_v == 1) {
        pottargsort = (fcomplex *)malloc(nd_v * nt_v * sizeof(fcomplex));
        gradtargsort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
        hesstargsort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
    } else if (ifpghtarg_v == 2) {
        pottargsort = (fcomplex *)malloc(nd_v * nt_v * sizeof(fcomplex));
        gradtargsort = (fcomplex *)malloc(nd_v * nt_v * sizeof(fcomplex));
        hesstargsort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
    } else if (ifpghtarg_v == 3) {
        pottargsort = (fcomplex *)malloc(nd_v * nt_v * sizeof(fcomplex));
        gradtargsort = (fcomplex *)malloc(nd_v * nt_v * sizeof(fcomplex));
        hesstargsort = (fcomplex *)malloc(nd_v * nt_v * sizeof(fcomplex));
    } else {
        pottargsort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
        gradtargsort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
        hesstargsort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
    }

    /* initialize potentials, gradients, hessians */
    if (ifpgh_v == 1) {
        for (i = 1; i <= ns_v; i++) {
            for (idim = 1; idim <= nd_v; idim++) {
                potsort[FA2(idim, i, nd_v)] = 0.0;
            }
        }
    }
    if (ifpgh_v == 2) {
        for (i = 1; i <= ns_v; i++) {
            for (idim = 1; idim <= nd_v; idim++) {
                potsort[FA2(idim, i, nd_v)] = 0.0;
                gradsort[FA2(idim, i, nd_v)] = 0.0;
            }
        }
    }
    if (ifpgh_v == 3) {
        for (i = 1; i <= ns_v; i++) {
            for (idim = 1; idim <= nd_v; idim++) {
                potsort[FA2(idim, i, nd_v)] = 0.0;
                gradsort[FA2(idim, i, nd_v)] = 0.0;
                hesssort[FA2(idim, i, nd_v)] = 0.0;
            }
        }
    }

    if (ifpghtarg_v == 1) {
        for (i = 1; i <= nt_v; i++) {
            for (idim = 1; idim <= nd_v; idim++) {
                pottarg[FA2(idim, i, nd_v)] = 0.0;
                pottargsort[FA2(idim, i, nd_v)] = 0.0;
            }
        }
    }
    if (ifpghtarg_v == 2) {
        for (i = 1; i <= nt_v; i++) {
            for (idim = 1; idim <= nd_v; idim++) {
                pottargsort[FA2(idim, i, nd_v)] = 0.0;
                gradtargsort[FA2(idim, i, nd_v)] = 0.0;
                gradtargsort[FA2(idim, i, nd_v)] = 0.0;
            }
        }
    }
    if (ifpghtarg_v == 3) {
        for (i = 1; i <= nt_v; i++) {
            for (idim = 1; idim <= nd_v; idim++) {
                pottargsort[FA2(idim, i, nd_v)] = 0.0;
                gradtargsort[FA2(idim, i, nd_v)] = 0.0;
                hesstargsort[FA2(idim, i, nd_v)] = 0.0;
            }
        }
    }

    /* compute scaling factors and lengths of mp/local expansions */
    rscales = (double *)malloc((nlevels + 1) * sizeof(double));
    nterms = (fint *)malloc((nlevels + 1) * sizeof(fint));

    nmax = 0;
    *ier = 0;
    for (i = 0; i <= nlevels; i++) {
        rscales[i] = boxsize[i];
        l2dterms_(eps, &nterms[i], ier);
        /* nterms(i) = nterms(i) - no-op in source */
        if (nterms[i] > nmax) nmax = nterms[i];
    }

    /* allocate iaddr and temporary arrays */
    iaddr = (fint *)malloc(2 * nboxes * sizeof(fint));
    lmptmp = (nmax + 1) * nd_v;
    mptemp = (fcomplex *)malloc(lmptmp * sizeof(fcomplex));

    /* reorder sources */
    {
        fint two = 2;
        fint twond = 2 * nd_v;
        dreorderf_(&two, ns, sources, sourcesort, isrc);
        if (ifcharge_v == 1) {
            dreorderf_(&twond, ns, (const double *)charge,
                       (double *)chargesort, isrc);
        }
        if (ifdipole_v == 1) {
            dreorderf_(&twond, ns, (const double *)dipstr,
                       (double *)dipstrsort, isrc);
        }
        /* reorder targets */
        dreorderf_(&two, nt, targ, targsort, itarg);
    }

    /* allocate workspace for multipole / local expansions */
    {
        /* l2dmpalloc takes laddr (which lives inside itree at iptr(1)). */
        const fint *laddr_p = &itree[iptr[0] - 1];
        l2dmpalloc_(&nd_v, laddr_p, iaddr, &nlevels, &lmptot, nterms);
    }

    rmlexp = (double *)malloc(lmptot * sizeof(double));

    /* call the main fmm routine */
    {
        const fint *laddr_p = &itree[iptr[0] - 1];
        cfmm2dmain_(&nd_v, eps,
                    ns, sourcesort,
                    &ifcharge_v, chargesort,
                    &ifdipole_v, dipstrsort,
                    nt, targsort, &nexpc, expc,
                    iaddr, rmlexp, mptemp, &lmptmp,
                    itree, &ltree, iptr, ndiv, &nlevels,
                    &nboxes, iper, boxsize, rscales, tcenters,
                    laddr_p,
                    isrcse, itargse, iexpcse, nterms, &ntj,
                    &ifpgh_v, potsort, gradsort, hesssort,
                    &ifpghtarg_v, pottargsort, gradtargsort,
                    hesstargsort, jexps, &scj, ifnear, timeinfo, ier);
    }

    /* resort the output arrays in input order */
    {
        fint twond = 2 * nd_v;
        if (ifpgh_v == 1) {
            dreorderi_(&twond, ns, (const double *)potsort,
                       (double *)pot, isrc);
        }
        if (ifpgh_v == 2) {
            dreorderi_(&twond, ns, (const double *)potsort,
                       (double *)pot, isrc);
            dreorderi_(&twond, ns, (const double *)gradsort,
                       (double *)grad, isrc);
        }
        if (ifpgh_v == 3) {
            dreorderi_(&twond, ns, (const double *)potsort,
                       (double *)pot, isrc);
            dreorderi_(&twond, ns, (const double *)gradsort,
                       (double *)grad, isrc);
            dreorderi_(&twond, ns, (const double *)hesssort,
                       (double *)hess, isrc);
        }

        if (ifpghtarg_v == 1) {
            dreorderi_(&twond, nt, (const double *)pottargsort,
                       (double *)pottarg, itarg);
        }
        if (ifpghtarg_v == 2) {
            dreorderi_(&twond, nt, (const double *)pottargsort,
                       (double *)pottarg, itarg);
            dreorderi_(&twond, nt, (const double *)gradtargsort,
                       (double *)gradtarg, itarg);
        }
        if (ifpghtarg_v == 3) {
            dreorderi_(&twond, nt, (const double *)pottargsort,
                       (double *)pottarg, itarg);
            dreorderi_(&twond, nt, (const double *)gradtargsort,
                       (double *)gradtarg, itarg);
            dreorderi_(&twond, nt, (const double *)hesstargsort,
                       (double *)hesstarg, itarg);
        }
    }

    free(itree);
    free(boxsize);
    free(tcenters);
    free(isrc);
    free(isrcse);
    free(itarg);
    free(itargse);
    free(iexpcse);
    free(sourcesort);
    free(targsort);
    free(chargesort);
    free(dipstrsort);
    free(potsort);
    free(gradsort);
    free(hesssort);
    free(pottargsort);
    free(gradtargsort);
    free(hesstargsort);
    free(rscales);
    free(nterms);
    free(iaddr);
    free(mptemp);
    free(rmlexp);
}
