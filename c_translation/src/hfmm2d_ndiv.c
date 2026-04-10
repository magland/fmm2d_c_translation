/*
 * hfmm2d_ndiv.c - C translation of src/helmholtz/hfmm2d_ndiv.f
 *
 * MATLAB entry point for the Helmholtz 2D FMM. Structurally identical
 * to hfmm2d but takes ndiv/idivflag/ifnear as explicit parameters.
 *
 * Cross-file calls use bare Fortran symbol names. hfmm2dmain_ is the
 * Fortran version from libfmm2d.a (or its C drop-in).
 */

#include "hfmm2d_ndiv.h"

/* Cross-file calls */
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

extern void dreorderf_(const fint *ndim, const fint *n, const double *arr,
    double *arrsort, const fint *iarr);
extern void dreorderi_(const fint *ndim, const fint *n, const double *arr,
    double *arrsort, const fint *iarr);

extern void h2dterms_(const double *bsize, const fcomplex *zk,
    const double *eps, fint *nterms, fint *ier);

extern void h2dmpalloc_(const fint *nd, const fint *laddr, fint *iaddr,
    const fint *nlevels, fint *lmptot, const fint *nterms);

extern void hfmm2dmain_(const fint *nd, const double *eps,
    const fcomplex *zk, const fint *ns, const double *sourcesort,
    const fint *ifcharge, const fcomplex *chargesort,
    const fint *ifdipole, const fcomplex *dipstrsort,
    const double *dipvecsort,
    const fint *nt, const double *targsort, const fint *nexpc,
    const double *expc,
    const fint *iaddr, double *rmlexp, fcomplex *mptemp,
    const fint *lmptmp,
    const fint *itree, const fint *ltree, const fint *iptr,
    const fint *ndiv, const fint *nlevels,
    const fint *nboxes, const fint *iper, const double *boxsize,
    const double *rscales, const double *centers,
    const fint *laddr,
    const fint *isrcse, const fint *itargse, const fint *iexpcse,
    const fint *nterms, const fint *ntj,
    const fint *ifpgh, fcomplex *pot, fcomplex *grad, fcomplex *hess,
    const fint *ifpghtarg, fcomplex *pottarg, fcomplex *gradtarg,
    fcomplex *hesstarg, const fcomplex *jexps, const double *scj,
    const fint *ifnear, double *timeinfo, fint *ier);

void FNAME(hfmm2d_ndiv)(const fint *nd, const double *eps,
    const fcomplex *zk, const fint *ns, const double *sources,
    const fint *ifcharge, const fcomplex *charge,
    const fint *ifdipole, const fcomplex *dipstr, const double *dipvec,
    fint *iper, const fint *ifpgh, fcomplex *pot, fcomplex *grad,
    fcomplex *hess, const fint *nt, const double *targ,
    const fint *ifpghtarg, fcomplex *pottarg, fcomplex *gradtarg,
    fcomplex *hesstarg, const fint *ndiv, const fint *idivflag,
    const fint *ifnear, double *timeinfo, fint *ier)
{
    fint nd_v = *nd;
    fint ns_v = *ns;
    fint nt_v = *nt;
    fint ifcharge_v = *ifcharge;
    fint ifdipole_v = *ifdipole;
    fint ifpgh_v = *ifpgh;
    fint ifpghtarg_v = *ifpghtarg;
    fint ndiv_v = *ndiv;
    fint idivflag_v = *idivflag;
    fint ifnear_v = *ifnear;

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
    fint *itarg_arr = NULL, *itargse = NULL, *iexpcse = NULL;
    double *sourcesort = NULL;
    double *targsort = NULL;
    fcomplex *chargesort = NULL, *dipstrsort = NULL;
    double *dipvecsort = NULL;
    fcomplex *potsort = NULL, *gradsort = NULL, *hesssort = NULL;
    fcomplex *pottargsort = NULL, *gradtargsort = NULL, *hesstargsort = NULL;

    fint lmptot;
    double *rscales = NULL;
    fint *nterms = NULL;
    fint *iaddr = NULL;
    double *rmlexp = NULL;
    fcomplex *mptemp = NULL;

    fint i, lmptmp, nmax, idim;
    double pi;

    pi = 4.0 * atan(1.0);

    nexpc = 0;
    nlevels = 0;
    nboxes = 0;
    ntj = 0;
    scj = 0.0;
    expc[0] = 0.0;
    expc[1] = 0.0;
    for (i = 0; i < 100; i++) jexps[i] = 0.0;

    ltree = 0;
    nlmin = 0;
    nlmax = 51;
    ifunif = 0;
    *iper = 0;

    for (i = 0; i < 8; i++) timeinfo[i] = 0.0;

    pts_tree_mem_(sources, ns, targ, nt, &idivflag_v, &ndiv_v, &nlmin, &nlmax,
                  &ifunif, iper, &nlevels, &nboxes, &ltree);

    itree = (fint *)malloc(ltree * sizeof(fint));
    boxsize = (double *)malloc((nlevels + 1) * sizeof(double));
    tcenters = (double *)malloc(2 * nboxes * sizeof(double));

    pts_tree_build_(sources, ns, targ, nt, &idivflag_v, &ndiv_v, &nlmin, &nlmax,
                    &ifunif, iper, &nlevels, &nboxes, &ltree, itree, iptr,
                    tcenters, boxsize);

    isrc = (fint *)malloc(ns_v * sizeof(fint));
    isrcse = (fint *)malloc(2 * nboxes * sizeof(fint));
    itarg_arr = (fint *)malloc(nt_v * sizeof(fint));
    itargse = (fint *)malloc(2 * nboxes * sizeof(fint));
    iexpcse = (fint *)malloc(2 * nboxes * sizeof(fint));

    for (i = 1; i <= nboxes; i++) {
        iexpcse[FA2(1, i, 2)] = 1;
        iexpcse[FA2(2, i, 2)] = 0;
    }

    pts_tree_sort_(ns, sources, itree, &ltree, &nboxes, &nlevels, iptr,
                   tcenters, isrc, isrcse);
    pts_tree_sort_(nt, targ, itree, &ltree, &nboxes, &nlevels, iptr,
                   tcenters, itarg_arr, itargse);

    sourcesort = (double *)malloc(2 * ns_v * sizeof(double));
    targsort = (double *)malloc(2 * nt_v * sizeof(double));

    if (ifcharge_v == 1 && ifdipole_v == 0) {
        chargesort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        dipstrsort = (fcomplex *)malloc(nd_v * sizeof(fcomplex));
        dipvecsort = (double *)malloc(nd_v * 2 * sizeof(double));
    }
    if (ifcharge_v == 0 && ifdipole_v == 1) {
        chargesort = (fcomplex *)malloc(nd_v * sizeof(fcomplex));
        dipstrsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        dipvecsort = (double *)malloc(nd_v * 2 * ns_v * sizeof(double));
    }
    if (ifcharge_v == 1 && ifdipole_v == 1) {
        chargesort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        dipstrsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        dipvecsort = (double *)malloc(nd_v * 2 * ns_v * sizeof(double));
    }
    if (!chargesort) chargesort = (fcomplex *)malloc(nd_v * sizeof(fcomplex));
    if (!dipstrsort) dipstrsort = (fcomplex *)malloc(nd_v * sizeof(fcomplex));
    if (!dipvecsort) dipvecsort = (double *)malloc(nd_v * 2 * sizeof(double));

    if (ifpgh_v == 1) {
        potsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        gradsort = (fcomplex *)malloc(nd_v * 2 * sizeof(fcomplex));
        hesssort = (fcomplex *)malloc(nd_v * 3 * sizeof(fcomplex));
    } else if (ifpgh_v == 2) {
        potsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        gradsort = (fcomplex *)malloc(nd_v * 2 * ns_v * sizeof(fcomplex));
        hesssort = (fcomplex *)malloc(nd_v * 3 * sizeof(fcomplex));
    } else if (ifpgh_v == 3) {
        potsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        gradsort = (fcomplex *)malloc(nd_v * 2 * ns_v * sizeof(fcomplex));
        hesssort = (fcomplex *)malloc(nd_v * 3 * ns_v * sizeof(fcomplex));
    } else {
        potsort = (fcomplex *)malloc(nd_v * sizeof(fcomplex));
        gradsort = (fcomplex *)malloc(nd_v * 2 * sizeof(fcomplex));
        hesssort = (fcomplex *)malloc(nd_v * 3 * sizeof(fcomplex));
    }

    if (ifpghtarg_v == 1) {
        pottargsort = (fcomplex *)malloc(nd_v * nt_v * sizeof(fcomplex));
        gradtargsort = (fcomplex *)malloc(nd_v * 2 * sizeof(fcomplex));
        hesstargsort = (fcomplex *)malloc(nd_v * 3 * sizeof(fcomplex));
    } else if (ifpghtarg_v == 2) {
        pottargsort = (fcomplex *)malloc(nd_v * nt_v * sizeof(fcomplex));
        gradtargsort = (fcomplex *)malloc(nd_v * 2 * nt_v * sizeof(fcomplex));
        hesstargsort = (fcomplex *)malloc(nd_v * 3 * sizeof(fcomplex));
    } else if (ifpghtarg_v == 3) {
        pottargsort = (fcomplex *)malloc(nd_v * nt_v * sizeof(fcomplex));
        gradtargsort = (fcomplex *)malloc(nd_v * 2 * nt_v * sizeof(fcomplex));
        hesstargsort = (fcomplex *)malloc(nd_v * 3 * nt_v * sizeof(fcomplex));
    } else {
        pottargsort = (fcomplex *)malloc(nd_v * sizeof(fcomplex));
        gradtargsort = (fcomplex *)malloc(nd_v * 2 * sizeof(fcomplex));
        hesstargsort = (fcomplex *)malloc(nd_v * 3 * sizeof(fcomplex));
    }

    /* initialize output arrays */
    if (ifpgh_v >= 1)
        for (i = 1; i <= ns_v; i++)
            for (idim = 1; idim <= nd_v; idim++)
                potsort[FA2(idim, i, nd_v)] = 0.0;
    if (ifpgh_v >= 2)
        for (i = 1; i <= ns_v; i++)
            for (idim = 1; idim <= nd_v; idim++) {
                gradsort[FA3(idim, 1, i, nd_v, 2)] = 0.0;
                gradsort[FA3(idim, 2, i, nd_v, 2)] = 0.0;
            }
    if (ifpgh_v >= 3)
        for (i = 1; i <= ns_v; i++)
            for (idim = 1; idim <= nd_v; idim++) {
                hesssort[FA3(idim, 1, i, nd_v, 3)] = 0.0;
                hesssort[FA3(idim, 2, i, nd_v, 3)] = 0.0;
                hesssort[FA3(idim, 3, i, nd_v, 3)] = 0.0;
            }

    if (ifpghtarg_v == 1)
        for (i = 1; i <= nt_v; i++)
            for (idim = 1; idim <= nd_v; idim++) {
                pottarg[FA2(idim, i, nd_v)] = 0.0;
                pottargsort[FA2(idim, i, nd_v)] = 0.0;
            }
    if (ifpghtarg_v >= 2)
        for (i = 1; i <= nt_v; i++)
            for (idim = 1; idim <= nd_v; idim++) {
                pottargsort[FA2(idim, i, nd_v)] = 0.0;
                gradtargsort[FA3(idim, 1, i, nd_v, 2)] = 0.0;
                gradtargsort[FA3(idim, 2, i, nd_v, 2)] = 0.0;
            }
    if (ifpghtarg_v >= 3)
        for (i = 1; i <= nt_v; i++)
            for (idim = 1; idim <= nd_v; idim++) {
                hesstargsort[FA3(idim, 1, i, nd_v, 3)] = 0.0;
                hesstargsort[FA3(idim, 2, i, nd_v, 3)] = 0.0;
                hesstargsort[FA3(idim, 3, i, nd_v, 3)] = 0.0;
            }

    rscales = (double *)malloc((nlevels + 1) * sizeof(double));
    nterms = (fint *)malloc((nlevels + 1) * sizeof(fint));

    nmax = 0;
    *ier = 0;
    for (i = 0; i <= nlevels; i++) {
        fcomplex ztmp = (*zk) * boxsize[i] / (2.0 * pi);
        rscales[i] = cabs(ztmp);
        if (rscales[i] > 1.0) rscales[i] = 1.0;
        h2dterms_(&boxsize[i], zk, eps, &nterms[i], ier);
        if (nterms[i] > nmax) nmax = nterms[i];
    }

    iaddr = (fint *)malloc(4 * nboxes * sizeof(fint));
    lmptmp = (2 * nmax + 1) * nd_v;
    mptemp = (fcomplex *)malloc(lmptmp * sizeof(fcomplex));

    {
        fint two = 2;
        fint twond = 2 * nd_v;
        dreorderf_(&two, ns, sources, sourcesort, isrc);
        if (ifcharge_v == 1)
            dreorderf_(&twond, ns, (const double *)charge,
                       (double *)chargesort, isrc);
        if (ifdipole_v == 1) {
            dreorderf_(&twond, ns, (const double *)dipstr,
                       (double *)dipstrsort, isrc);
            dreorderf_(&twond, ns, dipvec, dipvecsort, isrc);
        }
        dreorderf_(&two, nt, targ, targsort, itarg_arr);
    }

    {
        const fint *laddr_p = &itree[iptr[0] - 1];
        h2dmpalloc_(&nd_v, laddr_p, iaddr, &nlevels, &lmptot, nterms);
    }

    *ier = 0;
    rmlexp = (double *)malloc(lmptot * sizeof(double));
    if (!rmlexp) { *ier = 4; goto cleanup; }

    {
        const fint *laddr_p = &itree[iptr[0] - 1];
        hfmm2dmain_(&nd_v, eps,
            zk, ns, sourcesort,
            &ifcharge_v, chargesort,
            &ifdipole_v, dipstrsort, dipvecsort,
            nt, targsort, &nexpc, expc,
            iaddr, rmlexp, mptemp, &lmptmp,
            itree, &ltree, iptr, &ndiv_v, &nlevels,
            &nboxes, iper, boxsize, rscales, tcenters,
            laddr_p,
            isrcse, itargse, iexpcse, nterms, &ntj,
            &ifpgh_v, potsort, gradsort, hesssort,
            &ifpghtarg_v, pottargsort, gradtargsort,
            hesstargsort, jexps, &scj, &ifnear_v, timeinfo, ier);
    }

    if (*ier != 0) goto cleanup;

    {
        fint twond = 2 * nd_v, fournd = 4 * nd_v, sixnd = 6 * nd_v;
        if (ifpgh_v >= 1)
            dreorderi_(&twond, ns, (const double *)potsort, (double *)pot, isrc);
        if (ifpgh_v >= 2)
            dreorderi_(&fournd, ns, (const double *)gradsort, (double *)grad, isrc);
        if (ifpgh_v >= 3)
            dreorderi_(&sixnd, ns, (const double *)hesssort, (double *)hess, isrc);

        if (ifpghtarg_v >= 1)
            dreorderi_(&twond, nt, (const double *)pottargsort, (double *)pottarg, itarg_arr);
        if (ifpghtarg_v >= 2)
            dreorderi_(&fournd, nt, (const double *)gradtargsort, (double *)gradtarg, itarg_arr);
        if (ifpghtarg_v >= 3)
            dreorderi_(&sixnd, nt, (const double *)hesstargsort, (double *)hesstarg, itarg_arr);
    }

cleanup:
    free(itree); free(boxsize); free(tcenters);
    free(isrc); free(isrcse);
    free(itarg_arr); free(itargse); free(iexpcse);
    free(sourcesort); free(targsort);
    free(chargesort); free(dipstrsort); free(dipvecsort);
    free(potsort); free(gradsort); free(hesssort);
    free(pottargsort); free(gradtargsort); free(hesstargsort);
    free(rscales); free(nterms);
    free(iaddr); free(rmlexp); free(mptemp);
}
