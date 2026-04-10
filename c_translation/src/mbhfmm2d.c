/*
 * mbhfmm2d.c - C translation of src/modified-biharmonic/mbhfmm2d.f
 *
 * Top-level 2D modified biharmonic FMM driver and helpers.
 * Translated 1:1 from the Fortran reference: same control flow,
 * same allocations, same operation order, same loop bounds.
 * OpenMP directives are stripped (sequential C).
 * Print/timing calls (prinf, prin2, cpu_time, second) are omitted.
 */

#include "mbhfmm2d.h"

/* ==================== cross-file extern declarations ==================== */

/* pts_tree2d */
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

/* fmmcommon2d */
extern void dreorderf_(const fint *ndim, const fint *n, const double *arr,
    double *arrsort, const fint *iarr);
extern void dreorderi_(const fint *ndim, const fint *n, const double *arr,
    double *arrsort, const fint *iarr);

/* l2dterms */
extern void l2dterms_(const double *eps, fint *nterms, fint *ier);

/* tree_routs2d */
extern void computemnlists_(const fint *nlevels, const fint *nboxes,
    const fint *itree, const fint *ltree,
    const fint *iptr, const double *centers,
    const double *boxsize, const fint *iper,
    fint *mnlist1, fint *mnlist2,
    fint *mnlist3, fint *mnlist4);

extern void computelists_(const fint *nlevels, const fint *nboxes,
    const fint *itree, const fint *ltree,
    const fint *iptr, const double *centers,
    const double *boxsize, const fint *iper,
    const fint *mnlist1, fint *nlist1, fint *list1,
    const fint *mnlist2, fint *nlist2, fint *list2,
    const fint *mnlist3, fint *nlist3, fint *list3,
    const fint *mnlist4, fint *nlist4, fint *list4);

/* mbhrouts2d */
extern void mbh2d_init_carray_(double *carray, const fint *ldc);
extern void mbh2dconvtomp_vec_(const fint *nd, const double *beta,
    const fint *ns, const fint *ifcharge, const double *charge,
    const fint *ifdipole, const double *dipstr, const double *dipvec,
    const fint *ifquad, const double *quadstr, const double *quadvec,
    const fint *ifoct, const double *octstr, const double *octvec,
    const fint *nterms, fcomplex *mbhmpole, fcomplex *ympole);
extern void mbh2dformmpmp_vec_(const fint *nd, const double *beta,
    const double *rscale, const double *source, const fint *ns,
    fcomplex *mbhmpolesrc, fcomplex *ympolesrc,
    const fint *ntermsrc, const double *center,
    const fint *nterms, fcomplex *mbhmpole, fcomplex *ympole);
extern void mbh2dformtamp_vec_(const fint *nd, const double *beta,
    const double *rscale, const double *source, const fint *ns,
    fcomplex *mbhmpolesrc, fcomplex *ympolesrc,
    const fint *ntermsrc, const double *center,
    const fint *nterms, fcomplex *mbhloc, fcomplex *lloc);
extern void mbh2dmpmp_vec_(const fint *nd, const double *beta,
    const double *rscale1, const double *center1,
    const fcomplex *mbhmpole1, const fcomplex *ympole1, const fint *nterms1,
    const double *rscale2, const double *center2,
    fcomplex *mbhmpole2, fcomplex *ympole2, const fint *nterms2);
extern void mbh2dmploc_vec_(const fint *nd, const double *beta,
    const double *rscale1, const double *center1,
    const fcomplex *mbhmpole, const fcomplex *ympole, const fint *nterms1,
    const double *rscale2, const double *center2,
    fcomplex *mbhloc, fcomplex *lloc, const fint *nterms2);
extern void mbh2dlocloc_vec_(const fint *nd, const double *beta,
    const double *rscale1, const double *center1,
    const fcomplex *mbhloc1, const fcomplex *lloc1, const fint *nterms1,
    const double *rscale2, const double *center2,
    fcomplex *mbhloc2, fcomplex *lloc2, const fint *nterms2,
    const double *carray, const fint *ldc);
extern void mbh2dmpevalp_vec_(const fint *nd, const double *beta,
    const double *rscale, const double *center,
    const fcomplex *mbhmpole, const fcomplex *ympole,
    const fint *nterms, const double *ztarg, const fint *ntarg,
    double *pot);
extern void mbh2dmpevalg_vec_(const fint *nd, const double *beta,
    const double *rscale, const double *center,
    const fcomplex *mbhmpole, const fcomplex *ympole,
    const fint *nterms, const double *ztarg, const fint *ntarg,
    double *pot, double *grad);
extern void mbh2dmpevalh_vec_(const fint *nd, const double *beta,
    const double *rscale, const double *center,
    const fcomplex *mbhmpole, const fcomplex *ympole,
    const fint *nterms, const double *ztarg, const fint *ntarg,
    double *pot, double *grad, double *hess);
extern void mbh2dtaevalp_vec_(const fint *nd, const double *beta,
    const double *rscale, const double *center,
    const fcomplex *mbhloc, const fcomplex *lloc,
    const fint *nterms, const double *ztarg, const fint *ntarg,
    double *pot);
extern void mbh2dtaevalg_vec_(const fint *nd, const double *beta,
    const double *rscale, const double *center,
    const fcomplex *mbhloc, const fcomplex *lloc,
    const fint *nterms, const double *ztarg, const fint *ntarg,
    double *pot, double *grad);
extern void mbh2dtaevalh_vec_(const fint *nd, const double *beta,
    const double *rscale, const double *center,
    const fcomplex *mbhloc, const fcomplex *lloc,
    const fint *nterms, const double *ztarg, const fint *ntarg,
    double *pot, double *grad, double *hess);

/* mbhkernels2d */
extern void mbh2d_directmpsp_vec_(const fint *nd, const double *beta,
    const double *sources, const fint *ns,
    const fcomplex *mbhmps, const fcomplex *ymps, const fint *ntermsmps,
    const double *targ, const fint *nt,
    double *pot, const double *thresh);
extern void mbh2d_directmpsg_vec_(const fint *nd, const double *beta,
    const double *sources, const fint *ns,
    const fcomplex *mbhmps, const fcomplex *ymps, const fint *ntermsmps,
    const double *targ, const fint *nt,
    double *pot, double *grad, const double *thresh);
extern void mbh2d_directmpsh_vec_(const fint *nd, const double *beta,
    const double *sources, const fint *ns,
    const fcomplex *mbhmps, const fcomplex *ymps, const fint *ntermsmps,
    const double *targ, const fint *nt,
    double *pot, double *grad, double *hess, const double *thresh);


/* laddr(k, ilev): k is 1-based, ilev is 0-based, leading dim 2 */
#define LADDR(k, ilev) ((ilev) * 2 + ((k) - 1))


/* ================================================================
 * mbh2dmpzero - zero out one multipole/local expansion
 * ================================================================ */
void FNAME(mbh2dmpzero)(const fint *nd, fcomplex *mpole, const fint *nterms)
{
    fint nd_v = *nd, nt1 = *nterms + 1;
    for (fint n = 0; n < nt1; n++) {
        for (fint idim = 0; idim < nd_v; idim++) {
            mpole[n * nd_v + idim] = 0.0;
        }
    }
}


/* ================================================================
 * mbh2dmpalloc - compute workspace layout for expansions
 * iaddr(4, nboxes): 4 expansion pointers per box
 * ================================================================ */
void FNAME(mbh2dmpalloc)(const fint *nd_p, const fint *laddr, fint *iaddr,
    const fint *nlevels_p, fint *lmptot, const fint *nterms)
{
    fint nd = *nd_p;
    fint nlevels = *nlevels_p;

    /* nn = (nterms+1)*2*nd doubles per expansion set.
     * Each box gets 2 sets for multipole (mbh + y). */
    fint istart = 1;
    for (fint i = 0; i <= nlevels; i++) {
        fint nn = (nterms[i] + 1) * 2 * nd;
        for (fint ibox = laddr[LADDR(1, i)]; ibox <= laddr[LADDR(2, i)]; ibox++) {
            fint itmp = ibox - laddr[LADDR(1, i)];
            iaddr[FA2(1, ibox, 4)] = istart + itmp * nn * 2;
            iaddr[FA2(2, ibox, 4)] = istart + itmp * nn * 2 + nn;
        }
        istart = istart + (laddr[LADDR(2, i)] - laddr[LADDR(1, i)] + 1) * nn * 2;
    }

    /* Each box gets 2 sets for local (mbh + laplace). */
    for (fint i = 0; i <= nlevels; i++) {
        fint nn = (nterms[i] + 1) * 2 * nd;
        for (fint ibox = laddr[LADDR(1, i)]; ibox <= laddr[LADDR(2, i)]; ibox++) {
            fint itmp = ibox - laddr[LADDR(1, i)];
            iaddr[FA2(3, ibox, 4)] = istart + itmp * nn * 2;
            iaddr[FA2(4, ibox, 4)] = istart + itmp * nn * 2 + nn;
        }
        istart = istart + (laddr[LADDR(2, i)] - laddr[LADDR(1, i)] + 1) * nn * 2;
    }

    *lmptot = istart;
}


/* ================================================================
 * mbhfmm2dmps_direct_vec - direct multipolar-source evaluation
 * ================================================================ */
void FNAME(mbhfmm2dmps_direct_vec)(const fint *nd_p,
    const fint *istart_p, const fint *iend_p,
    const fint *jstart_p, const fint *jend_p,
    const double *beta,
    const double *source, const fcomplex *mbhmps, const fcomplex *ymps,
    const fint *ntermsmps_p,
    const double *targ, const fint *ifpgh_p,
    double *pot, double *grad, double *hess,
    const double *thresh)
{
    fint nd = *nd_p;
    fint istart = *istart_p, iend = *iend_p;
    fint jstart = *jstart_p, jend = *jend_p;
    fint ifpgh = *ifpgh_p;
    fint ntermsmps = *ntermsmps_p;

    fint ns = iend - istart + 1;
    fint nt = jend - jstart + 1;

    fint mpstride = nd * (ntermsmps + 1);

    if (ifpgh == 1) {
        mbh2d_directmpsp_vec_(&nd, beta,
            source + 2 * (istart - 1), &ns,
            mbhmps + mpstride * (istart - 1),
            ymps + mpstride * (istart - 1), &ntermsmps,
            targ + 2 * (jstart - 1), &nt,
            pot + nd * (jstart - 1), thresh);
    }
    if (ifpgh == 2) {
        mbh2d_directmpsg_vec_(&nd, beta,
            source + 2 * (istart - 1), &ns,
            mbhmps + mpstride * (istart - 1),
            ymps + mpstride * (istart - 1), &ntermsmps,
            targ + 2 * (jstart - 1), &nt,
            pot + nd * (jstart - 1),
            grad + 2 * nd * (jstart - 1), thresh);
    }
    if (ifpgh == 3) {
        mbh2d_directmpsh_vec_(&nd, beta,
            source + 2 * (istart - 1), &ns,
            mbhmps + mpstride * (istart - 1),
            ymps + mpstride * (istart - 1), &ntermsmps,
            targ + 2 * (jstart - 1), &nt,
            pot + nd * (jstart - 1),
            grad + 2 * nd * (jstart - 1),
            hess + 3 * nd * (jstart - 1), thresh);
    }
}


/* ================================================================
 * mbhfmm2dmain - main FMM driver (8 steps)
 * ================================================================ */
void FNAME(mbhfmm2dmain)(const fint *nd_p, const double *eps,
    const double *beta_p, const fint *nsource, const double *source,
    const fint *ntermsmps_p, const fcomplex *mbhmps, const fcomplex *ymps,
    const fint *ntarget, const double *target,
    const fint *iaddr, double *rmlexp,
    const double *carray, const fint *ldc_p,
    const fint *itree, const fint *ltree,
    const fint *iptr, const fint *ndiv, const fint *nlevels_p,
    const fint *nboxes_p, const fint *iper,
    const double *boxsize, const double *rscales,
    const double *centers, const fint *laddr,
    const fint *isrcse, const fint *itargse,
    const fint *nterms,
    const fint *ifpgh_p, double *pot, double *grad, double *hess,
    const fint *ifpghtarg_p, double *pottarg, double *gradtarg,
    double *hesstarg, fint *ier)
{
    fint nd = *nd_p;
    double beta = *beta_p;
    fint ntermsmps = *ntermsmps_p;
    fint nlevels = *nlevels_p;
    fint nboxes = *nboxes_p;
    fint ldc = *ldc_p;
    fint ifpgh = *ifpgh_p;
    fint ifpghtarg = *ifpghtarg_p;

    (void)eps; (void)nsource; (void)ntarget; (void)ltree; (void)ndiv;

    double zkiupbound = 40.0;
    double zi = beta;

    /* compute list info */
    fint mnlist1 = 0, mnlist2 = 0, mnlist3 = 0, mnlist4 = 0;
    fint iper_v = *iper;
    computemnlists_(&nlevels, &nboxes, itree, ltree, iptr, centers,
        boxsize, &iper_v, &mnlist1, &mnlist2, &mnlist3, &mnlist4);

    fint *nlist1s = (fint *)malloc(nboxes * sizeof(fint));
    fint *list1 = (fint *)malloc(mnlist1 * nboxes * sizeof(fint));
    fint *nlist2s = (fint *)malloc(nboxes * sizeof(fint));
    fint *list2 = (fint *)malloc(mnlist2 * nboxes * sizeof(fint));
    fint *nlist3s = (fint *)malloc(nboxes * sizeof(fint));
    fint *list3 = (fint *)malloc(mnlist3 * nboxes * sizeof(fint));
    fint *nlist4s = (fint *)malloc(nboxes * sizeof(fint));
    fint *list4 = (fint *)malloc(mnlist4 * nboxes * sizeof(fint));

    computelists_(&nlevels, &nboxes, itree, ltree, iptr, centers,
        boxsize, &iper_v,
        &mnlist1, nlist1s, list1,
        &mnlist2, nlist2s, list2,
        &mnlist3, nlist3s, list3,
        &mnlist4, nlist4s, list4);

    /* zero all multipole and local expansions */
    for (fint ilev = 0; ilev <= nlevels; ilev++) {
        for (fint ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            fint nt_ilev = nterms[ilev];
            FNAME(mbh2dmpzero)(&nd, (fcomplex *)(rmlexp + iaddr[FA2(1,ibox,4)] - 1), &nt_ilev);
            FNAME(mbh2dmpzero)(&nd, (fcomplex *)(rmlexp + iaddr[FA2(2,ibox,4)] - 1), &nt_ilev);
            FNAME(mbh2dmpzero)(&nd, (fcomplex *)(rmlexp + iaddr[FA2(3,ibox,4)] - 1), &nt_ilev);
            FNAME(mbh2dmpzero)(&nd, (fcomplex *)(rmlexp + iaddr[FA2(4,ibox,4)] - 1), &nt_ilev);
        }
    }

    /* STEP 1: form multipole expansions at leaf boxes */
    for (fint ilev = 2; ilev <= nlevels; ilev++) {
        if (zi * boxsize[ilev] < zkiupbound) {
            for (fint ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
                fint nchild = itree[iptr[3] + ibox - 1 - 1]; /* iptr(4) is 1-based */
                fint istart = isrcse[FA2(1, ibox, 2)];
                fint iend = isrcse[FA2(2, ibox, 2)];
                fint npts = iend - istart + 1;
                if (nchild == 0 && npts > 0) {
                    fint mpstride = nd * (ntermsmps + 1);
                    fint nt_ilev = nterms[ilev];
                    mbh2dformmpmp_vec_(&nd, &beta, &rscales[ilev],
                        source + 2 * (istart - 1), &npts,
                        (fcomplex *)(mbhmps + mpstride * (istart - 1)),
                        (fcomplex *)(ymps + mpstride * (istart - 1)),
                        &ntermsmps, centers + 2 * (ibox - 1),
                        &nt_ilev,
                        (fcomplex *)(rmlexp + iaddr[FA2(1,ibox,4)] - 1),
                        (fcomplex *)(rmlexp + iaddr[FA2(2,ibox,4)] - 1));
                }
            }
        }
    }

    /* STEP 2: form local expansions from list-4 (far field) */
    for (fint ilev = 2; ilev <= nlevels; ilev++) {
        if (zi * boxsize[ilev] < zkiupbound) {
            for (fint ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
                fint npts = 0;
                if (ifpghtarg > 0) {
                    fint istart = itargse[FA2(1, ibox, 2)];
                    fint iend = itargse[FA2(2, ibox, 2)];
                    npts += iend - istart + 1;
                }
                if (ifpgh > 0) {
                    fint istart = isrcse[FA2(1, ibox, 2)];
                    fint iend = isrcse[FA2(2, ibox, 2)];
                    npts += iend - istart + 1;
                }
                if (npts > 0) {
                    for (fint i = 1; i <= nlist4s[ibox - 1]; i++) {
                        fint jbox = list4[FA2(i, ibox, mnlist4)];
                        fint jstart = isrcse[FA2(1, jbox, 2)];
                        fint jend = isrcse[FA2(2, jbox, 2)];
                        fint jnpts = jend - jstart + 1;
                        fint mpstride = nd * (ntermsmps + 1);
                        fint nt_ilev = nterms[ilev];
                        mbh2dformtamp_vec_(&nd, &beta, &rscales[ilev],
                            source + 2 * (jstart - 1), &jnpts,
                            (fcomplex *)(mbhmps + mpstride * (jstart - 1)),
                            (fcomplex *)(ymps + mpstride * (jstart - 1)),
                            &ntermsmps, centers + 2 * (ibox - 1),
                            &nt_ilev,
                            (fcomplex *)(rmlexp + iaddr[FA2(3,ibox,4)] - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(4,ibox,4)] - 1));
                    }
                }
            }
        }
    }

    /* STEP 3: merge multipoles (upward pass) */
    for (fint ilev = nlevels - 1; ilev >= 1; ilev--) {
        if (zi * boxsize[ilev] < zkiupbound) {
            for (fint ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
                fint nchild = itree[iptr[3] + ibox - 1 - 1];
                for (fint i = 1; i <= nchild; i++) {
                    fint jbox = itree[iptr[4] + 4 * (ibox - 1) + i - 1 - 1];
                    fint jstart = isrcse[FA2(1, jbox, 2)];
                    fint jend = isrcse[FA2(2, jbox, 2)];
                    fint jnpts = jend - jstart + 1;
                    if (jnpts > 0) {
                        fint nt_child = nterms[ilev + 1];
                        fint nt_ilev = nterms[ilev];
                        mbh2dmpmp_vec_(&nd, &beta,
                            &rscales[ilev + 1], centers + 2 * (jbox - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(1,jbox,4)] - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(2,jbox,4)] - 1),
                            &nt_child,
                            &rscales[ilev], centers + 2 * (ibox - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(1,ibox,4)] - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(2,ibox,4)] - 1),
                            &nt_ilev);
                    }
                }
            }
        }
    }

    /* STEP 4: multipole-to-local (list-2) */
    for (fint ilev = 2; ilev <= nlevels; ilev++) {
        if (zi * boxsize[ilev] < zkiupbound) {
            for (fint ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
                fint npts = 0;
                if (ifpghtarg > 0) {
                    fint istart = itargse[FA2(1, ibox, 2)];
                    fint iend = itargse[FA2(2, ibox, 2)];
                    npts += iend - istart + 1;
                }
                if (ifpgh > 0) {
                    fint istart = isrcse[FA2(1, ibox, 2)];
                    fint iend = isrcse[FA2(2, ibox, 2)];
                    npts += iend - istart + 1;
                }
                if (npts > 0) {
                    for (fint i = 1; i <= nlist2s[ibox - 1]; i++) {
                        fint jbox = list2[FA2(i, ibox, mnlist2)];
                        fint nt_ilev = nterms[ilev];
                        mbh2dmploc_vec_(&nd, &beta,
                            &rscales[ilev], centers + 2 * (jbox - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(1,jbox,4)] - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(2,jbox,4)] - 1),
                            &nt_ilev,
                            &rscales[ilev], centers + 2 * (ibox - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(3,ibox,4)] - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(4,ibox,4)] - 1),
                            &nt_ilev);
                    }
                }
            }
        }
    }

    /* STEP 5: split local (downward pass) */
    for (fint ilev = 1; ilev <= nlevels - 1; ilev++) {
        if (zi * boxsize[ilev] < zkiupbound) {
            for (fint ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
                fint nchild = itree[iptr[3] + ibox - 1 - 1];
                fint npts = 0;
                if (ifpghtarg > 0) {
                    fint istart = itargse[FA2(1, ibox, 2)];
                    fint iend = itargse[FA2(2, ibox, 2)];
                    npts += iend - istart + 1;
                }
                if (ifpgh > 0) {
                    fint istart = isrcse[FA2(1, ibox, 2)];
                    fint iend = isrcse[FA2(2, ibox, 2)];
                    npts += iend - istart + 1;
                }
                if (npts > 0) {
                    for (fint i = 1; i <= nchild; i++) {
                        fint jbox = itree[iptr[4] + 4 * (ibox - 1) + i - 1 - 1];
                        fint nt_ilev = nterms[ilev];
                        fint nt_child = nterms[ilev + 1];
                        mbh2dlocloc_vec_(&nd, &beta,
                            &rscales[ilev], centers + 2 * (ibox - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(3,ibox,4)] - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(4,ibox,4)] - 1),
                            &nt_ilev,
                            &rscales[ilev + 1], centers + 2 * (jbox - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(3,jbox,4)] - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(4,jbox,4)] - 1),
                            &nt_child, carray, &ldc);
                    }
                }
            }
        }
    }

    /* STEP 6: evaluate multipoles at list-3 targets/sources */
    for (fint ilev = 1; ilev <= nlevels - 1; ilev++) {
        if (zi * boxsize[ilev + 1] < zkiupbound) {
            for (fint ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
                /* evaluate at targets */
                fint istart = itargse[FA2(1, ibox, 2)];
                fint iend = itargse[FA2(2, ibox, 2)];
                fint npts = iend - istart + 1;

                if (ifpghtarg == 1) {
                    for (fint i = 1; i <= nlist3s[ibox - 1]; i++) {
                        fint jbox = list3[FA2(i, ibox, mnlist3)];
                        fint nt_child = nterms[ilev + 1];
                        mbh2dmpevalp_vec_(&nd, &beta, &rscales[ilev + 1],
                            centers + 2 * (jbox - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(1,jbox,4)] - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(2,jbox,4)] - 1),
                            &nt_child, target + 2 * (istart - 1), &npts,
                            pottarg + nd * (istart - 1));
                    }
                }
                if (ifpghtarg == 2) {
                    for (fint i = 1; i <= nlist3s[ibox - 1]; i++) {
                        fint jbox = list3[FA2(i, ibox, mnlist3)];
                        fint nt_child = nterms[ilev + 1];
                        mbh2dmpevalg_vec_(&nd, &beta, &rscales[ilev + 1],
                            centers + 2 * (jbox - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(1,jbox,4)] - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(2,jbox,4)] - 1),
                            &nt_child, target + 2 * (istart - 1), &npts,
                            pottarg + nd * (istart - 1),
                            gradtarg + 2 * nd * (istart - 1));
                    }
                }
                if (ifpghtarg == 3) {
                    for (fint i = 1; i <= nlist3s[ibox - 1]; i++) {
                        fint jbox = list3[FA2(i, ibox, mnlist3)];
                        fint nt_child = nterms[ilev + 1];
                        mbh2dmpevalh_vec_(&nd, &beta, &rscales[ilev + 1],
                            centers + 2 * (jbox - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(1,jbox,4)] - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(2,jbox,4)] - 1),
                            &nt_child, target + 2 * (istart - 1), &npts,
                            pottarg + nd * (istart - 1),
                            gradtarg + 2 * nd * (istart - 1),
                            hesstarg + 3 * nd * (istart - 1));
                    }
                }

                /* evaluate at sources */
                istart = isrcse[FA2(1, ibox, 2)];
                iend = isrcse[FA2(2, ibox, 2)];
                npts = iend - istart + 1;

                if (ifpgh == 1) {
                    for (fint i = 1; i <= nlist3s[ibox - 1]; i++) {
                        fint jbox = list3[FA2(i, ibox, mnlist3)];
                        fint nt_child = nterms[ilev + 1];
                        mbh2dmpevalp_vec_(&nd, &beta, &rscales[ilev + 1],
                            centers + 2 * (jbox - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(1,jbox,4)] - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(2,jbox,4)] - 1),
                            &nt_child, source + 2 * (istart - 1), &npts,
                            pot + nd * (istart - 1));
                    }
                }
                if (ifpgh == 2) {
                    for (fint i = 1; i <= nlist3s[ibox - 1]; i++) {
                        fint jbox = list3[FA2(i, ibox, mnlist3)];
                        fint nt_child = nterms[ilev + 1];
                        mbh2dmpevalg_vec_(&nd, &beta, &rscales[ilev + 1],
                            centers + 2 * (jbox - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(1,jbox,4)] - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(2,jbox,4)] - 1),
                            &nt_child, source + 2 * (istart - 1), &npts,
                            pot + nd * (istart - 1),
                            grad + 2 * nd * (istart - 1));
                    }
                }
                if (ifpgh == 3) {
                    for (fint i = 1; i <= nlist3s[ibox - 1]; i++) {
                        fint jbox = list3[FA2(i, ibox, mnlist3)];
                        fint nt_child = nterms[ilev + 1];
                        mbh2dmpevalh_vec_(&nd, &beta, &rscales[ilev + 1],
                            centers + 2 * (jbox - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(1,jbox,4)] - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(2,jbox,4)] - 1),
                            &nt_child, source + 2 * (istart - 1), &npts,
                            pot + nd * (istart - 1),
                            grad + 2 * nd * (istart - 1),
                            hess + 3 * nd * (istart - 1));
                    }
                }
            }
        }
    }

    /* STEP 7: evaluate local expansions at leaf boxes */
    for (fint ilev = 0; ilev <= nlevels; ilev++) {
        if (zi * boxsize[ilev] < zkiupbound) {
            for (fint ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
                fint nchild = itree[iptr[3] + ibox - 1 - 1];
                if (nchild == 0) {
                    fint nt_ilev = nterms[ilev];

                    /* evaluate at targets */
                    fint istart = itargse[FA2(1, ibox, 2)];
                    fint iend = itargse[FA2(2, ibox, 2)];
                    fint npts = iend - istart + 1;
                    if (ifpghtarg == 1) {
                        mbh2dtaevalp_vec_(&nd, &beta, &rscales[ilev],
                            centers + 2 * (ibox - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(3,ibox,4)] - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(4,ibox,4)] - 1),
                            &nt_ilev, target + 2 * (istart - 1), &npts,
                            pottarg + nd * (istart - 1));
                    }
                    if (ifpghtarg == 2) {
                        mbh2dtaevalg_vec_(&nd, &beta, &rscales[ilev],
                            centers + 2 * (ibox - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(3,ibox,4)] - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(4,ibox,4)] - 1),
                            &nt_ilev, target + 2 * (istart - 1), &npts,
                            pottarg + nd * (istart - 1),
                            gradtarg + 2 * nd * (istart - 1));
                    }
                    if (ifpghtarg == 3) {
                        mbh2dtaevalh_vec_(&nd, &beta, &rscales[ilev],
                            centers + 2 * (ibox - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(3,ibox,4)] - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(4,ibox,4)] - 1),
                            &nt_ilev, target + 2 * (istart - 1), &npts,
                            pottarg + nd * (istart - 1),
                            gradtarg + 2 * nd * (istart - 1),
                            hesstarg + 3 * nd * (istart - 1));
                    }

                    /* evaluate at sources */
                    istart = isrcse[FA2(1, ibox, 2)];
                    iend = isrcse[FA2(2, ibox, 2)];
                    npts = iend - istart + 1;
                    if (ifpgh == 1) {
                        mbh2dtaevalp_vec_(&nd, &beta, &rscales[ilev],
                            centers + 2 * (ibox - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(3,ibox,4)] - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(4,ibox,4)] - 1),
                            &nt_ilev, source + 2 * (istart - 1), &npts,
                            pot + nd * (istart - 1));
                    }
                    if (ifpgh == 2) {
                        mbh2dtaevalg_vec_(&nd, &beta, &rscales[ilev],
                            centers + 2 * (ibox - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(3,ibox,4)] - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(4,ibox,4)] - 1),
                            &nt_ilev, source + 2 * (istart - 1), &npts,
                            pot + nd * (istart - 1),
                            grad + 2 * nd * (istart - 1));
                    }
                    if (ifpgh == 3) {
                        mbh2dtaevalh_vec_(&nd, &beta, &rscales[ilev],
                            centers + 2 * (ibox - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(3,ibox,4)] - 1),
                            (fcomplex *)(rmlexp + iaddr[FA2(4,ibox,4)] - 1),
                            &nt_ilev, source + 2 * (istart - 1), &npts,
                            pot + nd * (istart - 1),
                            grad + 2 * nd * (istart - 1),
                            hess + 3 * nd * (istart - 1));
                    }
                }
            }
        }
    }

    /* STEP 8: direct evaluation (list-1) */
    double thresh = boxsize[0] * 1.0e-16;
    for (fint ilev = 0; ilev <= nlevels; ilev++) {
        for (fint ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            fint istartt = itargse[FA2(1, ibox, 2)];
            fint iendt = itargse[FA2(2, ibox, 2)];
            fint istarts = isrcse[FA2(1, ibox, 2)];
            fint iends = isrcse[FA2(2, ibox, 2)];

            for (fint i = 1; i <= nlist1s[ibox - 1]; i++) {
                fint jbox = list1[FA2(i, ibox, mnlist1)];
                fint jstart = isrcse[FA2(1, jbox, 2)];
                fint jend = isrcse[FA2(2, jbox, 2)];

                FNAME(mbhfmm2dmps_direct_vec)(&nd, &jstart, &jend,
                    &istartt, &iendt, &beta, source,
                    mbhmps, ymps, &ntermsmps,
                    target, &ifpghtarg, pottarg, gradtarg, hesstarg,
                    &thresh);

                FNAME(mbhfmm2dmps_direct_vec)(&nd, &jstart, &jend,
                    &istarts, &iends, &beta, source,
                    mbhmps, ymps, &ntermsmps,
                    source, &ifpgh, pot, grad, hess,
                    &thresh);
            }
        }
    }

    free(nlist1s); free(list1);
    free(nlist2s); free(list2);
    free(nlist3s); free(list3);
    free(nlist4s); free(list4);
}


/* ================================================================
 * mbhfmm2d - top-level FMM wrapper
 * ================================================================ */
void FNAME(mbhfmm2d)(const fint *nd_p, const double *eps, const double *beta_p,
    const fint *ns_p, const double *sources,
    const fint *ifcharge, const double *charge,
    const fint *ifdipole, const double *dipstr, const double *dipvec,
    const fint *ifquadpole, const double *quadstr, const double *quadvec,
    const fint *ifoctpole, const double *octstr, const double *octvec,
    const fint *iper, const fint *ifpgh_p,
    double *pot, double *grad, double *hess,
    const fint *nt_p, const double *targ, const fint *ifpghtarg_p,
    double *pottarg, double *gradtarg, double *hesstarg,
    fint *ier)
{
    fint nd = *nd_p, ns = *ns_p, nt = *nt_p;
    double beta = *beta_p;
    fint ifpgh = *ifpgh_p, ifpghtarg = *ifpghtarg_p;
    double pi = 4.0 * atan(1.0);
    fcomplex zk = beta * I;

    fint idivflag = 0, ndiv = 20, nlmin = 0, nlmax = 51, ifunif = 0;
    fint iper_v = 0;
    fint nlevels = 0, nboxes = 0, ltree = 0;

    /* tree memory */
    pts_tree_mem_(sources, &ns, targ, &nt, &idivflag, &ndiv,
        &nlmin, &nlmax, &ifunif, &iper_v, &nlevels, &nboxes, &ltree);

    fint *itree = (fint *)malloc(ltree * sizeof(fint));
    double *boxsize = (double *)malloc((nlevels + 1) * sizeof(double));
    double *tcenters = (double *)malloc(2 * nboxes * sizeof(double));
    fint iptr[8];

    pts_tree_build_(sources, &ns, targ, &nt, &idivflag, &ndiv,
        &nlmin, &nlmax, &ifunif, &iper_v,
        &nlevels, &nboxes, &ltree, itree, iptr, tcenters, boxsize);

    fint *isrc = (fint *)malloc(ns * sizeof(fint));
    fint *isrcse = (fint *)malloc(2 * nboxes * sizeof(fint));
    fint *itarg_arr = (fint *)malloc(nt * sizeof(fint));
    fint *itargse = (fint *)malloc(2 * nboxes * sizeof(fint));

    pts_tree_sort_(&ns, sources, itree, &ltree, &nboxes, &nlevels, iptr,
        tcenters, isrc, isrcse);
    pts_tree_sort_(&nt, targ, itree, &ltree, &nboxes, &nlevels, iptr,
        tcenters, itarg_arr, itargse);

    double *sourcesort = (double *)malloc(2 * ns * sizeof(double));
    double *targsort = (double *)malloc(2 * nt * sizeof(double));

    /* compute ntermsmps */
    fint ntermsmps = 0;
    if (*ifdipole == 1) ntermsmps = 1;
    if (*ifquadpole == 1) ntermsmps = 2;
    if (*ifoctpole == 1) ntermsmps = 3;

    /* allocate and initialize multipolar sources */
    fcomplex *mbhmps = (fcomplex *)calloc(nd * (ntermsmps + 1) * ns, sizeof(fcomplex));
    fcomplex *ymps = (fcomplex *)calloc(nd * (ntermsmps + 1) * ns, sizeof(fcomplex));
    fcomplex *mbhmpssort = (fcomplex *)malloc(nd * (ntermsmps + 1) * ns * sizeof(fcomplex));
    fcomplex *ympssort = (fcomplex *)malloc(nd * (ntermsmps + 1) * ns * sizeof(fcomplex));

    mbh2dconvtomp_vec_(&nd, &beta, &ns, ifcharge, charge,
        ifdipole, dipstr, dipvec, ifquadpole, quadstr, quadvec,
        ifoctpole, octstr, octvec, &ntermsmps, mbhmps, ymps);

    /* allocate sorted output arrays */
    double *potsort = NULL, *gradsort = NULL, *hesssort = NULL;
    if (ifpgh == 1) {
        potsort = (double *)malloc(nd * ns * sizeof(double));
        gradsort = (double *)malloc(nd * 2 * sizeof(double));
        hesssort = (double *)malloc(nd * 3 * sizeof(double));
    } else if (ifpgh == 2) {
        potsort = (double *)malloc(nd * ns * sizeof(double));
        gradsort = (double *)malloc(nd * 2 * ns * sizeof(double));
        hesssort = (double *)malloc(nd * 3 * sizeof(double));
    } else if (ifpgh == 3) {
        potsort = (double *)malloc(nd * ns * sizeof(double));
        gradsort = (double *)malloc(nd * 2 * ns * sizeof(double));
        hesssort = (double *)malloc(nd * 3 * ns * sizeof(double));
    } else {
        potsort = (double *)malloc(nd * sizeof(double));
        gradsort = (double *)malloc(nd * 2 * sizeof(double));
        hesssort = (double *)malloc(nd * 3 * sizeof(double));
    }

    double *pottargsort = NULL, *gradtargsort = NULL, *hesstargsort = NULL;
    if (ifpghtarg == 1) {
        pottargsort = (double *)malloc(nd * nt * sizeof(double));
        gradtargsort = (double *)malloc(nd * 2 * sizeof(double));
        hesstargsort = (double *)malloc(nd * 3 * sizeof(double));
    } else if (ifpghtarg == 2) {
        pottargsort = (double *)malloc(nd * nt * sizeof(double));
        gradtargsort = (double *)malloc(nd * 2 * nt * sizeof(double));
        hesstargsort = (double *)malloc(nd * 3 * sizeof(double));
    } else if (ifpghtarg == 3) {
        pottargsort = (double *)malloc(nd * nt * sizeof(double));
        gradtargsort = (double *)malloc(nd * 2 * nt * sizeof(double));
        hesstargsort = (double *)malloc(nd * 3 * nt * sizeof(double));
    } else {
        pottargsort = (double *)malloc(nd * sizeof(double));
        gradtargsort = (double *)malloc(nd * 2 * sizeof(double));
        hesstargsort = (double *)malloc(nd * 3 * sizeof(double));
    }

    /* zero output arrays */
    if (ifpgh >= 1) {
        for (fint i = 0; i < nd * ns; i++) potsort[i] = 0.0;
    }
    if (ifpgh >= 2) {
        for (fint i = 0; i < nd * 2 * ns; i++) gradsort[i] = 0.0;
    }
    if (ifpgh >= 3) {
        for (fint i = 0; i < nd * 3 * ns; i++) hesssort[i] = 0.0;
    }

    if (ifpghtarg >= 1) {
        for (fint i = 0; i < nd * nt; i++) {
            pottarg[i] = 0.0;
            pottargsort[i] = 0.0;
        }
    }
    if (ifpghtarg >= 2) {
        for (fint i = 0; i < nd * 2 * nt; i++) gradtargsort[i] = 0.0;
    }
    if (ifpghtarg >= 3) {
        for (fint i = 0; i < nd * 3 * nt; i++) hesstargsort[i] = 0.0;
    }

    /* compute rscales and nterms */
    double *rscales = (double *)malloc((nlevels + 1) * sizeof(double));
    fint *nterms_arr = (fint *)malloc((nlevels + 1) * sizeof(fint));
    fint nmax = 0, nterms1;
    *ier = 0;
    for (fint i = 0; i <= nlevels; i++) {
        rscales[i] = fmin(cabs(zk * boxsize[i] / (2.0 * pi)), 1.0);
        l2dterms_(eps, &nterms1, ier);
        nterms_arr[i] = nterms1;
        if (nterms_arr[i] > nmax) nmax = nterms_arr[i];
    }

    fint ldc = nmax + 5;
    double *carray = (double *)malloc((ldc + 1) * (ldc + 1) * sizeof(double));
    mbh2d_init_carray_(carray, &ldc);

    fint *iaddr = (fint *)malloc(4 * nboxes * sizeof(fint));

    /* reorder sources */
    fint two = 2;
    dreorderf_(&two, &ns, sources, sourcesort, isrc);
    fint ndtmp = nd * 2 * (ntermsmps + 1);
    dreorderf_(&ndtmp, &ns, (const double *)mbhmps, (double *)mbhmpssort, isrc);
    dreorderf_(&ndtmp, &ns, (const double *)ymps, (double *)ympssort, isrc);

    /* reorder targets */
    dreorderf_(&two, &nt, targ, targsort, itarg_arr);

    /* allocate workspace */
    fint lmptot;
    FNAME(mbh2dmpalloc)(&nd, itree + iptr[0] - 1, iaddr, &nlevels,
        &lmptot, nterms_arr);

    *ier = 0;
    double *rmlexp = (double *)malloc(lmptot * sizeof(double));
    if (rmlexp == NULL) {
        *ier = 4;
        /* cleanup and return */
        free(iaddr); free(carray); free(nterms_arr); free(rscales);
        free(pottargsort); free(gradtargsort); free(hesstargsort);
        free(potsort); free(gradsort); free(hesssort);
        free(mbhmpssort); free(ympssort); free(mbhmps); free(ymps);
        free(targsort); free(sourcesort);
        free(itargse); free(itarg_arr); free(isrcse); free(isrc);
        free(tcenters); free(boxsize); free(itree);
        return;
    }

    /* call main FMM */
    FNAME(mbhfmm2dmain)(&nd, eps,
        &beta, &ns, sourcesort,
        &ntermsmps, mbhmpssort, ympssort,
        &nt, targsort,
        iaddr, rmlexp, carray, &ldc,
        itree, &ltree, iptr, &ndiv, &nlevels,
        &nboxes, &iper_v, boxsize, rscales, tcenters,
        itree + iptr[0] - 1,
        isrcse, itargse, nterms_arr,
        &ifpgh, potsort, gradsort, hesssort,
        &ifpghtarg, pottargsort, gradtargsort, hesstargsort,
        ier);

    if (*ier != 0) goto cleanup;

    /* reorder output back */
    if (ifpgh == 1) {
        dreorderi_(&nd, &ns, potsort, pot, isrc);
    }
    if (ifpgh == 2) {
        dreorderi_(&nd, &ns, potsort, pot, isrc);
        fint tnd = 2 * nd;
        dreorderi_(&tnd, &ns, gradsort, grad, isrc);
    }
    if (ifpgh == 3) {
        dreorderi_(&nd, &ns, potsort, pot, isrc);
        fint tnd = 2 * nd;
        dreorderi_(&tnd, &ns, gradsort, grad, isrc);
        tnd = 3 * nd;
        dreorderi_(&tnd, &ns, hesssort, hess, isrc);
    }

    if (ifpghtarg == 1) {
        dreorderi_(&nd, &nt, pottargsort, pottarg, itarg_arr);
    }
    if (ifpghtarg == 2) {
        dreorderi_(&nd, &nt, pottargsort, pottarg, itarg_arr);
        fint tnd = 2 * nd;
        dreorderi_(&tnd, &nt, gradtargsort, gradtarg, itarg_arr);
    }
    if (ifpghtarg == 3) {
        dreorderi_(&nd, &nt, pottargsort, pottarg, itarg_arr);
        fint tnd = 2 * nd;
        dreorderi_(&tnd, &nt, gradtargsort, gradtarg, itarg_arr);
        tnd = 3 * nd;
        dreorderi_(&tnd, &nt, hesstargsort, hesstarg, itarg_arr);
    }

cleanup:
    free(rmlexp);
    free(iaddr); free(carray); free(nterms_arr); free(rscales);
    free(pottargsort); free(gradtargsort); free(hesstargsort);
    free(potsort); free(gradsort); free(hesssort);
    free(mbhmpssort); free(ympssort); free(mbhmps); free(ymps);
    free(targsort); free(sourcesort);
    free(itargse); free(itarg_arr); free(isrcse); free(isrc);
    free(tcenters); free(boxsize); free(itree);
}
