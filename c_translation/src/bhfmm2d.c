/*
 * bhfmm2d.c - C translation of src/biharmonic/bhfmm2d.f
 *
 * Top-level user-facing 2D biharmonic (Stokes) FMM driver and the
 * four subroutines from the Fortran source. Translated 1:1 from the
 * Fortran reference: same control flow, same allocations, same
 * operation order, same loop bounds. OpenMP directives are stripped
 * (sequential C). Print/log/error functions (prinf/prin2) are guarded
 * in the Fortran source by ifprint=0 and never executed; they are
 * omitted entirely. Timing calls (cpu_time, omp_get_wtime, second)
 * are also stripped: time1/time2/tt1/tt2 are left at zero so the
 * timeinfo / timelev outputs become deterministic zero values.
 *
 * Cross-file dependencies (bhndiv2d, pts_tree_*, dreorder*, init_carray,
 * bh2dterms, computemnlists, computelists, bh2dmpzero, bh2dformmp*,
 * bh2dformta*, bh2dmpeval*, bh2dtaeval*, bh2dmpmp, bh2dmploc, bh2dlocloc,
 * bh2d_direct*) are issued via bare Fortran symbol names so the diff
 * test isolates this file. Same-file calls (bhfmm2d -> bhfmm2dmain /
 * bh2dmpalloc, bhfmm2dmain -> bhfmm2dpart_direct) use the FNAME()
 * dispatcher.
 *
 * Differences from cfmm2d.c:
 *   - charge has shape (nd, 2, *) — leading dim 2*nd
 *   - dip has shape (nd, 3, *)    — leading dim 3*nd
 *   - mpole/local have FIVE coefficient sets per density:
 *       size per box = (nterms+1) * 5 * nd * 2 doubles
 *     (the *2 is because complex *16 in rmlexp is two doubles)
 *   - pot has shape (nd, *), grad has shape (nd, 3, *)
 *   - bhfmm2d only supports ifpgh up to 2 (no hessian)
 *   - There are no expansion centers, no list-3 mp-eval into expansions,
 *     no jsort/scjsort beyond placeholders. The Fortran source still
 *     references them; we preserve them as zero-initialized scratch.
 */

#include "bhfmm2d.h"

/* ---------------- cross-file extern declarations ---------------- */

extern void bhndiv2d_(const double *eps, const fint *ns, const fint *nt,
                      const fint *ifcharge, const fint *ifdipole,
                      const fint *ifpgh, const fint *ifpghtarg,
                      fint *ndiv, fint *idivflag);

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
extern void init_carray_(double *carray, const fint *ldc);

/* bh2dterms. */
extern void bh2dterms_(const double *eps, fint *nterms, fint *ier);

/* tree_routs2d. */
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

/* bhrouts2d. */
extern void bh2dmpzero_(const fint *nd, fcomplex *mpole, const fint *nterms);

extern void bh2dformmpc_(const fint *nd, const double *rscale,
                         const double *source, const fint *ns,
                         const fcomplex *charge, const double *center,
                         const fint *nterms, fcomplex *mpole);
extern void bh2dformmpd_(const fint *nd, const double *rscale,
                         const double *source, const fint *ns,
                         const fcomplex *dip, const double *center,
                         const fint *nterms, fcomplex *mpole);
extern void bh2dformmpcd_(const fint *nd, const double *rscale,
                          const double *source, const fint *ns,
                          const fcomplex *charge, const fcomplex *dip,
                          const double *center, const fint *nterms,
                          fcomplex *mpole);

extern void bh2dformtac_(const fint *nd, const double *rscale,
                         const double *source, const fint *ns,
                         const fcomplex *charge, const double *center,
                         const fint *nterms, fcomplex *local);
extern void bh2dformtad_(const fint *nd, const double *rscale,
                         const double *source, const fint *ns,
                         const fcomplex *dip, const double *center,
                         const fint *nterms, fcomplex *local);
extern void bh2dformtacd_(const fint *nd, const double *rscale,
                          const double *source, const fint *ns,
                          const fcomplex *charge, const fcomplex *dip,
                          const double *center, const fint *nterms,
                          fcomplex *local);

extern void bh2dmpevalp_(const fint *nd, const double *rscale,
                         const double *center, const fcomplex *mpole,
                         const fint *nterms, const double *ztarg,
                         const fint *ntarg, fcomplex *vel);
extern void bh2dmpevalg_(const fint *nd, const double *rscale,
                         const double *center, const fcomplex *mpole,
                         const fint *nterms, const double *ztarg,
                         const fint *ntarg, fcomplex *vel, fcomplex *grad);

extern void bh2dtaevalp_(const fint *nd, const double *rscale,
                         const double *center, const fcomplex *local,
                         const fint *nterms, const double *ztarg,
                         const fint *ntarg, fcomplex *vel);
extern void bh2dtaevalg_(const fint *nd, const double *rscale,
                         const double *center, const fcomplex *local,
                         const fint *nterms, const double *ztarg,
                         const fint *ntarg, fcomplex *vel, fcomplex *grad);

extern void bh2dmpmp_(const fint *nd,
                      const double *rscale1, const double *c1,
                      const fcomplex *hexp, const fint *nterms1,
                      const double *rscale2, const double *c2,
                      fcomplex *jexp, const fint *nterms2,
                      const double *carray, const fint *ldc);

extern void bh2dmploc_(const fint *nd,
                       const double *rscale1, const double *c1,
                       const fcomplex *hexp, const fint *nterms1,
                       const double *rscale2, const double *c2,
                       fcomplex *jexp, const fint *nterms2,
                       const double *carray, const fint *ldc);

extern void bh2dlocloc_(const fint *nd,
                        const double *rscale1, const double *c1,
                        const fcomplex *hexp, const fint *nterms1,
                        const double *rscale2, const double *c2,
                        fcomplex *jexp, const fint *nterms2,
                        const double *carray, const fint *ldc);

/* bhkernels2d. */
extern void bh2d_directcp_(const fint *nd, const double *sources,
                           const fint *ns, const fcomplex *charge,
                           const double *targ, const fint *nt,
                           fcomplex *vel, const double *thresh);
extern void bh2d_directcg_(const fint *nd, const double *sources,
                           const fint *ns, const fcomplex *charge,
                           const double *targ, const fint *nt,
                           fcomplex *vel, fcomplex *grad,
                           const double *thresh);
extern void bh2d_directdp_(const fint *nd, const double *sources,
                           const fint *ns, const fcomplex *dip,
                           const double *targ, const fint *nt,
                           fcomplex *vel, const double *thresh);
extern void bh2d_directdg_(const fint *nd, const double *sources,
                           const fint *ns, const fcomplex *dip,
                           const double *targ, const fint *nt,
                           fcomplex *vel, fcomplex *grad,
                           const double *thresh);
extern void bh2d_directcdp_(const fint *nd, const double *sources,
                            const fint *ns, const fcomplex *charge,
                            const fcomplex *dip, const double *targ,
                            const fint *nt, fcomplex *vel,
                            const double *thresh);
extern void bh2d_directcdg_(const fint *nd, const double *sources,
                            const fint *ns, const fcomplex *charge,
                            const fcomplex *dip, const double *targ,
                            const fint *nt, fcomplex *vel, fcomplex *grad,
                            const double *thresh);


/*
 * laddr(k, ilev): k is 1-based, ilev is 0-based, leading dim 2.
 */
#define LADDR(k, ilev) ((ilev) * 2 + ((k) - 1))


/* ---------------------------------------------------------------- */
/* bh2dmpalloc - lay out workspace for multipole and local expansions */
/* ---------------------------------------------------------------- */
void FNAME(bh2dmpalloc)(const fint *nd, const fint *laddr, fint *iaddr,
                        const fint *nlevels, fint *lmptot,
                        const fint *nterms)
{
    fint nd_v = *nd;
    fint nlevels_v = *nlevels;
    fint istart, i, ibox, nn, itmp;

    /* Each multipole/local expansion has 5 coefficient sets per
     * density. complex *16 = 2 doubles. So nn doubles per box. */
    istart = 1;
    for (i = 0; i <= nlevels_v; i++) {
        nn = (nterms[i] + 1) * 2 * nd_v * 5;
        for (ibox = laddr[LADDR(1, i)]; ibox <= laddr[LADDR(2, i)]; ibox++) {
            itmp = ibox - laddr[LADDR(1, i)];
            iaddr[FA2(1, ibox, 2)] = istart + itmp * nn;
        }
        istart = istart + (laddr[LADDR(2, i)] - laddr[LADDR(1, i)] + 1) * nn;
    }

    for (i = 0; i <= nlevels_v; i++) {
        nn = (nterms[i] + 1) * 2 * nd_v * 5;
        for (ibox = laddr[LADDR(1, i)]; ibox <= laddr[LADDR(2, i)]; ibox++) {
            itmp = ibox - laddr[LADDR(1, i)];
            iaddr[FA2(2, ibox, 2)] = istart + itmp * nn;
        }
        istart = istart + (laddr[LADDR(2, i)] - laddr[LADDR(1, i)] + 1) * nn;
    }

    *lmptot = istart;
}


/* ---------------------------------------------------------------- */
/* bhfmm2dpart_direct - direct particle-to-particle dispatcher       */
/* ---------------------------------------------------------------- */
void FNAME(bhfmm2dpart_direct)(const fint *nd,
                               const fint *istart, const fint *iend,
                               const fint *jstart, const fint *jend,
                               const double *source,
                               const fint *ifcharge, const fcomplex *charge,
                               const fint *ifdipole, const fcomplex *dip,
                               const double *targ, const fint *ifpgh,
                               fcomplex *pot, fcomplex *grad, fcomplex *hess,
                               const double *thresh)
{
    fint nd_v = *nd;
    fint istart_v = *istart;
    fint iend_v = *iend;
    fint jstart_v = *jstart;
    fint jend_v = *jend;
    fint ifcharge_v = *ifcharge;
    fint ifdipole_v = *ifdipole;
    fint ifpgh_v = *ifpgh;
    fint ns, nt;

    /* charge(nd, 2, *): leading dim is 2*nd
     * dip(nd, 3, *):    leading dim is 3*nd
     * pot(nd, *):       leading dim is nd
     * grad(nd, 3, *):   leading dim is 3*nd
     * source(2, *):     leading dim is 2
     * targ(2, *):       leading dim is 2 */
    const double *src_p = &source[(istart_v - 1) * 2];
    const fcomplex *chg_p = &charge[(istart_v - 1) * 2 * nd_v];
    const fcomplex *dip_p = &dip[(istart_v - 1) * 3 * nd_v];
    const double *tg_p = &targ[(jstart_v - 1) * 2];
    fcomplex *pot_p = &pot[(jstart_v - 1) * nd_v];
    fcomplex *grad_p = &grad[(jstart_v - 1) * 3 * nd_v];
    (void)hess;

    ns = iend_v - istart_v + 1;
    nt = jend_v - jstart_v + 1;

    if (ifcharge_v == 1 && ifdipole_v == 0) {
        if (ifpgh_v == 1) {
            bh2d_directcp_(&nd_v, src_p, &ns, chg_p, tg_p, &nt,
                           pot_p, thresh);
        }
        if (ifpgh_v == 2) {
            bh2d_directcg_(&nd_v, src_p, &ns, chg_p, tg_p, &nt,
                           pot_p, grad_p, thresh);
        }
    }

    if (ifcharge_v == 0 && ifdipole_v == 1) {
        if (ifpgh_v == 1) {
            bh2d_directdp_(&nd_v, src_p, &ns, dip_p, tg_p, &nt,
                           pot_p, thresh);
        }
        if (ifpgh_v == 2) {
            bh2d_directdg_(&nd_v, src_p, &ns, dip_p, tg_p, &nt,
                           pot_p, grad_p, thresh);
        }
    }

    if (ifcharge_v == 1 && ifdipole_v == 1) {
        if (ifpgh_v == 1) {
            bh2d_directcdp_(&nd_v, src_p, &ns, chg_p, dip_p, tg_p, &nt,
                            pot_p, thresh);
        }
        if (ifpgh_v == 2) {
            bh2d_directcdg_(&nd_v, src_p, &ns, chg_p, dip_p, tg_p, &nt,
                            pot_p, grad_p, thresh);
        }
    }
}


/* ---------------------------------------------------------------- */
/* bhfmm2dmain - the main FMM engine (8 steps)                       */
/* ---------------------------------------------------------------- */
void FNAME(bhfmm2dmain)(const fint *nd, const double *eps,
                        const fint *nsource, const double *sourcesort,
                        const fint *ifcharge, const fcomplex *chargesort,
                        const fint *ifdipole, const fcomplex *dipsort,
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
                        fcomplex *jsort, double *scjsort)
{
    fint nd_v = *nd;
    fint nexpc_v = *nexpc;
    fint nlevels_v = *nlevels;
    fint nboxes_v = *nboxes;
    fint iper_v = *iper;
    fint ifcharge_v = *ifcharge;
    fint ifdipole_v = *ifdipole;
    fint ntj_v = *ntj;
    fint ifpgh_v = *ifpgh;
    fint ifpghtarg_v = *ifpghtarg;

    fint i, j, ilev, idim;
    fint ibox, jbox, nchild, npts;
    fint istart, iend;
    fint istarts, iends, istartt, iendt;
    fint jstart, jend;

    fint mnlist1 = 0, mnlist2 = 0, mnlist3 = 0, mnlist4 = 0;

    fint *list1 = NULL, *list2 = NULL, *list3 = NULL, *list4 = NULL;
    fint *nlist1s = NULL, *nlist2s = NULL, *nlist3s = NULL, *nlist4s = NULL;

    double *carray = NULL;
    fint ldc;

    double timelev[201];
    double thresh;
    double timeinfo[10];

    /* Suppress unused warnings on parameters that are passed through. */
    (void)eps;
    (void)nsource;
    (void)ntarget;
    (void)mptemp;
    (void)lmptmp;
    (void)ltree;
    (void)ndiv;
    (void)expcsort;
    (void)hess;
    (void)hesstarg;
    (void)timelev;
    (void)timeinfo;

    /* time1 / time2 / tt1 / tt2: kept as zero so timeinfo entries are zero. */

    /* timelev(0:nlevels) zeroed. */
    for (i = 0; i <= nlevels_v; i++) {
        timelev[i] = 0.0;
    }

    ldc = 100;
    /* carray(0:ldc, 0:ldc) - leading dim (ldc+1), total (ldc+1)^2. */
    carray = (double *)malloc((ldc + 1) * (ldc + 1) * sizeof(double));
    init_carray_(carray, &ldc);

    /* compute list info */
    computemnlists_(&nlevels_v, &nboxes_v, itree, ltree, iptr, centers,
                    boxsize, &iper_v, &mnlist1, &mnlist2, &mnlist3, &mnlist4);

    nlist1s = (fint *)malloc(nboxes_v * sizeof(fint));
    list1 = (fint *)malloc(mnlist1 * nboxes_v * sizeof(fint));
    nlist2s = (fint *)malloc(nboxes_v * sizeof(fint));
    list2 = (fint *)malloc(mnlist2 * nboxes_v * sizeof(fint));
    nlist3s = (fint *)malloc(nboxes_v * sizeof(fint));
    list3 = (fint *)malloc(mnlist3 * nboxes_v * sizeof(fint));
    nlist4s = (fint *)malloc(nboxes_v * sizeof(fint));
    list4 = (fint *)malloc(mnlist4 * nboxes_v * sizeof(fint));

    computelists_(&nlevels_v, &nboxes_v, itree, ltree, iptr, centers,
                  boxsize, &iper_v,
                  &mnlist1, nlist1s, list1,
                  &mnlist2, nlist2s, list2,
                  &mnlist3, nlist3s, list3,
                  &mnlist4, nlist4s, list4);

    /*
     * Set the expansion coefficients (jsort) to zero.
     * jsort(nd, 5, 0:ntj, *) — but in this code path nexpc=0, so the
     * loop body never runs. Preserve the structure for fidelity.
     */
    for (i = 1; i <= nexpc_v; i++) {
        for (j = 0; j <= ntj_v; j++) {
            for (idim = 1; idim <= nd_v; idim++) {
                /* jsort(idim,1,j,i) ... jsort(idim,5,j,i) all zero */
                fint base = (((i - 1) * (ntj_v + 1) + j) * 5) * nd_v
                          + (idim - 1);
                jsort[base + 0 * nd_v] = 0;
                jsort[base + 1 * nd_v] = 0;
                jsort[base + 2 * nd_v] = 0;
                jsort[base + 3 * nd_v] = 0;
                jsort[base + 4 * nd_v] = 0;
            }
        }
    }

    for (i = 0; i < 10; i++) {
        timeinfo[i] = 0.0;
    }

    /* set all multipole and local expansions to zero */
    for (ilev = 0; ilev <= nlevels_v; ilev++) {
        for (ibox = laddr[LADDR(1, ilev)];
             ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            bh2dmpzero_(&nd_v,
                        (fcomplex *)&rmlexp[iaddr[FA2(1, ibox, 2)] - 1],
                        &nterms[ilev]);
            bh2dmpzero_(&nd_v,
                        (fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 2)] - 1],
                        &nterms[ilev]);
        }
    }

    /* Set scjsort - but nexpc = 0 so the inner loop never runs.
     * Preserve the structure for fidelity. */
    for (ilev = 0; ilev <= nlevels_v; ilev++) {
        for (ibox = laddr[LADDR(1, ilev)];
             ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            nchild = itree[iptr[3] + ibox - 1 - 1];
            if (nchild == 0) {
                istart = iexpcse[FA2(1, ibox, 2)];
                iend = iexpcse[FA2(2, ibox, 2)];
                for (i = istart; i <= iend; i++) {
                    scjsort[i - 1] = rscales[ilev];
                }
            }
        }
    }

    /* ============================================================ */
    /* Step 1: form mp                                              */
    /* ============================================================ */
    for (ilev = 2; ilev <= nlevels_v; ilev++) {
        if (ifcharge_v == 1 && ifdipole_v == 0) {
            for (ibox = laddr[LADDR(1, ilev)];
                 ibox <= laddr[LADDR(2, ilev)]; ibox++) {
                nchild = itree[iptr[3] + ibox - 1 - 1];
                istart = isrcse[FA2(1, ibox, 2)];
                iend = isrcse[FA2(2, ibox, 2)];
                npts = iend - istart + 1;
                if (nchild == 0 && npts > 0) {
                    bh2dformmpc_(&nd_v, &rscales[ilev],
                                 &sourcesort[(istart - 1) * 2], &npts,
                                 &chargesort[(istart - 1) * 2 * nd_v],
                                 &centers[(ibox - 1) * 2], &nterms[ilev],
                                 (fcomplex *)&rmlexp[iaddr[FA2(1, ibox, 2)] - 1]);
                }
            }
        }

        if (ifdipole_v == 1 && ifcharge_v == 0) {
            for (ibox = laddr[LADDR(1, ilev)];
                 ibox <= laddr[LADDR(2, ilev)]; ibox++) {
                nchild = itree[iptr[3] + ibox - 1 - 1];
                istart = isrcse[FA2(1, ibox, 2)];
                iend = isrcse[FA2(2, ibox, 2)];
                npts = iend - istart + 1;
                if (nchild == 0 && npts > 0) {
                    bh2dformmpd_(&nd_v, &rscales[ilev],
                                 &sourcesort[(istart - 1) * 2], &npts,
                                 &dipsort[(istart - 1) * 3 * nd_v],
                                 &centers[(ibox - 1) * 2], &nterms[ilev],
                                 (fcomplex *)&rmlexp[iaddr[FA2(1, ibox, 2)] - 1]);
                }
            }
        }

        if (ifdipole_v == 1 && ifcharge_v == 1) {
            for (ibox = laddr[LADDR(1, ilev)];
                 ibox <= laddr[LADDR(2, ilev)]; ibox++) {
                nchild = itree[iptr[3] + ibox - 1 - 1];
                istart = isrcse[FA2(1, ibox, 2)];
                iend = isrcse[FA2(2, ibox, 2)];
                npts = iend - istart + 1;
                if (nchild == 0 && npts > 0) {
                    bh2dformmpcd_(&nd_v, &rscales[ilev],
                                  &sourcesort[(istart - 1) * 2], &npts,
                                  &chargesort[(istart - 1) * 2 * nd_v],
                                  &dipsort[(istart - 1) * 3 * nd_v],
                                  &centers[(ibox - 1) * 2], &nterms[ilev],
                                  (fcomplex *)&rmlexp[iaddr[FA2(1, ibox, 2)] - 1]);
                }
            }
        }
    }
    timeinfo[0] = 0.0;

    /* ============================================================ */
    /* Step 2: form lo (form local expansions from list 4)          */
    /* ============================================================ */
    for (ilev = 2; ilev <= nlevels_v; ilev++) {
        if (ifcharge_v == 1 && ifdipole_v == 0) {
            for (ibox = laddr[LADDR(1, ilev)];
                 ibox <= laddr[LADDR(2, ilev)]; ibox++) {
                npts = 0;
                if (ifpghtarg_v > 0) {
                    istart = itargse[FA2(1, ibox, 2)];
                    iend = itargse[FA2(2, ibox, 2)];
                    npts = npts + iend - istart + 1;
                }
                istart = iexpcse[FA2(1, ibox, 2)];
                iend = iexpcse[FA2(2, ibox, 2)];
                npts = npts + iend - istart + 1;
                if (ifpgh_v > 0) {
                    istart = isrcse[FA2(1, ibox, 2)];
                    iend = isrcse[FA2(2, ibox, 2)];
                    npts = npts + iend - istart + 1;
                }

                if (npts > 0) {
                    for (i = 1; i <= nlist4s[ibox - 1]; i++) {
                        jbox = list4[FA2(i, ibox, mnlist4)];
                        istart = isrcse[FA2(1, jbox, 2)];
                        iend = isrcse[FA2(2, jbox, 2)];
                        npts = iend - istart + 1;
                        bh2dformtac_(&nd_v, &rscales[ilev],
                                     &sourcesort[(istart - 1) * 2], &npts,
                                     &chargesort[(istart - 1) * 2 * nd_v],
                                     &centers[(ibox - 1) * 2], &nterms[ilev],
                                     (fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 2)] - 1]);
                    }
                }
            }
        }

        if (ifcharge_v == 0 && ifdipole_v == 1) {
            for (ibox = laddr[LADDR(1, ilev)];
                 ibox <= laddr[LADDR(2, ilev)]; ibox++) {
                npts = 0;
                if (ifpghtarg_v > 0) {
                    istart = itargse[FA2(1, ibox, 2)];
                    iend = itargse[FA2(2, ibox, 2)];
                    npts = npts + iend - istart + 1;
                }
                istart = iexpcse[FA2(1, ibox, 2)];
                iend = iexpcse[FA2(2, ibox, 2)];
                npts = npts + iend - istart + 1;
                if (ifpgh_v > 0) {
                    istart = isrcse[FA2(1, ibox, 2)];
                    iend = isrcse[FA2(2, ibox, 2)];
                    npts = npts + iend - istart + 1;
                }

                if (npts > 0) {
                    for (i = 1; i <= nlist4s[ibox - 1]; i++) {
                        jbox = list4[FA2(i, ibox, mnlist4)];
                        istart = isrcse[FA2(1, jbox, 2)];
                        iend = isrcse[FA2(2, jbox, 2)];
                        npts = iend - istart + 1;
                        bh2dformtad_(&nd_v, &rscales[ilev],
                                     &sourcesort[(istart - 1) * 2], &npts,
                                     &dipsort[(istart - 1) * 3 * nd_v],
                                     &centers[(ibox - 1) * 2], &nterms[ilev],
                                     (fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 2)] - 1]);
                    }
                }
            }
        }

        if (ifcharge_v == 1 && ifdipole_v == 1) {
            for (ibox = laddr[LADDR(1, ilev)];
                 ibox <= laddr[LADDR(2, ilev)]; ibox++) {
                npts = 0;
                if (ifpghtarg_v > 0) {
                    istart = itargse[FA2(1, ibox, 2)];
                    iend = itargse[FA2(2, ibox, 2)];
                    npts = npts + iend - istart + 1;
                }
                istart = iexpcse[FA2(1, ibox, 2)];
                iend = iexpcse[FA2(2, ibox, 2)];
                npts = npts + iend - istart + 1;
                if (ifpgh_v > 0) {
                    istart = isrcse[FA2(1, ibox, 2)];
                    iend = isrcse[FA2(2, ibox, 2)];
                    npts = npts + iend - istart + 1;
                }

                if (npts > 0) {
                    for (i = 1; i <= nlist4s[ibox - 1]; i++) {
                        jbox = list4[FA2(i, ibox, mnlist4)];
                        istart = isrcse[FA2(1, jbox, 2)];
                        iend = isrcse[FA2(2, jbox, 2)];
                        npts = iend - istart + 1;
                        bh2dformtacd_(&nd_v, &rscales[ilev],
                                      &sourcesort[(istart - 1) * 2], &npts,
                                      &chargesort[(istart - 1) * 2 * nd_v],
                                      &dipsort[(istart - 1) * 3 * nd_v],
                                      &centers[(ibox - 1) * 2], &nterms[ilev],
                                      (fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 2)] - 1]);
                    }
                }
            }
        }
    }
    timeinfo[1] = 0.0;

    /* ============================================================ */
    /* Step 3: merge mp (M2M up the tree)                           */
    /* ============================================================ */
    for (ilev = nlevels_v - 1; ilev >= 1; ilev--) {
        for (ibox = laddr[LADDR(1, ilev)];
             ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            nchild = itree[iptr[3] + ibox - 1 - 1];
            for (i = 1; i <= nchild; i++) {
                /* itree(iptr(5) + 4*(ibox-1) + i - 1) */
                jbox = itree[iptr[4] + 4 * (ibox - 1) + i - 1 - 1];
                istart = isrcse[FA2(1, jbox, 2)];
                iend = isrcse[FA2(2, jbox, 2)];
                npts = iend - istart + 1;
                if (npts > 0) {
                    fint nterms_p1 = nterms[ilev + 1];
                    fint nterms_l = nterms[ilev];
                    bh2dmpmp_(&nd_v, &rscales[ilev + 1],
                              &centers[(jbox - 1) * 2],
                              (fcomplex *)&rmlexp[iaddr[FA2(1, jbox, 2)] - 1],
                              &nterms_p1, &rscales[ilev],
                              &centers[(ibox - 1) * 2],
                              (fcomplex *)&rmlexp[iaddr[FA2(1, ibox, 2)] - 1],
                              &nterms_l, carray, &ldc);
                }
            }
        }
    }
    timeinfo[2] = 0.0;

    /* ============================================================ */
    /* Step 4: mp to loc (list 2 / M2L)                             */
    /* ============================================================ */
    for (ilev = 2; ilev <= nlevels_v; ilev++) {
        for (ibox = laddr[LADDR(1, ilev)];
             ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            npts = 0;
            if (ifpghtarg_v > 0) {
                istart = itargse[FA2(1, ibox, 2)];
                iend = itargse[FA2(2, ibox, 2)];
                npts = npts + iend - istart + 1;
            }
            istart = iexpcse[FA2(1, ibox, 2)];
            iend = iexpcse[FA2(2, ibox, 2)];
            npts = npts + iend - istart + 1;
            if (ifpgh_v > 0) {
                istart = isrcse[FA2(1, ibox, 2)];
                iend = isrcse[FA2(2, ibox, 2)];
                npts = npts + iend - istart + 1;
            }
            if (npts > 0) {
                for (i = 1; i <= nlist2s[ibox - 1]; i++) {
                    fint nterms_l = nterms[ilev];
                    jbox = list2[FA2(i, ibox, mnlist2)];
                    bh2dmploc_(&nd_v, &rscales[ilev],
                               &centers[(jbox - 1) * 2],
                               (fcomplex *)&rmlexp[iaddr[FA2(1, jbox, 2)] - 1],
                               &nterms_l, &rscales[ilev],
                               &centers[(ibox - 1) * 2],
                               (fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 2)] - 1],
                               &nterms_l, carray, &ldc);
                }
            }
        }
        timelev[ilev] = 0.0;
    }
    timeinfo[3] = 0.0;

    /* ============================================================ */
    /* Step 5: split loc (L2L down the tree)                        */
    /* ============================================================ */
    for (ilev = 1; ilev <= nlevels_v - 1; ilev++) {
        for (ibox = laddr[LADDR(1, ilev)];
             ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            nchild = itree[iptr[3] + ibox - 1 - 1];
            istart = iexpcse[FA2(1, ibox, 2)];
            iend = iexpcse[FA2(2, ibox, 2)];
            npts = iend - istart + 1;

            if (ifpghtarg_v > 0) {
                istart = itargse[FA2(1, ibox, 2)];
                iend = itargse[FA2(2, ibox, 2)];
                npts = npts + iend - istart + 1;
            }
            if (ifpgh_v > 0) {
                istart = isrcse[FA2(1, ibox, 2)];
                iend = isrcse[FA2(2, ibox, 2)];
                npts = npts + iend - istart + 1;
            }

            if (npts > 0) {
                for (i = 1; i <= nchild; i++) {
                    fint nterms_l = nterms[ilev];
                    fint nterms_p1 = nterms[ilev + 1];
                    jbox = itree[iptr[4] + 4 * (ibox - 1) + i - 1 - 1];
                    bh2dlocloc_(&nd_v, &rscales[ilev],
                                &centers[(ibox - 1) * 2],
                                (fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 2)] - 1],
                                &nterms_l, &rscales[ilev + 1],
                                &centers[(jbox - 1) * 2],
                                (fcomplex *)&rmlexp[iaddr[FA2(2, jbox, 2)] - 1],
                                &nterms_p1, carray, &ldc);
                }
            }
        }
    }
    timeinfo[4] = 0.0;

    /* ============================================================ */
    /* Step 6: mp eval (list 3)                                     */
    /* ============================================================ */
    for (ilev = 1; ilev <= nlevels_v - 1; ilev++) {
        for (ibox = laddr[LADDR(1, ilev)];
             ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            /* evaluate multipole expansion at all targets */
            istart = itargse[FA2(1, ibox, 2)];
            iend = itargse[FA2(2, ibox, 2)];
            npts = iend - istart + 1;

            if (ifpghtarg_v == 1) {
                for (i = 1; i <= nlist3s[ibox - 1]; i++) {
                    fint nterms_p1 = nterms[ilev + 1];
                    jbox = list3[FA2(i, ibox, mnlist3)];
                    bh2dmpevalp_(&nd_v, &rscales[ilev + 1],
                                 &centers[(jbox - 1) * 2],
                                 (const fcomplex *)&rmlexp[iaddr[FA2(1, jbox, 2)] - 1],
                                 &nterms_p1, &targetsort[(istart - 1) * 2],
                                 &npts,
                                 &pottarg[(istart - 1) * nd_v]);
                }
            }
            if (ifpghtarg_v == 2) {
                for (i = 1; i <= nlist3s[ibox - 1]; i++) {
                    fint nterms_p1 = nterms[ilev + 1];
                    jbox = list3[FA2(i, ibox, mnlist3)];
                    bh2dmpevalg_(&nd_v, &rscales[ilev + 1],
                                 &centers[(jbox - 1) * 2],
                                 (const fcomplex *)&rmlexp[iaddr[FA2(1, jbox, 2)] - 1],
                                 &nterms_p1, &targetsort[(istart - 1) * 2],
                                 &npts,
                                 &pottarg[(istart - 1) * nd_v],
                                 &gradtarg[(istart - 1) * 3 * nd_v]);
                }
            }

            /* evaluate multipole expansion at all sources */
            istart = isrcse[FA2(1, ibox, 2)];
            iend = isrcse[FA2(2, ibox, 2)];
            npts = iend - istart + 1;

            if (ifpgh_v == 1) {
                for (i = 1; i <= nlist3s[ibox - 1]; i++) {
                    fint nterms_p1 = nterms[ilev + 1];
                    jbox = list3[FA2(i, ibox, mnlist3)];
                    bh2dmpevalp_(&nd_v, &rscales[ilev + 1],
                                 &centers[(jbox - 1) * 2],
                                 (const fcomplex *)&rmlexp[iaddr[FA2(1, jbox, 2)] - 1],
                                 &nterms_p1, &sourcesort[(istart - 1) * 2],
                                 &npts,
                                 &pot[(istart - 1) * nd_v]);
                }
            }
            if (ifpgh_v == 2) {
                for (i = 1; i <= nlist3s[ibox - 1]; i++) {
                    fint nterms_p1 = nterms[ilev + 1];
                    jbox = list3[FA2(i, ibox, mnlist3)];
                    bh2dmpevalg_(&nd_v, &rscales[ilev + 1],
                                 &centers[(jbox - 1) * 2],
                                 (const fcomplex *)&rmlexp[iaddr[FA2(1, jbox, 2)] - 1],
                                 &nterms_p1, &sourcesort[(istart - 1) * 2],
                                 &npts,
                                 &pot[(istart - 1) * nd_v],
                                 &grad[(istart - 1) * 3 * nd_v]);
                }
            }
        }
    }
    /* label 1000 (unused goto target in Fortran source) */
    timeinfo[5] = 0.0;

    /* ============================================================ */
    /* Step 7: eval lo (evaluate local expansions)                  */
    /* ============================================================ */
    for (ilev = 0; ilev <= nlevels_v; ilev++) {
        for (ibox = laddr[LADDR(1, ilev)];
             ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            nchild = itree[iptr[3] + ibox - 1 - 1];
            if (nchild == 0) {
                /* evaluate local expansion at targets */
                istart = itargse[FA2(1, ibox, 2)];
                iend = itargse[FA2(2, ibox, 2)];
                npts = iend - istart + 1;
                if (ifpghtarg_v == 1) {
                    fint nterms_l = nterms[ilev];
                    bh2dtaevalp_(&nd_v, &rscales[ilev],
                                 &centers[(ibox - 1) * 2],
                                 (const fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 2)] - 1],
                                 &nterms_l, &targetsort[(istart - 1) * 2],
                                 &npts,
                                 &pottarg[(istart - 1) * nd_v]);
                }
                if (ifpghtarg_v == 2) {
                    fint nterms_l = nterms[ilev];
                    bh2dtaevalg_(&nd_v, &rscales[ilev],
                                 &centers[(ibox - 1) * 2],
                                 (const fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 2)] - 1],
                                 &nterms_l, &targetsort[(istart - 1) * 2],
                                 &npts,
                                 &pottarg[(istart - 1) * nd_v],
                                 &gradtarg[(istart - 1) * 3 * nd_v]);
                }

                /* evaluate local expansion at sources */
                istart = isrcse[FA2(1, ibox, 2)];
                iend = isrcse[FA2(2, ibox, 2)];
                npts = iend - istart + 1;
                if (ifpgh_v == 1) {
                    fint nterms_l = nterms[ilev];
                    bh2dtaevalp_(&nd_v, &rscales[ilev],
                                 &centers[(ibox - 1) * 2],
                                 (const fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 2)] - 1],
                                 &nterms_l, &sourcesort[(istart - 1) * 2],
                                 &npts,
                                 &pot[(istart - 1) * nd_v]);
                }
                if (ifpgh_v == 2) {
                    fint nterms_l = nterms[ilev];
                    bh2dtaevalg_(&nd_v, &rscales[ilev],
                                 &centers[(ibox - 1) * 2],
                                 (const fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 2)] - 1],
                                 &nterms_l, &sourcesort[(istart - 1) * 2],
                                 &npts,
                                 &pot[(istart - 1) * nd_v],
                                 &grad[(istart - 1) * 3 * nd_v]);
                }
            }
        }
    }
    timeinfo[6] = 0.0;

    /* ============================================================ */
    /* Step 8: direct                                               */
    /* ============================================================ */
    /* thresh = boxsize(0) * 2^(-51) */
    thresh = boxsize[0] * 0x1p-51;

    for (ilev = 0; ilev <= nlevels_v; ilev++) {
        for (ibox = laddr[LADDR(1, ilev)];
             ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            istartt = itargse[FA2(1, ibox, 2)];
            iendt = itargse[FA2(2, ibox, 2)];

            istarts = isrcse[FA2(1, ibox, 2)];
            iends = isrcse[FA2(2, ibox, 2)];

            for (i = 1; i <= nlist1s[ibox - 1]; i++) {
                jbox = list1[FA2(i, ibox, mnlist1)];

                jstart = isrcse[FA2(1, jbox, 2)];
                jend = isrcse[FA2(2, jbox, 2)];

                FNAME(bhfmm2dpart_direct)(&nd_v, &jstart, &jend, &istartt,
                                          &iendt, sourcesort, &ifcharge_v,
                                          chargesort, &ifdipole_v, dipsort,
                                          targetsort, &ifpghtarg_v, pottarg,
                                          gradtarg, hesstarg, &thresh);

                FNAME(bhfmm2dpart_direct)(&nd_v, &jstart, &jend, &istarts,
                                          &iends, sourcesort, &ifcharge_v,
                                          chargesort, &ifdipole_v, dipsort,
                                          sourcesort, &ifpgh_v, pot, grad,
                                          hess, &thresh);
            }
        }
    }
    timeinfo[7] = 0.0;

    free(carray);
    free(nlist1s);
    free(list1);
    free(nlist2s);
    free(list2);
    free(nlist3s);
    free(list3);
    free(nlist4s);
    free(list4);
}


/* ---------------------------------------------------------------- */
/* bhfmm2d - top-level user-callable driver                          */
/* ---------------------------------------------------------------- */
void FNAME(bhfmm2d)(const fint *nd, const double *eps,
                    const fint *ns, const double *sources,
                    const fint *ifcharge, const fcomplex *charge,
                    const fint *ifdipole, const fcomplex *dip,
                    fint *iper, const fint *ifpgh,
                    fcomplex *pot, fcomplex *grad, fcomplex *hess,
                    const fint *nt, const double *targ,
                    const fint *ifpghtarg,
                    fcomplex *pottarg, fcomplex *gradtarg,
                    fcomplex *hesstarg, fint *ier)
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
    fint idivflag, nlevels, nboxes, ndiv;
    fint ltree;

    fint *isrc = NULL, *isrcse = NULL;
    fint *itarg = NULL, *itargse = NULL, *iexpcse = NULL;
    double *sourcesort = NULL;
    double *targsort = NULL;
    fcomplex *chargesort = NULL, *dipsort = NULL;
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

    *ier = 0;

    bhndiv2d_(eps, ns, nt, ifcharge, ifdipole, ifpgh, ifpghtarg,
              &ndiv, &idivflag);

    ltree = 0;
    nlmin = 0;
    nlmax = 51;
    ifunif = 0;
    *iper = 0;

    /* tree memory sizing */
    pts_tree_mem_(sources, ns, targ, nt, &idivflag, &ndiv, &nlmin, &nlmax,
                  &ifunif, iper, &nlevels, &nboxes, &ltree);

    itree = (fint *)malloc(ltree * sizeof(fint));
    boxsize = (double *)malloc((nlevels + 1) * sizeof(double));
    tcenters = (double *)malloc(2 * nboxes * sizeof(double));

    /* build the tree */
    pts_tree_build_(sources, ns, targ, nt, &idivflag, &ndiv, &nlmin, &nlmax,
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

    /* charge has shape (nd, 2, ns); dip has shape (nd, 3, ns) */
    if (ifcharge_v == 1 && ifdipole_v == 0) {
        chargesort = (fcomplex *)malloc(2 * nd_v * ns_v * sizeof(fcomplex));
        dipsort = (fcomplex *)malloc(3 * nd_v * 1 * sizeof(fcomplex));
    }
    if (ifcharge_v == 0 && ifdipole_v == 1) {
        chargesort = (fcomplex *)malloc(2 * nd_v * 1 * sizeof(fcomplex));
        dipsort = (fcomplex *)malloc(3 * nd_v * ns_v * sizeof(fcomplex));
    }
    if (ifcharge_v == 1 && ifdipole_v == 1) {
        chargesort = (fcomplex *)malloc(2 * nd_v * ns_v * sizeof(fcomplex));
        dipsort = (fcomplex *)malloc(3 * nd_v * ns_v * sizeof(fcomplex));
    }
    if (chargesort == NULL) {
        chargesort = (fcomplex *)malloc(2 * nd_v * 1 * sizeof(fcomplex));
    }
    if (dipsort == NULL) {
        dipsort = (fcomplex *)malloc(3 * nd_v * 1 * sizeof(fcomplex));
    }

    if (ifpgh_v == 1) {
        potsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        gradsort = (fcomplex *)malloc(3 * nd_v * 1 * sizeof(fcomplex));
        hesssort = (fcomplex *)malloc(3 * nd_v * 1 * sizeof(fcomplex));
    } else if (ifpgh_v == 2) {
        potsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        gradsort = (fcomplex *)malloc(3 * nd_v * ns_v * sizeof(fcomplex));
        hesssort = (fcomplex *)malloc(3 * nd_v * 1 * sizeof(fcomplex));
    } else if (ifpgh_v == 3) {
        potsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        gradsort = (fcomplex *)malloc(3 * nd_v * ns_v * sizeof(fcomplex));
        hesssort = (fcomplex *)malloc(3 * nd_v * ns_v * sizeof(fcomplex));
    } else {
        potsort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
        gradsort = (fcomplex *)malloc(3 * nd_v * 1 * sizeof(fcomplex));
        hesssort = (fcomplex *)malloc(3 * nd_v * 1 * sizeof(fcomplex));
    }

    if (ifpghtarg_v == 1) {
        pottargsort = (fcomplex *)malloc(nd_v * nt_v * sizeof(fcomplex));
        gradtargsort = (fcomplex *)malloc(3 * nd_v * 1 * sizeof(fcomplex));
        hesstargsort = (fcomplex *)malloc(3 * nd_v * 1 * sizeof(fcomplex));
    } else if (ifpghtarg_v == 2) {
        pottargsort = (fcomplex *)malloc(nd_v * nt_v * sizeof(fcomplex));
        gradtargsort = (fcomplex *)malloc(3 * nd_v * nt_v * sizeof(fcomplex));
        hesstargsort = (fcomplex *)malloc(3 * nd_v * 1 * sizeof(fcomplex));
    } else if (ifpghtarg_v == 3) {
        pottargsort = (fcomplex *)malloc(nd_v * nt_v * sizeof(fcomplex));
        gradtargsort = (fcomplex *)malloc(3 * nd_v * nt_v * sizeof(fcomplex));
        hesstargsort = (fcomplex *)malloc(3 * nd_v * nt_v * sizeof(fcomplex));
    } else {
        pottargsort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
        gradtargsort = (fcomplex *)malloc(3 * nd_v * 1 * sizeof(fcomplex));
        hesstargsort = (fcomplex *)malloc(3 * nd_v * 1 * sizeof(fcomplex));
    }

    /* initialize potentials, gradients, hessians */
    if (ifpgh_v == 1) {
        for (i = 1; i <= ns_v; i++) {
            for (idim = 1; idim <= nd_v; idim++) {
                potsort[FA2(idim, i, nd_v)] = 0;
            }
        }
    }
    if (ifpgh_v == 2) {
        for (i = 1; i <= ns_v; i++) {
            for (idim = 1; idim <= nd_v; idim++) {
                potsort[FA2(idim, i, nd_v)] = 0;
                gradsort[FA3(idim, 1, i, nd_v, 3)] = 0;
                gradsort[FA3(idim, 2, i, nd_v, 3)] = 0;
                gradsort[FA3(idim, 3, i, nd_v, 3)] = 0;
            }
        }
    }
    if (ifpgh_v == 3) {
        for (i = 1; i <= ns_v; i++) {
            for (idim = 1; idim <= nd_v; idim++) {
                potsort[FA2(idim, i, nd_v)] = 0;
                gradsort[FA3(idim, 1, i, nd_v, 3)] = 0;
                gradsort[FA3(idim, 2, i, nd_v, 3)] = 0;
                gradsort[FA3(idim, 3, i, nd_v, 3)] = 0;
                hesssort[FA3(idim, 1, i, nd_v, 3)] = 0;
                hesssort[FA3(idim, 2, i, nd_v, 3)] = 0;
                hesssort[FA3(idim, 3, i, nd_v, 3)] = 0;
            }
        }
    }

    if (ifpghtarg_v == 1) {
        for (i = 1; i <= nt_v; i++) {
            for (idim = 1; idim <= nd_v; idim++) {
                pottargsort[FA2(idim, i, nd_v)] = 0;
            }
        }
    }
    if (ifpghtarg_v == 2) {
        for (i = 1; i <= nt_v; i++) {
            for (idim = 1; idim <= nd_v; idim++) {
                pottargsort[FA2(idim, i, nd_v)] = 0;
                gradtargsort[FA3(idim, 1, i, nd_v, 3)] = 0;
                gradtargsort[FA3(idim, 2, i, nd_v, 3)] = 0;
                gradtargsort[FA3(idim, 3, i, nd_v, 3)] = 0;
            }
        }
    }
    if (ifpghtarg_v == 3) {
        for (i = 1; i <= nt_v; i++) {
            for (idim = 1; idim <= nd_v; idim++) {
                pottargsort[FA2(idim, i, nd_v)] = 0;
                gradtargsort[FA3(idim, 1, i, nd_v, 3)] = 0;
                gradtargsort[FA3(idim, 2, i, nd_v, 3)] = 0;
                gradtargsort[FA3(idim, 3, i, nd_v, 3)] = 0;
                hesstargsort[FA3(idim, 1, i, nd_v, 3)] = 0;
                hesstargsort[FA3(idim, 2, i, nd_v, 3)] = 0;
                hesstargsort[FA3(idim, 3, i, nd_v, 3)] = 0;
            }
        }
    }

    /* compute scaling factors and lengths of mp/local expansions */
    rscales = (double *)malloc((nlevels + 1) * sizeof(double));
    nterms = (fint *)malloc((nlevels + 1) * sizeof(fint));

    nmax = 0;
    for (i = 0; i <= nlevels; i++) {
        rscales[i] = boxsize[i];
        bh2dterms_(eps, &nterms[i], ier);
        /* nterms(i) = nterms(i) - no-op in source */
        if (nterms[i] > nmax) nmax = nterms[i];
    }

    /* allocate iaddr and temporary arrays */
    iaddr = (fint *)malloc(2 * nboxes * sizeof(fint));
    lmptmp = (nmax + 1) * nd_v * 5;
    mptemp = (fcomplex *)malloc(lmptmp * sizeof(fcomplex));

    /* reorder sources. Strides:
     *   sources(2,*)  -> dim = 2
     *   charge(nd,2,*) is complex => 4*nd doubles per source
     *   dip(nd,3,*) is complex    => 6*nd doubles per source
     */
    {
        fint two = 2;
        fint chg_dim = 4 * nd_v;
        fint dip_dim = 6 * nd_v;
        dreorderf_(&two, ns, sources, sourcesort, isrc);
        if (ifcharge_v == 1) {
            dreorderf_(&chg_dim, ns, (const double *)charge,
                       (double *)chargesort, isrc);
        }
        if (ifdipole_v == 1) {
            dreorderf_(&dip_dim, ns, (const double *)dip,
                       (double *)dipsort, isrc);
        }
        /* reorder targets */
        dreorderf_(&two, nt, targ, targsort, itarg);
    }

    /* allocate workspace for multipole / local expansions */
    {
        const fint *laddr_p = &itree[iptr[0] - 1];
        FNAME(bh2dmpalloc)(&nd_v, laddr_p, iaddr, &nlevels, &lmptot, nterms);
    }

    rmlexp = (double *)malloc(lmptot * sizeof(double));

    /* call the main fmm routine */
    {
        const fint *laddr_p = &itree[iptr[0] - 1];
        FNAME(bhfmm2dmain)(&nd_v, eps,
                           ns, sourcesort,
                           &ifcharge_v, chargesort,
                           &ifdipole_v, dipsort,
                           nt, targsort, &nexpc, expc,
                           iaddr, rmlexp, mptemp, &lmptmp,
                           itree, &ltree, iptr, &ndiv, &nlevels,
                           &nboxes, iper, boxsize, rscales, tcenters,
                           laddr_p,
                           isrcse, itargse, iexpcse, nterms, &ntj,
                           &ifpgh_v, potsort, gradsort, hesssort,
                           &ifpghtarg_v, pottargsort, gradtargsort,
                           hesstargsort, jexps, &scj);
    }

    /* resort the output arrays in input order. Strides:
     *   pot(nd,*)     is complex => 2*nd doubles per element
     *   grad(nd,3,*)  is complex => 6*nd doubles per element
     *   hess(nd,3,*)  is complex => 6*nd doubles per element
     */
    {
        fint twond = 2 * nd_v;
        fint sixnd = 6 * nd_v;
        if (ifpgh_v == 1) {
            dreorderi_(&twond, ns, (const double *)potsort,
                       (double *)pot, isrc);
        }
        if (ifpgh_v == 2) {
            dreorderi_(&twond, ns, (const double *)potsort,
                       (double *)pot, isrc);
            dreorderi_(&sixnd, ns, (const double *)gradsort,
                       (double *)grad, isrc);
        }
        if (ifpgh_v == 3) {
            dreorderi_(&twond, ns, (const double *)potsort,
                       (double *)pot, isrc);
            dreorderi_(&sixnd, ns, (const double *)gradsort,
                       (double *)grad, isrc);
            dreorderi_(&sixnd, ns, (const double *)hesssort,
                       (double *)hess, isrc);
        }

        if (ifpghtarg_v == 1) {
            dreorderi_(&twond, nt, (const double *)pottargsort,
                       (double *)pottarg, itarg);
        }
        if (ifpghtarg_v == 2) {
            dreorderi_(&twond, nt, (const double *)pottargsort,
                       (double *)pottarg, itarg);
            dreorderi_(&sixnd, nt, (const double *)gradtargsort,
                       (double *)gradtarg, itarg);
        }
        if (ifpghtarg_v == 3) {
            dreorderi_(&twond, nt, (const double *)pottargsort,
                       (double *)pottarg, itarg);
            dreorderi_(&sixnd, nt, (const double *)gradtargsort,
                       (double *)gradtarg, itarg);
            dreorderi_(&sixnd, nt, (const double *)hesstargsort,
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
    free(dipsort);
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
