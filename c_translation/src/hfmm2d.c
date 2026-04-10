/*
 * hfmm2d.c - C translation of src/helmholtz/hfmm2d.f
 *
 * Top-level user-facing 2D Helmholtz FMM driver and the
 * five subroutines from the Fortran source. Translated 1:1 from the
 * Fortran reference: same control flow, same allocations, same
 * operation order, same loop bounds. OpenMP directives are stripped
 * (sequential C). Print/log/error functions (prinf/prin2) are guarded
 * in the Fortran source by ifprint=0 and never executed; they are
 * omitted entirely. Timing calls (cpu_time, omp_get_wtime) are also
 * stripped: time1/time2/tt1/tt2 are left at zero so the timeinfo
 * outputs become deterministic zero values.
 *
 * Cross-file dependencies (hndiv2d, pts_tree_*, dreorder*,
 * h2dterms, computemnlists, computelists, h2dmpzero, h2dsigzero,
 * h2dformmp*, h2dformta*, h2dmpeval*, h2dtaeval*, h2dmpmp, h2dmploc,
 * h2dlocloc, h2dmpmphf, h2dloclochf, h2dmplochf, h2d_mptosig,
 * h2d_sig2exp, h2d_mkmpshift, h2d_mkm2ltrans, h2d_diagtrans,
 * h2d_direct*, next235, zffti) are issued via bare Fortran symbol
 * names so the diff test isolates this file. Same-file calls
 * (hfmm2d -> hfmm2dmain / h2dmpalloc, hfmm2dmain ->
 * hfmm2dexpc_direct / hfmm2dpart_direct) use the FNAME() dispatcher.
 */

#include "hfmm2d.h"

/* ================================================================ */
/* Cross-file calls: bare Fortran symbol names.                      */
/* ================================================================ */

/* hndiv2d (src/helmholtz/hndiv2d.f). */
extern void hndiv2d_(const double *eps, const fint *ns, const fint *nt,
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

/* h2dterms. */
extern void h2dterms_(const double *bsize, const fcomplex *zk,
                      const double *eps, fint *nterms, fint *ier);

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

/* h2dcommon. */
extern void h2dmpzero_(const fint *nd, fcomplex *mpole, const fint *nterms);
extern void h2dsigzero_(const fint *nd, fcomplex *sig, const fint *nsig);

/* helmrouts2d: form multipole. */
extern void h2dformmpc_(const fint *nd, const fcomplex *zk,
                        const double *rscale, const double *source,
                        const fint *ns, const fcomplex *charge,
                        const double *center, const fint *nterms,
                        fcomplex *mpole);
extern void h2dformmpd_(const fint *nd, const fcomplex *zk,
                        const double *rscale, const double *source,
                        const fint *ns, const fcomplex *dipstr,
                        const double *dipvec, const double *center,
                        const fint *nterms, fcomplex *mpole);
extern void h2dformmpcd_(const fint *nd, const fcomplex *zk,
                         const double *rscale, const double *source,
                         const fint *ns, const fcomplex *charge,
                         const fcomplex *dipstr, const double *dipvec,
                         const double *center, const fint *nterms,
                         fcomplex *mpole);

/* helmrouts2d: form local. */
extern void h2dformtac_(const fint *nd, const fcomplex *zk,
                        const double *rscale, const double *source,
                        const fint *ns, const fcomplex *charge,
                        const double *center, const fint *nterms,
                        fcomplex *local);
extern void h2dformtad_(const fint *nd, const fcomplex *zk,
                        const double *rscale, const double *source,
                        const fint *ns, const fcomplex *dipstr,
                        const double *dipvec, const double *center,
                        const fint *nterms, fcomplex *local);
extern void h2dformtacd_(const fint *nd, const fcomplex *zk,
                         const double *rscale, const double *source,
                         const fint *ns, const fcomplex *charge,
                         const fcomplex *dipstr, const double *dipvec,
                         const double *center, const fint *nterms,
                         fcomplex *local);

/* helmrouts2d: mp eval. */
extern void h2dmpevalp_(const fint *nd, const fcomplex *zk,
                        const double *rscale, const double *center,
                        const fcomplex *mpole, const fint *nterms,
                        const double *ztarg, const fint *ntarg,
                        fcomplex *pot);
extern void h2dmpevalg_(const fint *nd, const fcomplex *zk,
                        const double *rscale, const double *center,
                        const fcomplex *mpole, const fint *nterms,
                        const double *ztarg, const fint *ntarg,
                        fcomplex *pot, fcomplex *grad);
extern void h2dmpevalh_(const fint *nd, const fcomplex *zk,
                        const double *rscale, const double *center,
                        const fcomplex *mpole, const fint *nterms,
                        const double *ztarg, const fint *ntarg,
                        fcomplex *pot, fcomplex *grad, fcomplex *hess);

/* helmrouts2d: ta eval. */
extern void h2dtaevalp_(const fint *nd, const fcomplex *zk,
                        const double *rscale, const double *center,
                        const fcomplex *local, const fint *nterms,
                        const double *ztarg, const fint *ntarg,
                        fcomplex *pot);
extern void h2dtaevalg_(const fint *nd, const fcomplex *zk,
                        const double *rscale, const double *center,
                        const fcomplex *local, const fint *nterms,
                        const double *ztarg, const fint *ntarg,
                        fcomplex *pot, fcomplex *grad);
extern void h2dtaevalh_(const fint *nd, const fcomplex *zk,
                        const double *rscale, const double *center,
                        const fcomplex *local, const fint *nterms,
                        const double *ztarg, const fint *ntarg,
                        fcomplex *pot, fcomplex *grad, fcomplex *hess);

/* helmrouts2d: translation operators. */
extern void h2dmpmp_(const fint *nd, const fcomplex *zk,
                     const double *rscale1, const double *center1,
                     const fcomplex *hexp1, const fint *nterms1,
                     const double *rscale2, const double *center2,
                     fcomplex *hexp2, const fint *nterms2);

extern void h2dlocloc_(const fint *nd, const fcomplex *zk,
                       const double *rscale1, const double *center1,
                       const fcomplex *jexp1, const fint *nterms1,
                       const double *rscale2, const double *center2,
                       fcomplex *jexp2, const fint *nterms2);

extern void h2dmploc_(const fint *nd, const fcomplex *zk,
                      const double *rscale1, const double *center1,
                      const fcomplex *hexp, const fint *nterms1,
                      const double *rscale2, const double *center2,
                      fcomplex *jexp, const fint *nterms2);

/* wideband2d: high-frequency routines. */
extern void h2dmpmphf_(const fint *nd, const fcomplex *zk,
                       const double *rscale1, const double *center1,
                       const fcomplex *hexp1, const fint *nterms1,
                       const double *rscale2, const double *center2,
                       fcomplex *sig2, const fint *nterms2,
                       const fint *nsig, const fcomplex *wsave,
                       const fcomplex *transvec);

extern void h2dloclochf_(const fint *nd, const fcomplex *zk,
                         const double *rscale1, const double *center1,
                         const fcomplex *sig, const fint *nterms1,
                         const fint *nsig,
                         const double *rscale2, const double *center2,
                         fcomplex *hexp2, const fint *nterms2,
                         const fcomplex *transvec, const fcomplex *wsave);

extern void h2dmplochf_(const fint *nd, const fcomplex *zk,
                        const double *rscale1, const double *center1,
                        const fcomplex *sig, const fint *nterms1,
                        const double *rscale2, const double *center2,
                        fcomplex *sig2, const fint *nterms2,
                        const fint *nsig, const fcomplex *wsave,
                        const fcomplex *tvec);

extern void h2d_mptosig_(const fint *nd, const fint *nterms1,
                         const fint *nsig, const fcomplex *hexp,
                         fcomplex *sig, const fcomplex *wsave);

extern void h2d_sig2exp_(const fint *nd, const fint *nsig,
                         const fcomplex *sig, const fcomplex *wsave,
                         const fint *nterms, fcomplex *expans);

extern void h2d_mkmpshift_(const fcomplex *zk, const double *center1,
                           const fint *nterms1, const double *center2,
                           const fint *nterms2, const fint *nsig,
                           const fcomplex *wsave, fcomplex *transvec);

extern void h2d_mkm2ltrans_(const fcomplex *zk, const double *center1,
                            const fint *nterms1, const double *center2,
                            const fint *nterms2, const fint *nsig,
                            const fcomplex *wsave, fcomplex *transvec);

extern void h2d_diagtrans_(const fint *nd, const fint *nsig,
                           const fcomplex *sig, const fcomplex *transvec,
                           fcomplex *sig2);

/* helmkernels2d: direct evaluation. */
extern void h2d_directcp_(const fint *nd, const fcomplex *zk,
                          const double *sources, const fint *ns,
                          const fcomplex *charge,
                          const double *targ, const fint *nt,
                          fcomplex *pot, const double *thresh);
extern void h2d_directcg_(const fint *nd, const fcomplex *zk,
                          const double *sources, const fint *ns,
                          const fcomplex *charge,
                          const double *targ, const fint *nt,
                          fcomplex *pot, fcomplex *grad,
                          const double *thresh);
extern void h2d_directch_(const fint *nd, const fcomplex *zk,
                          const double *sources, const fint *ns,
                          const fcomplex *charge,
                          const double *targ, const fint *nt,
                          fcomplex *pot, fcomplex *grad, fcomplex *hess,
                          const double *thresh);

extern void h2d_directdp_(const fint *nd, const fcomplex *zk,
                          const double *sources, const fint *ns,
                          const fcomplex *dipstr, const double *dipvec,
                          const double *targ, const fint *nt,
                          fcomplex *pot, const double *thresh);
extern void h2d_directdg_(const fint *nd, const fcomplex *zk,
                          const double *sources, const fint *ns,
                          const fcomplex *dipstr, const double *dipvec,
                          const double *targ, const fint *nt,
                          fcomplex *pot, fcomplex *grad,
                          const double *thresh);
extern void h2d_directdh_(const fint *nd, const fcomplex *zk,
                          const double *sources, const fint *ns,
                          const fcomplex *dipstr, const double *dipvec,
                          const double *targ, const fint *nt,
                          fcomplex *pot, fcomplex *grad, fcomplex *hess,
                          const double *thresh);

extern void h2d_directcdp_(const fint *nd, const fcomplex *zk,
                           const double *sources, const fint *ns,
                           const fcomplex *charge, const fcomplex *dipstr,
                           const double *dipvec,
                           const double *targ, const fint *nt,
                           fcomplex *pot, const double *thresh);
extern void h2d_directcdg_(const fint *nd, const fcomplex *zk,
                           const double *sources, const fint *ns,
                           const fcomplex *charge, const fcomplex *dipstr,
                           const double *dipvec,
                           const double *targ, const fint *nt,
                           fcomplex *pot, fcomplex *grad,
                           const double *thresh);
extern void h2d_directcdh_(const fint *nd, const fcomplex *zk,
                           const double *sources, const fint *ns,
                           const fcomplex *charge, const fcomplex *dipstr,
                           const double *dipvec,
                           const double *targ, const fint *nt,
                           fcomplex *pot, fcomplex *grad, fcomplex *hess,
                           const double *thresh);

/* next235 (common/next235.f) - integer function, takes double. */
extern fint next235_(const double *base);

/* zffti - from FFTPACK (complex FFT init). */
extern void zffti_(const fint *n, fcomplex *wsave);

/* ================================================================ */
/* Forward declarations for same-file helpers.                       */
/* ================================================================ */

/* Forward declaration for same-file hfmm2dmain (not in header). */
void FNAME(hfmm2dmain)(const fint *nd, const double *eps,
    const fcomplex *zk, const fint *nsource, const double *sourcesort,
    const fint *ifcharge, const fcomplex *chargesort,
    const fint *ifdipole, const fcomplex *dipstrsort,
    const double *dipvecsort,
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
    const fint *isrcse, const fint *itargse, const fint *iexpcse,
    const fint *nterms, const fint *ntj,
    const fint *ifpgh, fcomplex *pot, fcomplex *grad, fcomplex *hess,
    const fint *ifpghtarg, fcomplex *pottarg,
    fcomplex *gradtarg, fcomplex *hesstarg,
    fcomplex *jsort, double *scjsort,
    const fint *ifnear, double *timeinfo, fint *ier);

/* Same-file helpers called via FNAME. */
void FNAME(hfmm2dexpc_direct)(const fint *nd,
    const fint *istart, const fint *iend,
    const fint *jstart, const fint *jend,
    const fcomplex *zk,
    const double *rscales, const fint *nlevels,
    const double *source,
    const fint *ifcharge, const fcomplex *charge,
    const fint *ifdipole, const fcomplex *dipstr,
    const double *dipvec,
    const double *targ, fcomplex *jexps,
    const double *scj, const fint *ntj);

void FNAME(hfmm2dpart_direct)(const fint *nd,
    const fint *istart, const fint *iend,
    const fint *jstart, const fint *jend,
    const fcomplex *zk,
    const double *source,
    const fint *ifcharge, const fcomplex *charge,
    const fint *ifdipole, const fcomplex *dipstr,
    const double *dipvec,
    const double *targ, const fint *ifpgh,
    fcomplex *pot, fcomplex *grad, fcomplex *hess,
    const double *thresh);


/*
 * laddr(k, ilev): k is 1-based, ilev is 0-based, leading dim 2.
 * Used for laddr(2, 0:nlevels) inside hfmm2dmain and h2dmpalloc.
 */
#define LADDR(k, ilev) ((ilev) * 2 + ((k) - 1))


/* ---------------------------------------------------------------- */
/* h2dmpalloc - lay out workspace for multipole, local, and diag    */
/*              form expansions                                     */
/* ---------------------------------------------------------------- */
void FNAME(h2dmpalloc)(const fint *nd, const fint *laddr, fint *iaddr,
                       const fint *nlevels, fint *lmptot,
                       const fint *nterms)
{
    fint nd_v = *nd;
    fint nlevels_v = *nlevels;
    fint istart, i, ibox, nn, itmp;
    double dn;

    /* multipole expansions: iaddr(1,ibox) */
    istart = 1;
    for (i = 0; i <= nlevels_v; i++) {
        nn = (2 * nterms[i] + 1) * 2 * nd_v;
        for (ibox = laddr[LADDR(1, i)]; ibox <= laddr[LADDR(2, i)]; ibox++) {
            itmp = ibox - laddr[LADDR(1, i)];
            iaddr[FA2(1, ibox, 4)] = istart + itmp * nn;
        }
        istart = istart + (laddr[LADDR(2, i)] - laddr[LADDR(1, i)] + 1) * nn;
    }

    /* local expansions: iaddr(2,ibox) */
    for (i = 0; i <= nlevels_v; i++) {
        nn = (2 * nterms[i] + 1) * 2 * nd_v;
        for (ibox = laddr[LADDR(1, i)]; ibox <= laddr[LADDR(2, i)]; ibox++) {
            itmp = ibox - laddr[LADDR(1, i)];
            iaddr[FA2(2, ibox, 4)] = istart + itmp * nn;
        }
        istart = istart + (laddr[LADDR(2, i)] - laddr[LADDR(1, i)] + 1) * nn;
    }

    /* outgoing diag form: iaddr(3,ibox) */
    for (i = 0; i <= nlevels_v; i++) {
        dn = 2 * (nterms[i] + nterms[i]) + 1;
        nn = 2 * nd_v * next235_(&dn);
        for (ibox = laddr[LADDR(1, i)]; ibox <= laddr[LADDR(2, i)]; ibox++) {
            itmp = ibox - laddr[LADDR(1, i)];
            iaddr[FA2(3, ibox, 4)] = istart + itmp * nn;
        }
        istart = istart + (laddr[LADDR(2, i)] - laddr[LADDR(1, i)] + 1) * nn;
    }

    /* incoming diag form: iaddr(4,ibox) */
    for (i = 0; i <= nlevels_v; i++) {
        dn = 2 * (nterms[i] + nterms[i]) + 1;
        nn = 2 * nd_v * next235_(&dn);
        for (ibox = laddr[LADDR(1, i)]; ibox <= laddr[LADDR(2, i)]; ibox++) {
            itmp = ibox - laddr[LADDR(1, i)];
            iaddr[FA2(4, ibox, 4)] = istart + itmp * nn;
        }
        istart = istart + (laddr[LADDR(2, i)] - laddr[LADDR(1, i)] + 1) * nn;
    }

    *lmptot = istart;
}


/* ---------------------------------------------------------------- */
/* hfmm2dexpc_direct - direct sources -> expansion-center contribs  */
/* ---------------------------------------------------------------- */
void FNAME(hfmm2dexpc_direct)(const fint *nd,
    const fint *istart, const fint *iend,
    const fint *jstart, const fint *jend,
    const fcomplex *zk,
    const double *rscales, const fint *nlevels,
    const double *source,
    const fint *ifcharge, const fcomplex *charge,
    const fint *ifdipole, const fcomplex *dipstr,
    const double *dipvec,
    const double *targ, fcomplex *jexps,
    const double *scj, const fint *ntj)
{
    fint nd_v = *nd;
    fint istart_v = *istart;
    fint iend_v = *iend;
    fint jstart_v = *jstart;
    fint jend_v = *jend;
    fint ifcharge_v = *ifcharge;
    fint ifdipole_v = *ifdipole;
    fint ntj_v = *ntj;
    fint ns, j;
    (void)rscales;
    (void)nlevels;

    ns = iend_v - istart_v + 1;
    for (j = jstart_v; j <= jend_v; j++) {
        /*
         * jexps(nd, -ntj:ntj, *): leading dim nd, second dim (2*ntj+1).
         * jexps(1, -ntj, j) is at offset (j - 1) * nd * (2*ntj + 1).
         * scj is 1D, scj(j) at offset (j - 1).
         * targ is (2, *), targ(1, j) at offset (j - 1) * 2.
         * source(1, istart) at offset (istart - 1) * 2.
         * charge(nd, *), charge(1, istart) at offset (istart - 1) * nd.
         * dipstr(nd, *), dipstr(1, istart) at offset (istart - 1) * nd.
         * dipvec(nd, 2, *), dipvec(1, 1, istart) at offset (istart - 1) * nd * 2.
         */
        const double *src_p = &source[(istart_v - 1) * 2];
        const fcomplex *chg_p = &charge[(istart_v - 1) * nd_v];
        const fcomplex *dip_p = &dipstr[(istart_v - 1) * nd_v];
        const double *dvec_p = &dipvec[(istart_v - 1) * nd_v * 2];
        const double *tg_p = &targ[(j - 1) * 2];
        fcomplex *jexp_p = &jexps[(j - 1) * nd_v * (2 * ntj_v + 1)];
        const double *scj_p = &scj[j - 1];

        /*
         * 1:1 translation of the Fortran call order in hfmm2dexpc_direct.
         * Note: the Fortran source passes arguments in the order shown below,
         * which matches the call as written in hfmm2d.f lines 1677-1692.
         */
        /* Note: the Fortran source has arguments charge/ns swapped in these
           calls, but the code is dead (nexpc=0 always). We use the correct
           argument order so the C compiles cleanly. */
        if (ifcharge_v == 1 && ifdipole_v == 0) {
            h2dformtac_(&nd_v, zk, scj_p,
                src_p, &ns, chg_p, tg_p,
                &ntj_v, jexp_p);
        }

        if (ifdipole_v == 1 && ifcharge_v == 0) {
            h2dformtad_(&nd_v, zk, scj_p,
                src_p, &ns, dip_p, dvec_p,
                tg_p, &ntj_v, jexp_p);
        }
        if (ifdipole_v == 1 && ifcharge_v == 1) {
            h2dformtacd_(&nd_v, zk, scj_p,
                src_p, &ns, chg_p, dip_p,
                dvec_p,
                tg_p, &ntj_v, jexp_p);
        }
    }
}


/* ---------------------------------------------------------------- */
/* hfmm2dpart_direct - direct particle-to-particle contributions    */
/* ---------------------------------------------------------------- */
void FNAME(hfmm2dpart_direct)(const fint *nd,
    const fint *istart, const fint *iend,
    const fint *jstart, const fint *jend,
    const fcomplex *zk,
    const double *source,
    const fint *ifcharge, const fcomplex *charge,
    const fint *ifdipole, const fcomplex *dipstr,
    const double *dipvec,
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
    fint ns, j;
    fint one = 1;

    const double *src_p = &source[(istart_v - 1) * 2];
    const fcomplex *chg_p = &charge[(istart_v - 1) * nd_v];
    const fcomplex *dip_p = &dipstr[(istart_v - 1) * nd_v];
    const double *dvec_p = &dipvec[(istart_v - 1) * nd_v * 2];

    ns = iend_v - istart_v + 1;

    if (ifcharge_v == 1 && ifdipole_v == 0) {
        if (ifpgh_v == 1) {
            for (j = jstart_v; j <= jend_v; j++) {
                h2d_directcp_(&nd_v, zk, src_p, &ns,
                    chg_p, &targ[(j - 1) * 2], &one,
                    &pot[(j - 1) * nd_v], thresh);
            }
        }
        if (ifpgh_v == 2) {
            for (j = jstart_v; j <= jend_v; j++) {
                h2d_directcg_(&nd_v, zk, src_p, &ns,
                    chg_p, &targ[(j - 1) * 2], &one,
                    &pot[(j - 1) * nd_v],
                    &grad[(j - 1) * nd_v * 2],
                    thresh);
            }
        }
        if (ifpgh_v == 3) {
            for (j = jstart_v; j <= jend_v; j++) {
                h2d_directch_(&nd_v, zk, src_p, &ns,
                    chg_p, &targ[(j - 1) * 2], &one,
                    &pot[(j - 1) * nd_v],
                    &grad[(j - 1) * nd_v * 2],
                    &hess[(j - 1) * nd_v * 3],
                    thresh);
            }
        }
    }

    if (ifcharge_v == 0 && ifdipole_v == 1) {
        if (ifpgh_v == 1) {
            for (j = jstart_v; j <= jend_v; j++) {
                h2d_directdp_(&nd_v, zk, src_p, &ns,
                    dip_p, dvec_p,
                    &targ[(j - 1) * 2], &one,
                    &pot[(j - 1) * nd_v], thresh);
            }
        }
        if (ifpgh_v == 2) {
            for (j = jstart_v; j <= jend_v; j++) {
                h2d_directdg_(&nd_v, zk, src_p, &ns,
                    dip_p, dvec_p,
                    &targ[(j - 1) * 2], &one,
                    &pot[(j - 1) * nd_v],
                    &grad[(j - 1) * nd_v * 2],
                    thresh);
            }
        }
        if (ifpgh_v == 3) {
            for (j = jstart_v; j <= jend_v; j++) {
                h2d_directdh_(&nd_v, zk, src_p, &ns,
                    dip_p, dvec_p,
                    &targ[(j - 1) * 2], &one,
                    &pot[(j - 1) * nd_v],
                    &grad[(j - 1) * nd_v * 2],
                    &hess[(j - 1) * nd_v * 3],
                    thresh);
            }
        }
    }

    if (ifcharge_v == 1 && ifdipole_v == 1) {
        if (ifpgh_v == 1) {
            for (j = jstart_v; j <= jend_v; j++) {
                h2d_directcdp_(&nd_v, zk, src_p, &ns,
                    chg_p, dip_p, dvec_p,
                    &targ[(j - 1) * 2], &one,
                    &pot[(j - 1) * nd_v], thresh);
            }
        }
        if (ifpgh_v == 2) {
            for (j = jstart_v; j <= jend_v; j++) {
                h2d_directcdg_(&nd_v, zk, src_p, &ns,
                    chg_p, dip_p, dvec_p,
                    &targ[(j - 1) * 2], &one,
                    &pot[(j - 1) * nd_v],
                    &grad[(j - 1) * nd_v * 2],
                    thresh);
            }
        }
        if (ifpgh_v == 3) {
            for (j = jstart_v; j <= jend_v; j++) {
                h2d_directcdh_(&nd_v, zk, src_p, &ns,
                    chg_p, dip_p, dvec_p,
                    &targ[(j - 1) * 2], &one,
                    &pot[(j - 1) * nd_v],
                    &grad[(j - 1) * nd_v * 2],
                    &hess[(j - 1) * nd_v * 3],
                    thresh);
            }
        }
    }
}


/* ---------------------------------------------------------------- */
/* hfmm2dmain - the main FMM engine (8 steps)                       */
/* ---------------------------------------------------------------- */
void FNAME(hfmm2dmain)(const fint *nd, const double *eps,
    const fcomplex *zk, const fint *nsource, const double *sourcesort,
    const fint *ifcharge, const fcomplex *chargesort,
    const fint *ifdipole, const fcomplex *dipstrsort,
    const double *dipvecsort,
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
    const fint *isrcse, const fint *itargse, const fint *iexpcse,
    const fint *nterms, const fint *ntj,
    const fint *ifpgh, fcomplex *pot, fcomplex *grad, fcomplex *hess,
    const fint *ifpghtarg, fcomplex *pottarg,
    fcomplex *gradtarg, fcomplex *hesstarg,
    fcomplex *jsort, double *scjsort,
    const fint *ifnear, double *timeinfo, fint *ier)
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
    fint ifnear_v = *ifnear;

    fint i, j, k, ilev, idim;
    fint ibox, jbox, nchild, npts;
    fint istart, iend;
    fint istarts, iends, istartt, iendt, istarte, iende;
    fint jstart, jend;
    fint ix, iy;
    fint nn;

    fint mnlist1 = 0, mnlist2 = 0, mnlist3 = 0, mnlist4 = 0;

    fint *list1 = NULL, *list2 = NULL, *list3 = NULL, *list4 = NULL;
    fint *nlist1s = NULL, *nlist2s = NULL, *nlist3s = NULL, *nlist4s = NULL;

    double timelev[201];
    double thresh;
    double zkiupbound, zi;
    double pi;
    double dn, dx, dy;
    double dlam, boxlam;
    double c1[2], c2[2];

    fint ilevhf;
    fint ni, nsig;
    fcomplex *wsave = NULL;
    fcomplex *transvecall = NULL;
    fcomplex *transvecmpmp = NULL;
    fcomplex *sig = NULL;

    /* Suppress unused-warnings on parameters that are passed through. */
    (void)eps;
    (void)nsource;
    (void)ntarget;
    (void)mptemp;
    (void)lmptmp;
    (void)ltree;
    (void)ndiv;

    pi = 4 * atan(1.0);

    /* upper limit for zk along imaginary axis */
    zkiupbound = 40.0;
    zi = cimag(*zk);

    *ier = 0;

    /* timelev(0:nlevels) zeroed. */
    for (i = 0; i <= nlevels_v; i++) {
        timelev[i] = 0.0;
    }

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

    /* set the expansion coefficients (jsort) to zero. */
    /* jsort(nd, -ntj:ntj, *): for i=1..nexpc, j=-ntj..ntj, idim=1..nd. */
    for (i = 1; i <= nexpc_v; i++) {
        for (j = -ntj_v; j <= ntj_v; j++) {
            for (idim = 1; idim <= nd_v; idim++) {
                /* offset = ((i-1)*(2*ntj+1) + (j+ntj))*nd + (idim-1) */
                jsort[((i - 1) * (2 * ntj_v + 1) + (j + ntj_v)) * nd_v + (idim - 1)] = 0.0;
            }
        }
    }

    for (i = 0; i < 8; i++) {
        timeinfo[i] = 0.0;
    }

    /* set all multipole and local expansions to zero */
    for (ilev = 0; ilev <= nlevels_v; ilev++) {
        for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            h2dmpzero_(&nd_v,
                       (fcomplex *)&rmlexp[iaddr[FA2(1, ibox, 4)] - 1],
                       &nterms[ilev]);
            h2dmpzero_(&nd_v,
                       (fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 4)] - 1],
                       &nterms[ilev]);
            dn = 2 * (nterms[ilev] + nterms[ilev]) + 1;
            nn = next235_(&dn);
            h2dsigzero_(&nd_v,
                        (fcomplex *)&rmlexp[iaddr[FA2(3, ibox, 4)] - 1],
                        &nn);
            h2dsigzero_(&nd_v,
                        (fcomplex *)&rmlexp[iaddr[FA2(4, ibox, 4)] - 1],
                        &nn);
        }
    }

    /* Set scjsort */
    for (ilev = 0; ilev <= nlevels_v; ilev++) {
        for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            nchild = itree[iptr[3] + ibox - 1 - 1]; /* itree(iptr(4)+ibox-1) */
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
      if (zi * boxsize[ilev] < zkiupbound) {

        if (ifcharge_v == 1 && ifdipole_v == 0) {
            for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
                nchild = itree[iptr[3] + ibox - 1 - 1];
                istart = isrcse[FA2(1, ibox, 2)];
                iend = isrcse[FA2(2, ibox, 2)];
                npts = iend - istart + 1;
                if (nchild == 0 && npts > 0) {
                    h2dformmpc_(&nd_v, zk, &rscales[ilev],
                                &sourcesort[(istart - 1) * 2], &npts,
                                &chargesort[(istart - 1) * nd_v],
                                &centers[(ibox - 1) * 2], &nterms[ilev],
                                (fcomplex *)&rmlexp[iaddr[FA2(1, ibox, 4)] - 1]);
                }
            }
        }

        if (ifdipole_v == 1 && ifcharge_v == 0) {
            for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
                nchild = itree[iptr[3] + ibox - 1 - 1];
                istart = isrcse[FA2(1, ibox, 2)];
                iend = isrcse[FA2(2, ibox, 2)];
                npts = iend - istart + 1;
                if (nchild == 0 && npts > 0) {
                    h2dformmpd_(&nd_v, zk, &rscales[ilev],
                                &sourcesort[(istart - 1) * 2], &npts,
                                &dipstrsort[(istart - 1) * nd_v],
                                &dipvecsort[(istart - 1) * nd_v * 2],
                                &centers[(ibox - 1) * 2], &nterms[ilev],
                                (fcomplex *)&rmlexp[iaddr[FA2(1, ibox, 4)] - 1]);
                }
            }
        }

        if (ifdipole_v == 1 && ifcharge_v == 1) {
            for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
                nchild = itree[iptr[3] + ibox - 1 - 1];
                istart = isrcse[FA2(1, ibox, 2)];
                iend = isrcse[FA2(2, ibox, 2)];
                npts = iend - istart + 1;
                if (nchild == 0 && npts > 0) {
                    h2dformmpcd_(&nd_v, zk, &rscales[ilev],
                                 &sourcesort[(istart - 1) * 2], &npts,
                                 &chargesort[(istart - 1) * nd_v],
                                 &dipstrsort[(istart - 1) * nd_v],
                                 &dipvecsort[(istart - 1) * nd_v * 2],
                                 &centers[(ibox - 1) * 2], &nterms[ilev],
                                 (fcomplex *)&rmlexp[iaddr[FA2(1, ibox, 4)] - 1]);
                }
            }
        }
      } /* zi*boxsize check */
    }
    timeinfo[0] = 0.0;

    /* ============================================================ */
    /* Step 2: form lo (form local expansions from list 4)          */
    /* ============================================================ */
    for (ilev = 2; ilev <= nlevels_v; ilev++) {
      if (zi * boxsize[ilev] < zkiupbound) {
        if (ifcharge_v == 1 && ifdipole_v == 0) {
            for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
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
                        h2dformtac_(&nd_v, zk, &rscales[ilev],
                                    &sourcesort[(istart - 1) * 2], &npts,
                                    &chargesort[(istart - 1) * nd_v],
                                    &centers[(ibox - 1) * 2], &nterms[ilev],
                                    (fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 4)] - 1]);
                    }
                }
            }
        }

        if (ifcharge_v == 0 && ifdipole_v == 1) {
            for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
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
                        h2dformtad_(&nd_v, zk, &rscales[ilev],
                                    &sourcesort[(istart - 1) * 2], &npts,
                                    &dipstrsort[(istart - 1) * nd_v],
                                    &dipvecsort[(istart - 1) * nd_v * 2],
                                    &centers[(ibox - 1) * 2], &nterms[ilev],
                                    (fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 4)] - 1]);
                    }
                }
            }
        }

        if (ifcharge_v == 1 && ifdipole_v == 1) {
            for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
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
                        h2dformtacd_(&nd_v, zk, &rscales[ilev],
                                     &sourcesort[(istart - 1) * 2], &npts,
                                     &chargesort[(istart - 1) * nd_v],
                                     &dipstrsort[(istart - 1) * nd_v],
                                     &dipvecsort[(istart - 1) * nd_v * 2],
                                     &centers[(ibox - 1) * 2], &nterms[ilev],
                                     (fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 4)] - 1]);
                    }
                }
            }
        }
      } /* zi*boxsize check */
    }
    timeinfo[1] = 0.0;

    /* ============================================================ */
    /* Identify level ilevhf where HF regime begins                 */
    /* ============================================================ */
    ilevhf = 0;
    for (ilev = nlevels_v - 1; ilev >= 1; ilev--) {
        if (zi * boxsize[ilev] < zkiupbound) {
            dlam = creal(*zk);
            /* dlam = zk; dlam = 1/(dlam/(2*pi)) */
            /* In Fortran: dlam = zk assigns real part; dlam = 1/(dlam/(2*pi)) */
            dlam = 1.0 / (dlam / (2 * pi));
            boxlam = boxsize[ilev] / dlam;
            if (boxlam > 16.0) {
                ilevhf = ilev;
                goto label_111;
            }
        }
    }
label_111: ;

    /* ============================================================ */
    /* Step 3a: merge mp low freq (M2M up the tree, low freq part)  */
    /* ============================================================ */
    for (ilev = nlevels_v - 1; ilev >= ilevhf + 1; ilev--) {
      if (zi * boxsize[ilev] < zkiupbound) {
        for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            nchild = itree[iptr[3] + ibox - 1 - 1];
            for (i = 1; i <= nchild; i++) {
                jbox = itree[iptr[4] + 4 * (ibox - 1) + i - 1 - 1];
                istart = isrcse[FA2(1, jbox, 2)];
                iend = isrcse[FA2(2, jbox, 2)];
                npts = iend - istart + 1;
                if (npts > 0) {
                    fint nterms_p1 = nterms[ilev + 1];
                    fint nterms_l = nterms[ilev];
                    h2dmpmp_(&nd_v, zk, &rscales[ilev + 1],
                             &centers[(jbox - 1) * 2],
                             (fcomplex *)&rmlexp[iaddr[FA2(1, jbox, 4)] - 1],
                             &nterms_p1, &rscales[ilev],
                             &centers[(ibox - 1) * 2],
                             (fcomplex *)&rmlexp[iaddr[FA2(1, ibox, 4)] - 1],
                             &nterms_l);
                }
            }
        }
      }
    }

    /* ============================================================ */
    /* Step 3b: merge mp high freq (M2M up the tree, HF part)       */
    /* ============================================================ */
    for (ilev = ilevhf; ilev >= 1; ilev--) {
      if (zi * boxsize[ilev] < zkiupbound) {
        dn = 2 * (nterms[ilev] + nterms[ilev + 1]) + 1;
        nsig = next235_(&dn);
        wsave = (fcomplex *)malloc((4 * nsig + 100) * sizeof(fcomplex));
        transvecmpmp = (fcomplex *)malloc(nsig * 4 * sizeof(fcomplex));
        zffti_(&nsig, wsave);
        c2[0] = 0.0;
        c2[1] = 0.0;
        for (jbox = 1; jbox <= 4; jbox++) {
            k = 2;
            if (jbox <= 2) k = 1;
            c1[0] = 0.25 * boxsize[ilev] * (jbox % 2 == 0 ? 1 : -1);
            c1[1] = 0.25 * boxsize[ilev] * (k % 2 == 0 ? 1 : -1);
            {
                fint nterms_p1 = nterms[ilev + 1];
                fint nterms_l = nterms[ilev];
                h2d_mkmpshift_(zk, c1, &nterms_p1,
                    c2, &nterms_l, &nsig, wsave,
                    &transvecmpmp[(jbox - 1) * nsig]);
            }
        }

        sig = (fcomplex *)malloc(nd_v * nsig * sizeof(fcomplex));

        for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            nchild = itree[iptr[3] + ibox - 1 - 1];
            h2dsigzero_(&nd_v, sig, &nsig);
            for (i = 1; i <= nchild; i++) {
                jbox = itree[iptr[4] + 4 * (ibox - 1) + i - 1 - 1];
                istart = isrcse[FA2(1, jbox, 2)];
                iend = isrcse[FA2(2, jbox, 2)];
                npts = iend - istart + 1;
                if (npts > 0) {
                    fint nterms_p1 = nterms[ilev + 1];
                    fint nterms_l = nterms[ilev];
                    h2dmpmphf_(&nd_v, zk, &rscales[ilev + 1],
                        &centers[(jbox - 1) * 2],
                        (fcomplex *)&rmlexp[iaddr[FA2(1, jbox, 4)] - 1],
                        &nterms_p1, &rscales[ilev],
                        &centers[(ibox - 1) * 2],
                        sig, &nterms_l, &nsig, wsave,
                        &transvecmpmp[(i - 1) * nsig]);
                }
            }
            {
                fint nterms_l = nterms[ilev];
                h2d_sig2exp_(&nd_v, &nsig, sig, wsave, &nterms_l,
                    (fcomplex *)&rmlexp[iaddr[FA2(1, ibox, 4)] - 1]);
            }
        }
        free(wsave);
        free(transvecmpmp);
        free(sig);
        wsave = NULL;
        transvecmpmp = NULL;
        sig = NULL;
      }
    }
    timeinfo[2] = 0.0;

    /* ============================================================ */
    /* Step 4: mp to loc (list 2 / M2L)                             */
    /* ============================================================ */
    for (ilev = 2; ilev <= nlevels_v; ilev++) {
        ni = nterms[ilev];
        dn = 2 * (ni + ni) + 1;
        nsig = next235_(&dn);
        wsave = (fcomplex *)malloc((4 * nsig + 100) * sizeof(fcomplex));
        /* transvecall(nsig, -3:3, -3:3) = nsig * 7 * 7 */
        transvecall = (fcomplex *)malloc(nsig * 7 * 7 * sizeof(fcomplex));
        dlam = creal(*zk);
        dlam = 1.0 / (dlam / (2 * pi));
        boxlam = boxsize[ilev] / dlam;
        zffti_(&nsig, wsave);

        /* precompute mp to sig for all boxes */
        if (boxlam > 16.0) {
            for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
                fint nterms_l = nterms[ilev];
                h2d_mptosig_(&nd_v, &nterms_l, &nsig,
                    (fcomplex *)&rmlexp[iaddr[FA2(1, ibox, 4)] - 1],
                    (fcomplex *)&rmlexp[iaddr[FA2(3, ibox, 4)] - 1],
                    wsave);
            }

            /* evaluate diagonal shifts for translation operators */
            c1[0] = 0.0;
            c1[1] = 0.0;
            for (ix = -3; ix <= 3; ix++) {
                for (iy = -3; iy <= 3; iy++) {
                    fint nterms_l = nterms[ilev];
                    c2[0] = ix * boxsize[ilev];
                    c2[1] = iy * boxsize[ilev];
                    /* transvecall(1, ix, iy): offset = ((iy+3)*7 + (ix+3)) * nsig */
                    h2d_mkm2ltrans_(zk, c1, &nterms_l,
                        c2, &nterms_l, &nsig, wsave,
                        &transvecall[((iy + 3) * 7 + (ix + 3)) * nsig]);
                }
            }
        }

      if (zi * boxsize[ilev] < zkiupbound) {
        for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
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
                    jbox = list2[FA2(i, ibox, mnlist2)];

                    if (boxlam > 16.0) {
                        dx = centers[(ibox - 1) * 2] - centers[(jbox - 1) * 2];
                        ix = (fint)(dx / boxsize[ilev] + (dx / boxsize[ilev] >= 0 ? 0.5 : -0.5));
                        dy = centers[(ibox - 1) * 2 + 1] - centers[(jbox - 1) * 2 + 1];
                        iy = (fint)(dy / boxsize[ilev] + (dy / boxsize[ilev] >= 0 ? 0.5 : -0.5));
                        h2d_diagtrans_(&nd_v, &nsig,
                            (fcomplex *)&rmlexp[iaddr[FA2(3, jbox, 4)] - 1],
                            &transvecall[((iy + 3) * 7 + (ix + 3)) * nsig],
                            (fcomplex *)&rmlexp[iaddr[FA2(4, ibox, 4)] - 1]);
                    } else {
                        fint nterms_l = nterms[ilev];
                        h2dmploc_(&nd_v, zk, &rscales[ilev],
                                  &centers[(jbox - 1) * 2],
                                  (fcomplex *)&rmlexp[iaddr[FA2(1, jbox, 4)] - 1],
                                  &nterms_l, &rscales[ilev],
                                  &centers[(ibox - 1) * 2],
                                  (fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 4)] - 1],
                                  &nterms_l);
                    }
                }
            }
        }

        /* convert diagonal form back to local expansion if HF */
        if (boxlam > 16.0) {
            for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
                fint nterms_l = nterms[ilev];
                h2d_sig2exp_(&nd_v, &nsig,
                    (fcomplex *)&rmlexp[iaddr[FA2(4, ibox, 4)] - 1],
                    wsave, &nterms_l,
                    (fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 4)] - 1]);
            }
        }
      }

        timelev[ilev] = 0.0;
        free(wsave);
        free(transvecall);
        wsave = NULL;
        transvecall = NULL;
    }
    timeinfo[3] = 0.0;

    /* ============================================================ */
    /* Step 5: split loc (L2L down the tree)                        */
    /* ============================================================ */
    for (ilev = 1; ilev <= nlevels_v - 1; ilev++) {
        dn = 2 * (nterms[ilev] + nterms[ilev + 1]) + 1;
        nsig = next235_(&dn);
        wsave = (fcomplex *)malloc((4 * nsig + 100) * sizeof(fcomplex));
        sig = (fcomplex *)malloc(nd_v * nsig * sizeof(fcomplex));
        transvecmpmp = (fcomplex *)malloc(nsig * 4 * sizeof(fcomplex));
        zffti_(&nsig, wsave);
        c1[0] = 0.0;
        c1[1] = 0.0;
        for (jbox = 1; jbox <= 4; jbox++) {
            k = 2;
            if (jbox <= 2) k = 1;
            c2[0] = 0.25 * boxsize[ilev] * (jbox % 2 == 0 ? 1 : -1);
            c2[1] = 0.25 * boxsize[ilev] * (k % 2 == 0 ? 1 : -1);
            {
                fint nterms_p1 = nterms[ilev + 1];
                fint nterms_l = nterms[ilev];
                h2d_mkmpshift_(zk, c1, &nterms_p1,
                    c2, &nterms_l, &nsig, wsave,
                    &transvecmpmp[(jbox - 1) * nsig]);
            }
        }
        dlam = creal(*zk);
        dlam = 1.0 / (dlam / (2 * pi));
        boxlam = boxsize[ilev] / dlam;

        /* convert local expansion to diagonal form */
        for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            fint nterms_l = nterms[ilev];
            h2d_mptosig_(&nd_v, &nterms_l, &nsig,
                (fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 4)] - 1],
                (fcomplex *)&rmlexp[iaddr[FA2(4, ibox, 4)] - 1],
                wsave);
        }

      if (zi * boxsize[ilev] < zkiupbound) {
        for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
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
                    jbox = itree[iptr[4] + 4 * (ibox - 1) + i - 1 - 1];

                    if (boxlam > 16.0) {
                        fint nterms_l = nterms[ilev];
                        fint nterms_p1 = nterms[ilev + 1];
                        h2dloclochf_(&nd_v, zk, &rscales[ilev],
                            &centers[(ibox - 1) * 2],
                            (fcomplex *)&rmlexp[iaddr[FA2(4, ibox, 4)] - 1],
                            &nterms_l, &nsig, &rscales[ilev + 1],
                            &centers[(jbox - 1) * 2],
                            (fcomplex *)&rmlexp[iaddr[FA2(2, jbox, 4)] - 1],
                            &nterms_p1,
                            &transvecmpmp[(i - 1) * nsig], wsave);
                    } else {
                        fint nterms_l = nterms[ilev];
                        fint nterms_p1 = nterms[ilev + 1];
                        h2dlocloc_(&nd_v, zk, &rscales[ilev],
                            &centers[(ibox - 1) * 2],
                            (fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 4)] - 1],
                            &nterms_l, &rscales[ilev + 1],
                            &centers[(jbox - 1) * 2],
                            (fcomplex *)&rmlexp[iaddr[FA2(2, jbox, 4)] - 1],
                            &nterms_p1);
                    }
                }
            }
        }
      }
        free(wsave);
        free(transvecmpmp);
        free(sig);
        wsave = NULL;
        transvecmpmp = NULL;
        sig = NULL;
    }
    timeinfo[4] = 0.0;

    /* ============================================================ */
    /* Step 6: mp eval (list 3)                                     */
    /* ============================================================ */
    for (ilev = 1; ilev <= nlevels_v - 1; ilev++) {
      if (zi * boxsize[ilev + 1] < zkiupbound) {
        for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            istart = iexpcse[FA2(1, ibox, 2)];
            iend = iexpcse[FA2(2, ibox, 2)];
            for (j = istart; j <= iend; j++) {
                for (i = 1; i <= nlist3s[ibox - 1]; i++) {
                    fint nterms_p1 = nterms[ilev + 1];
                    jbox = list3[FA2(i, ibox, mnlist3)];
                    dlam = creal(*zk);
                    dlam = 1.0 / (dlam / (2 * pi));
                    boxlam = boxsize[ilev] / dlam;
                    h2dmploc_(&nd_v, zk, &rscales[ilev + 1],
                              &centers[(jbox - 1) * 2],
                              (fcomplex *)&rmlexp[iaddr[FA2(1, jbox, 4)] - 1],
                              &nterms_p1, &scjsort[j - 1],
                              &expcsort[(j - 1) * 2],
                              &jsort[(j - 1) * nd_v * (2 * ntj_v + 1)],
                              &ntj_v);
                }
            }

            /* evaluate multipole expansion at all targets */
            istart = itargse[FA2(1, ibox, 2)];
            iend = itargse[FA2(2, ibox, 2)];
            npts = iend - istart + 1;

            if (ifpghtarg_v == 1) {
                for (i = 1; i <= nlist3s[ibox - 1]; i++) {
                    fint nterms_p1 = nterms[ilev + 1];
                    jbox = list3[FA2(i, ibox, mnlist3)];
                    h2dmpevalp_(&nd_v, zk, &rscales[ilev + 1],
                                &centers[(jbox - 1) * 2],
                                (const fcomplex *)&rmlexp[iaddr[FA2(1, jbox, 4)] - 1],
                                &nterms_p1, &targetsort[(istart - 1) * 2], &npts,
                                &pottarg[(istart - 1) * nd_v]);
                }
            }
            if (ifpghtarg_v == 2) {
                for (i = 1; i <= nlist3s[ibox - 1]; i++) {
                    fint nterms_p1 = nterms[ilev + 1];
                    jbox = list3[FA2(i, ibox, mnlist3)];
                    h2dmpevalg_(&nd_v, zk, &rscales[ilev + 1],
                                &centers[(jbox - 1) * 2],
                                (const fcomplex *)&rmlexp[iaddr[FA2(1, jbox, 4)] - 1],
                                &nterms_p1, &targetsort[(istart - 1) * 2], &npts,
                                &pottarg[(istart - 1) * nd_v],
                                &gradtarg[(istart - 1) * nd_v * 2]);
                }
            }
            if (ifpghtarg_v == 3) {
                for (i = 1; i <= nlist3s[ibox - 1]; i++) {
                    fint nterms_p1 = nterms[ilev + 1];
                    jbox = list3[FA2(i, ibox, mnlist3)];
                    h2dmpevalh_(&nd_v, zk, &rscales[ilev + 1],
                                &centers[(jbox - 1) * 2],
                                (const fcomplex *)&rmlexp[iaddr[FA2(1, jbox, 4)] - 1],
                                &nterms_p1, &targetsort[(istart - 1) * 2], &npts,
                                &pottarg[(istart - 1) * nd_v],
                                &gradtarg[(istart - 1) * nd_v * 2],
                                &hesstarg[(istart - 1) * nd_v * 3]);
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
                    h2dmpevalp_(&nd_v, zk, &rscales[ilev + 1],
                                &centers[(jbox - 1) * 2],
                                (const fcomplex *)&rmlexp[iaddr[FA2(1, jbox, 4)] - 1],
                                &nterms_p1, &sourcesort[(istart - 1) * 2], &npts,
                                &pot[(istart - 1) * nd_v]);
                }
            }
            if (ifpgh_v == 2) {
                for (i = 1; i <= nlist3s[ibox - 1]; i++) {
                    fint nterms_p1 = nterms[ilev + 1];
                    jbox = list3[FA2(i, ibox, mnlist3)];
                    h2dmpevalg_(&nd_v, zk, &rscales[ilev + 1],
                                &centers[(jbox - 1) * 2],
                                (const fcomplex *)&rmlexp[iaddr[FA2(1, jbox, 4)] - 1],
                                &nterms_p1, &sourcesort[(istart - 1) * 2], &npts,
                                &pot[(istart - 1) * nd_v],
                                &grad[(istart - 1) * nd_v * 2]);
                }
            }
            if (ifpgh_v == 3) {
                for (i = 1; i <= nlist3s[ibox - 1]; i++) {
                    fint nterms_p1 = nterms[ilev + 1];
                    jbox = list3[FA2(i, ibox, mnlist3)];
                    h2dmpevalh_(&nd_v, zk, &rscales[ilev + 1],
                                &centers[(jbox - 1) * 2],
                                (const fcomplex *)&rmlexp[iaddr[FA2(1, jbox, 4)] - 1],
                                &nterms_p1, &sourcesort[(istart - 1) * 2], &npts,
                                &pot[(istart - 1) * nd_v],
                                &grad[(istart - 1) * nd_v * 2],
                                &hess[(istart - 1) * nd_v * 3]);
                }
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
      if (zi * boxsize[ilev] < zkiupbound) {
        for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            nchild = itree[iptr[3] + ibox - 1 - 1];
            if (nchild == 0) {
                istart = iexpcse[FA2(1, ibox, 2)];
                iend = iexpcse[FA2(2, ibox, 2)];
                for (i = istart; i <= iend; i++) {
                    fint nterms_l = nterms[ilev];
                    h2dlocloc_(&nd_v, zk, &rscales[ilev],
                               &centers[(ibox - 1) * 2],
                               (fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 4)] - 1],
                               &nterms_l, &scjsort[i - 1],
                               &expcsort[(i - 1) * 2],
                               &jsort[(i - 1) * nd_v * (2 * ntj_v + 1)],
                               &ntj_v);
                }

                /* evaluate local expansion at targets */
                istart = itargse[FA2(1, ibox, 2)];
                iend = itargse[FA2(2, ibox, 2)];
                npts = iend - istart + 1;
                if (ifpghtarg_v == 1) {
                    fint nterms_l = nterms[ilev];
                    h2dtaevalp_(&nd_v, zk, &rscales[ilev],
                                &centers[(ibox - 1) * 2],
                                (const fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 4)] - 1],
                                &nterms_l, &targetsort[(istart - 1) * 2], &npts,
                                &pottarg[(istart - 1) * nd_v]);
                }
                if (ifpghtarg_v == 2) {
                    fint nterms_l = nterms[ilev];
                    h2dtaevalg_(&nd_v, zk, &rscales[ilev],
                                &centers[(ibox - 1) * 2],
                                (const fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 4)] - 1],
                                &nterms_l, &targetsort[(istart - 1) * 2], &npts,
                                &pottarg[(istart - 1) * nd_v],
                                &gradtarg[(istart - 1) * nd_v * 2]);
                }
                if (ifpghtarg_v == 3) {
                    fint nterms_l = nterms[ilev];
                    h2dtaevalh_(&nd_v, zk, &rscales[ilev],
                                &centers[(ibox - 1) * 2],
                                (const fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 4)] - 1],
                                &nterms_l, &targetsort[(istart - 1) * 2], &npts,
                                &pottarg[(istart - 1) * nd_v],
                                &gradtarg[(istart - 1) * nd_v * 2],
                                &hesstarg[(istart - 1) * nd_v * 3]);
                }

                /* evaluate local expansion at sources */
                istart = isrcse[FA2(1, ibox, 2)];
                iend = isrcse[FA2(2, ibox, 2)];
                npts = iend - istart + 1;
                if (ifpgh_v == 1) {
                    fint nterms_l = nterms[ilev];
                    h2dtaevalp_(&nd_v, zk, &rscales[ilev],
                                &centers[(ibox - 1) * 2],
                                (const fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 4)] - 1],
                                &nterms_l, &sourcesort[(istart - 1) * 2], &npts,
                                &pot[(istart - 1) * nd_v]);
                }
                if (ifpgh_v == 2) {
                    fint nterms_l = nterms[ilev];
                    h2dtaevalg_(&nd_v, zk, &rscales[ilev],
                                &centers[(ibox - 1) * 2],
                                (const fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 4)] - 1],
                                &nterms_l, &sourcesort[(istart - 1) * 2], &npts,
                                &pot[(istart - 1) * nd_v],
                                &grad[(istart - 1) * nd_v * 2]);
                }
                if (ifpgh_v == 3) {
                    fint nterms_l = nterms[ilev];
                    h2dtaevalh_(&nd_v, zk, &rscales[ilev],
                                &centers[(ibox - 1) * 2],
                                (const fcomplex *)&rmlexp[iaddr[FA2(2, ibox, 4)] - 1],
                                &nterms_l, &sourcesort[(istart - 1) * 2], &npts,
                                &pot[(istart - 1) * nd_v],
                                &grad[(istart - 1) * nd_v * 2],
                                &hess[(istart - 1) * nd_v * 3]);
                }
            }
        }
      }
    }

    timeinfo[6] = 0.0;

    /* ============================================================ */
    /* Step 8: direct                                               */
    /* ============================================================ */
    /* thresh = abs(zk)*boxsize(0)*2.0d0**(-51) */
    thresh = cabs(*zk) * boxsize[0] * 0x1p-51;

    if (ifnear_v == 0) goto skip_near;

    for (ilev = 0; ilev <= nlevels_v; ilev++) {
        for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {

            istartt = itargse[FA2(1, ibox, 2)];
            iendt = itargse[FA2(2, ibox, 2)];

            istarte = iexpcse[FA2(1, ibox, 2)];
            iende = iexpcse[FA2(2, ibox, 2)];

            istarts = isrcse[FA2(1, ibox, 2)];
            iends = isrcse[FA2(2, ibox, 2)];

            for (i = 1; i <= nlist1s[ibox - 1]; i++) {
                jbox = list1[FA2(i, ibox, mnlist1)];

                jstart = isrcse[FA2(1, jbox, 2)];
                jend = isrcse[FA2(2, jbox, 2)];

                FNAME(hfmm2dexpc_direct)(&nd_v, &jstart, &jend, &istarte,
                                         &iende, zk, rscales, &nlevels_v,
                                         sourcesort, &ifcharge_v, chargesort,
                                         &ifdipole_v, dipstrsort,
                                         dipvecsort,
                                         expcsort, jsort, scjsort, &ntj_v);

                FNAME(hfmm2dpart_direct)(&nd_v, &jstart, &jend, &istartt,
                                         &iendt, zk,
                                         sourcesort, &ifcharge_v,
                                         chargesort, &ifdipole_v, dipstrsort,
                                         dipvecsort,
                                         targetsort, &ifpghtarg_v, pottarg,
                                         gradtarg, hesstarg, &thresh);

                FNAME(hfmm2dpart_direct)(&nd_v, &jstart, &jend, &istarts,
                                         &iends, zk,
                                         sourcesort, &ifcharge_v,
                                         chargesort, &ifdipole_v, dipstrsort,
                                         dipvecsort,
                                         sourcesort, &ifpgh_v, pot, grad,
                                         hess, &thresh);
            }
        }
    }
skip_near:;
    timeinfo[7] = 0.0;

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
/* hfmm2d - top-level user-callable driver                          */
/* ---------------------------------------------------------------- */
void FNAME(hfmm2d)(const fint *nd, const double *eps, const fcomplex *zk,
    const fint *ns, const double *sources,
    const fint *ifcharge, const fcomplex *charge,
    const fint *ifdipole, const fcomplex *dipstr, const double *dipvec,
    fint *iper, const fint *ifpgh, fcomplex *pot, fcomplex *grad,
    fcomplex *hess, const fint *nt, const double *targ,
    const fint *ifpghtarg, fcomplex *pottarg, fcomplex *gradtarg,
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
    fint ifnear;
    double timeinfo[8];

    fint *isrc = NULL, *isrcse = NULL;
    fint *itarg = NULL, *itargse = NULL, *iexpcse = NULL;
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
    fint iert;
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

    hndiv2d_(eps, ns, nt, ifcharge, ifdipole, ifpgh, ifpghtarg,
             &ndiv, &idivflag);

    ltree = 0;
    nlmin = 0;
    nlmax = 51;
    ifunif = 0;
    *iper = 0;

    ifnear = 1;

    /* initialize timeinfo */
    for (i = 0; i < 8; i++) {
        timeinfo[i] = 0.0;
    }

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

    if (ifcharge_v == 1 && ifdipole_v == 0) {
        chargesort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        dipstrsort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
        dipvecsort = (double *)malloc(nd_v * 2 * 1 * sizeof(double));
    }
    if (ifcharge_v == 0 && ifdipole_v == 1) {
        chargesort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
        dipstrsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        dipvecsort = (double *)malloc(nd_v * 2 * ns_v * sizeof(double));
    }
    if (ifcharge_v == 1 && ifdipole_v == 1) {
        chargesort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        dipstrsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        dipvecsort = (double *)malloc(nd_v * 2 * ns_v * sizeof(double));
    }
    /* If both ifcharge=0 and ifdipole=0, allocate placeholders. */
    if (chargesort == NULL) {
        chargesort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
    }
    if (dipstrsort == NULL) {
        dipstrsort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
    }
    if (dipvecsort == NULL) {
        dipvecsort = (double *)malloc(nd_v * 2 * 1 * sizeof(double));
    }

    if (ifpgh_v == 1) {
        potsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        gradsort = (fcomplex *)malloc(nd_v * 2 * 1 * sizeof(fcomplex));
        hesssort = (fcomplex *)malloc(nd_v * 3 * 1 * sizeof(fcomplex));
    } else if (ifpgh_v == 2) {
        potsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        gradsort = (fcomplex *)malloc(nd_v * 2 * ns_v * sizeof(fcomplex));
        hesssort = (fcomplex *)malloc(nd_v * 3 * 1 * sizeof(fcomplex));
    } else if (ifpgh_v == 3) {
        potsort = (fcomplex *)malloc(nd_v * ns_v * sizeof(fcomplex));
        gradsort = (fcomplex *)malloc(nd_v * 2 * ns_v * sizeof(fcomplex));
        hesssort = (fcomplex *)malloc(nd_v * 3 * ns_v * sizeof(fcomplex));
    } else {
        potsort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
        gradsort = (fcomplex *)malloc(nd_v * 2 * 1 * sizeof(fcomplex));
        hesssort = (fcomplex *)malloc(nd_v * 3 * 1 * sizeof(fcomplex));
    }

    if (ifpghtarg_v == 1) {
        pottargsort = (fcomplex *)malloc(nd_v * nt_v * sizeof(fcomplex));
        gradtargsort = (fcomplex *)malloc(nd_v * 2 * 1 * sizeof(fcomplex));
        hesstargsort = (fcomplex *)malloc(nd_v * 3 * 1 * sizeof(fcomplex));
    } else if (ifpghtarg_v == 2) {
        pottargsort = (fcomplex *)malloc(nd_v * nt_v * sizeof(fcomplex));
        gradtargsort = (fcomplex *)malloc(nd_v * 2 * nt_v * sizeof(fcomplex));
        hesstargsort = (fcomplex *)malloc(nd_v * 3 * 1 * sizeof(fcomplex));
    } else if (ifpghtarg_v == 3) {
        pottargsort = (fcomplex *)malloc(nd_v * nt_v * sizeof(fcomplex));
        gradtargsort = (fcomplex *)malloc(nd_v * 2 * nt_v * sizeof(fcomplex));
        hesstargsort = (fcomplex *)malloc(nd_v * 3 * nt_v * sizeof(fcomplex));
    } else {
        pottargsort = (fcomplex *)malloc(nd_v * 1 * sizeof(fcomplex));
        gradtargsort = (fcomplex *)malloc(nd_v * 2 * 1 * sizeof(fcomplex));
        hesstargsort = (fcomplex *)malloc(nd_v * 3 * 1 * sizeof(fcomplex));
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
                gradsort[FA3(idim, 1, i, nd_v, 2)] = 0.0;
                gradsort[FA3(idim, 2, i, nd_v, 2)] = 0.0;
            }
        }
    }
    if (ifpgh_v == 3) {
        for (i = 1; i <= ns_v; i++) {
            for (idim = 1; idim <= nd_v; idim++) {
                potsort[FA2(idim, i, nd_v)] = 0.0;
                gradsort[FA3(idim, 1, i, nd_v, 2)] = 0.0;
                gradsort[FA3(idim, 2, i, nd_v, 2)] = 0.0;
                hesssort[FA3(idim, 1, i, nd_v, 3)] = 0.0;
                hesssort[FA3(idim, 2, i, nd_v, 3)] = 0.0;
                hesssort[FA3(idim, 3, i, nd_v, 3)] = 0.0;
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
                gradtargsort[FA3(idim, 1, i, nd_v, 2)] = 0.0;
                gradtargsort[FA3(idim, 2, i, nd_v, 2)] = 0.0;
            }
        }
    }
    if (ifpghtarg_v == 3) {
        for (i = 1; i <= nt_v; i++) {
            for (idim = 1; idim <= nd_v; idim++) {
                pottargsort[FA2(idim, i, nd_v)] = 0.0;
                gradtargsort[FA3(idim, 1, i, nd_v, 2)] = 0.0;
                gradtargsort[FA3(idim, 2, i, nd_v, 2)] = 0.0;
                hesstargsort[FA3(idim, 1, i, nd_v, 3)] = 0.0;
                hesstargsort[FA3(idim, 2, i, nd_v, 3)] = 0.0;
                hesstargsort[FA3(idim, 3, i, nd_v, 3)] = 0.0;
            }
        }
    }

    /* compute scaling factors and lengths of mp/local expansions */
    rscales = (double *)malloc((nlevels + 1) * sizeof(double));
    nterms = (fint *)malloc((nlevels + 1) * sizeof(fint));

    nmax = 0;
    *ier = 0;
    for (i = 0; i <= nlevels; i++) {
        /* Must match Fortran: min(abs(zk*boxsize(i)/(2*pi)), 1) with
           complex arithmetic. cabs(zk)*boxsize/(2*pi) gives different
           last-bit rounding. */
        {
            fcomplex ztmp = (*zk) * boxsize[i] / (2.0 * pi);
            rscales[i] = cabs(ztmp);
            if (rscales[i] > 1.0) rscales[i] = 1.0;
        }
        h2dterms_(&boxsize[i], zk, eps, &nterms[i], ier);
        if (nterms[i] > nmax) nmax = nterms[i];
    }

    /* allocate iaddr and temporary arrays */
    iaddr = (fint *)malloc(4 * nboxes * sizeof(fint));
    lmptmp = (2 * nmax + 1) * nd_v;
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
            dreorderf_(&twond, ns, dipvec, dipvecsort, isrc);
        }
        /* reorder targets */
        dreorderf_(&two, nt, targ, targsort, itarg);
    }

    /* allocate workspace for multipole / local / diag form expansions */
    {
        const fint *laddr_p = &itree[iptr[0] - 1];
        FNAME(h2dmpalloc)(&nd_v, laddr_p, iaddr, &nlevels, &lmptot, nterms);
    }

    *ier = 0;
    rmlexp = (double *)malloc(lmptot * sizeof(double));
    if (rmlexp == NULL) {
        *ier = 4;
        goto cleanup;
    }

    /* call the main fmm routine */
    {
        const fint *laddr_p = &itree[iptr[0] - 1];
        FNAME(hfmm2dmain)(&nd_v, eps,
                   zk, ns, sourcesort,
                   &ifcharge_v, chargesort,
                   &ifdipole_v, dipstrsort, dipvecsort,
                   nt, targsort, &nexpc, expc,
                   iaddr, rmlexp, mptemp, &lmptmp,
                   itree, &ltree, iptr, &ndiv, &nlevels,
                   &nboxes, iper, boxsize, rscales, tcenters,
                   laddr_p,
                   isrcse, itargse, iexpcse, nterms, &ntj,
                   &ifpgh_v, potsort, gradsort, hesssort,
                   &ifpghtarg_v, pottargsort, gradtargsort,
                   hesstargsort, jexps, &scj, &ifnear, timeinfo, ier);
    }

    if (*ier != 0) goto cleanup;

    /* resort the output arrays in input order */
    {
        fint twond = 2 * nd_v;
        fint fournd = 4 * nd_v;
        fint sixnd = 6 * nd_v;

        if (ifpgh_v == 1) {
            dreorderi_(&twond, ns, (const double *)potsort,
                       (double *)pot, isrc);
        }
        if (ifpgh_v == 2) {
            dreorderi_(&twond, ns, (const double *)potsort,
                       (double *)pot, isrc);
            dreorderi_(&fournd, ns, (const double *)gradsort,
                       (double *)grad, isrc);
        }
        if (ifpgh_v == 3) {
            dreorderi_(&twond, ns, (const double *)potsort,
                       (double *)pot, isrc);
            dreorderi_(&fournd, ns, (const double *)gradsort,
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
            dreorderi_(&fournd, nt, (const double *)gradtargsort,
                       (double *)gradtarg, itarg);
        }
        if (ifpghtarg_v == 3) {
            dreorderi_(&twond, nt, (const double *)pottargsort,
                       (double *)pottarg, itarg);
            dreorderi_(&fournd, nt, (const double *)gradtargsort,
                       (double *)gradtarg, itarg);
            dreorderi_(&sixnd, nt, (const double *)hesstargsort,
                       (double *)hesstarg, itarg);
        }
    }

cleanup:
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
    free(dipvecsort);
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
