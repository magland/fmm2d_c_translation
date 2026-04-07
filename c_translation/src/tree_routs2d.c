/*
 * tree_routs2d.c - C translation of src/common/tree_routs2d.f
 *
 * Common routines for generating and processing a level-restricted
 * quad tree in 2D for the FMM. These manage integer-index bookkeeping
 * (parent / child / neighbor / list arrays) on top of an explicit
 * tree representation.
 *
 * The print_tree routine is intentionally not translated (file I/O,
 * not reachable from any caller). All other 7 routines are present.
 *
 * All OpenMP directives (C$OMP, C$OMP$) are comments in the Fortran
 * source and are ignored here; the result is sequential C. Control
 * flow, operation order, and array indexing are preserved 1:1 from
 * the Fortran source so integer outputs match bit-for-bit.
 *
 * Cross-file dependency: tree_refine_boxes and tree_refine_boxes_flag
 * call cumsum. We call the BARE Fortran symbol `cumsum_` via an
 * extern forward declaration so that the diff tests isolate this
 * file: when linked alongside libfmm2d.a, cumsum_ resolves to the
 * gfortran-compiled reference; when the full library is built with
 * -DFMM2D_DROP_IN, cumsum_ resolves to the C port's drop-in symbol.
 */

#include "tree_routs2d.h"

/* Cross-file call: use bare Fortran symbol name (see header comment). */
extern void cumsum_(const fint *n, const fint *a, fint *b);

/*
 * laddr(k, ilev) where k is 1-based and ilev is 0-based; leading dim 2.
 */
#define LADDR(k, ilev) ((ilev) * 2 + ((k) - 1))


/*
 * tree_refine_boxes: refine boxes flagged in irefinebox. Each flagged
 * box gets 4 children appended past the current tail (nbctr). Updates
 * centers, iparent, nchild, ichild, ilevel.
 *
 * Fortran:
 *   integer irefinebox(nbloc), ichild(4,nboxes), nchild(nboxes), ...
 *   real *8 centers(2,nboxes), bs
 *   allocate isum(nbloc); call cumsum(nbloc, irefinebox, isum)
 *   for each refined box i: nbl = nbctr + (isum(i)-1)*4; add 4 kids at nbl+1..nbl+4
 *   after loop: nbctr += isum(nbloc)*4
 */
void FNAME(tree_refine_boxes)(const fint *irefinebox, const fint *nboxes,
                              const fint *ifirstbox, const fint *nbloc,
                              double *centers, const double *bs,
                              fint *nbctr, const fint *nlctr,
                              fint *ilevel, fint *iparent,
                              fint *nchild, fint *ichild)
{
    fint nboxes_v = *nboxes;
    fint nbloc_v = *nbloc;
    fint ifirstbox_v = *ifirstbox;
    double bs_v = *bs;
    fint nlctr_v = *nlctr;
    fint i, ibox, nbl, j, l, jbox, ii;
    fint *isum;

    (void)nboxes_v;

    isum = (fint *)malloc(nbloc_v * sizeof(fint));

    if (nbloc_v > 0) {
        cumsum_(&nbloc_v, irefinebox, isum);
    }

    for (i = 1; i <= nbloc_v; i++) {
        ibox = ifirstbox_v + i - 1;
        if (irefinebox[i - 1] == 1) {
            nbl = (*nbctr) + (isum[i - 1] - 1) * 4;

            nchild[ibox - 1] = 4;
            for (j = 1; j <= 4; j++) {
                ii = 2;
                if (j <= 2) ii = 1;
                jbox = nbl + j;
                centers[FA2(1, jbox, 2)] =
                    centers[FA2(1, ibox, 2)] + ((j % 2 == 0) ? 1 : -1) * bs_v / 2;
                centers[FA2(2, jbox, 2)] =
                    centers[FA2(2, ibox, 2)] + ((ii % 2 == 0) ? 1 : -1) * bs_v / 2;
                iparent[jbox - 1] = ibox;
                nchild[jbox - 1] = 0;
                for (l = 1; l <= 4; l++) {
                    ichild[FA2(l, jbox, 4)] = -1;
                }
                ichild[FA2(j, ibox, 4)] = jbox;
                ilevel[jbox - 1] = nlctr_v;
            }
        }
    }

    if (nbloc_v > 0) {
        *nbctr = (*nbctr) + isum[nbloc_v - 1] * 4;
    }

    free(isum);
}


/*
 * tree_copy: copy tree fields from source arrays to destination
 * arrays. Straight element-wise copy.
 */
void FNAME(tree_copy)(const fint *nb,
                      const double *centers, const fint *ilevel,
                      const fint *iparent, const fint *nchild,
                      const fint *ichild,
                      double *centers2, fint *ilevel2,
                      fint *iparent2, fint *nchild2,
                      fint *ichild2)
{
    fint nb_v = *nb;
    fint i;

    for (i = 1; i <= nb_v; i++) {
        centers2[FA2(1, i, 2)] = centers[FA2(1, i, 2)];
        centers2[FA2(2, i, 2)] = centers[FA2(2, i, 2)];
        ilevel2[i - 1] = ilevel[i - 1];
        iparent2[i - 1] = iparent[i - 1];
        nchild2[i - 1] = nchild[i - 1];
        ichild2[FA2(1, i, 4)] = ichild[FA2(1, i, 4)];
        ichild2[FA2(2, i, 4)] = ichild[FA2(2, i, 4)];
        ichild2[FA2(3, i, 4)] = ichild[FA2(3, i, 4)];
        ichild2[FA2(4, i, 4)] = ichild[FA2(4, i, 4)];
    }
}


/*
 * computecoll: compute colleague (same-level neighbor) lists. A box
 * jbox is a colleague of ibox if they share an edge or a vertex and
 * both are at the same level. The algorithm loops levels top-down
 * and, for each box, inspects the 4*(nnbors of parent) candidates
 * that are children of the parent's colleagues.
 */
void FNAME(computecoll)(const fint *nlevels, const fint *nboxes,
                        const fint *laddr, const double *boxsize,
                        const double *centers, const fint *iparent,
                        const fint *nchild, const fint *ichild,
                        const fint *iper,
                        fint *nnbors, fint *nbors)
{
    fint nlevels_v = *nlevels;
    fint nboxes_v = *nboxes;
    fint ilev, ibox, jbox, kbox, dad;
    fint i, j, ifirstbox, ilastbox;

    (void)nboxes_v;
    (void)iper;
    (void)nchild;

    /* Setting parameters for level = 0 */
    nnbors[0] = 1;
    nbors[FA2(1, 1, 9)] = 1;

    for (ilev = 1; ilev <= nlevels_v; ilev++) {
        ifirstbox = laddr[LADDR(1, ilev)];
        ilastbox = laddr[LADDR(2, ilev)];

        for (ibox = ifirstbox; ibox <= ilastbox; ibox++) {
            dad = iparent[ibox - 1];
            for (i = 1; i <= nnbors[dad - 1]; i++) {
                jbox = nbors[FA2(i, dad, 9)];
                for (j = 1; j <= 4; j++) {
                    kbox = ichild[FA2(j, jbox, 4)];
                    if (kbox > 0) {
                        if ((fabs(centers[FA2(1, kbox, 2)] -
                                  centers[FA2(1, ibox, 2)]) <=
                             1.05 * boxsize[ilev]) &&
                            (fabs(centers[FA2(2, kbox, 2)] -
                                  centers[FA2(2, ibox, 2)]) <=
                             1.05 * boxsize[ilev])) {
                            nnbors[ibox - 1] = nnbors[ibox - 1] + 1;
                            nbors[FA2(nnbors[ibox - 1], ibox, 9)] = kbox;
                        }
                    }
                }
            }
        }
    }
}


/*
 * updateflags: boxes currently flagged with iflag(ibox)==3 (flag++)
 * are re-checked; those that need subdivision are set to 1, those
 * that do not are set to 0. A box needs subdivision if any child of
 * one of its colleagues lies within a distance test of it.
 *
 * Fortran uses `goto 1111` to break out of the inner loops once the
 * box has been set to 1. We preserve that exactly with a C goto.
 */
void FNAME(updateflags)(const fint *curlev, const fint *nboxes,
                        const fint *nlevels, const fint *laddr,
                        const fint *nchild, const fint *ichild,
                        const fint *nnbors, const fint *nbors,
                        const double *centers, const double *boxsize,
                        fint *iflag)
{
    fint curlev_v = *curlev;
    fint nboxes_v = *nboxes;
    fint nlevels_v = *nlevels;
    fint i, j, ibox, jbox, kbox, ict;
    double distest, xdis, ydis;

    (void)nboxes_v;
    (void)nlevels_v;

    distest = 1.05 * (boxsize[curlev_v] + boxsize[curlev_v + 1]) / 2.0;

    for (ibox = laddr[LADDR(1, curlev_v)]; ibox <= laddr[LADDR(2, curlev_v)]; ibox++) {
        if (iflag[ibox - 1] == 3) {
            iflag[ibox - 1] = 0;
            for (i = 1; i <= nnbors[ibox - 1]; i++) {
                jbox = nbors[FA2(i, ibox, 9)];
                for (j = 1; j <= 4; j++) {
                    kbox = ichild[FA2(j, jbox, 4)];
                    if (kbox > 0) {
                        if (nchild[kbox - 1] > 0) {
                            xdis = centers[FA2(1, kbox, 2)] -
                                   centers[FA2(1, ibox, 2)];
                            ydis = centers[FA2(2, kbox, 2)] -
                                   centers[FA2(2, ibox, 2)];
                            ict = 0;
                            if (fabs(xdis) <= distest) ict = ict + 1;
                            if (fabs(ydis) <= distest) ict = ict + 1;
                            if (ict == 2) {
                                iflag[ibox - 1] = 1;
                                goto label_1111;
                            }
                        }
                    }
                }
            }
        label_1111:
            ;
        }
    }
}


/*
 * tree_refine_boxes_flag: like tree_refine_boxes, but driven by
 * iflag instead of an explicit irefinebox array. Any box with
 * iflag > 0 gets refined; the new children inherit iflag=3 if the
 * parent was iflag==1, or iflag=0 if the parent was iflag==2.
 * Sets ilevel(jbox)=nlctr+1 (NOTE: this differs from tree_refine_boxes
 * which sets nlctr; the Fortran source has the same distinction).
 */
void FNAME(tree_refine_boxes_flag)(fint *iflag, const fint *nboxes,
                                   const fint *ifirstbox, const fint *nbloc,
                                   double *centers, const double *bs,
                                   fint *nbctr, const fint *nlctr,
                                   fint *ilevel, fint *iparent,
                                   fint *nchild, fint *ichild)
{
    fint nboxes_v = *nboxes;
    fint nbloc_v = *nbloc;
    fint ifirstbox_v = *ifirstbox;
    double bs_v = *bs;
    fint nlctr_v = *nlctr;
    fint i, ibox, nbl, j, l, jbox, ii;
    fint *isum;
    fint *itmp;

    (void)nboxes_v;

    isum = (fint *)malloc(nbloc_v * sizeof(fint));
    itmp = (fint *)malloc(nbloc_v * sizeof(fint));

    for (i = 1; i <= nbloc_v; i++) {
        ibox = ifirstbox_v + i - 1;
        itmp[i - 1] = 0;
        if (iflag[ibox - 1] > 0) itmp[i - 1] = 1;
    }

    if (nbloc_v > 0) {
        cumsum_(&nbloc_v, itmp, isum);
    }

    for (i = 1; i <= nbloc_v; i++) {
        ibox = ifirstbox_v + i - 1;
        if (iflag[ibox - 1] > 0) {
            nbl = (*nbctr) + (isum[i - 1] - 1) * 4;

            nchild[ibox - 1] = 4;
            for (j = 1; j <= 4; j++) {
                ii = 2;
                if (j <= 2) ii = 1;
                jbox = nbl + j;
                centers[FA2(1, jbox, 2)] =
                    centers[FA2(1, ibox, 2)] + ((j % 2 == 0) ? 1 : -1) * bs_v / 2;
                centers[FA2(2, jbox, 2)] =
                    centers[FA2(2, ibox, 2)] + ((ii % 2 == 0) ? 1 : -1) * bs_v / 2;
                iparent[jbox - 1] = ibox;
                nchild[jbox - 1] = 0;
                for (l = 1; l <= 4; l++) {
                    ichild[FA2(l, jbox, 4)] = -1;
                }
                ichild[FA2(j, ibox, 4)] = jbox;
                ilevel[jbox - 1] = nlctr_v + 1;
                if (iflag[ibox - 1] == 1) iflag[jbox - 1] = 3;
                if (iflag[ibox - 1] == 2) iflag[jbox - 1] = 0;
            }
        }
    }

    if (nbloc_v > 0) {
        *nbctr = (*nbctr) + isum[nbloc_v - 1] * 4;
    }

    free(isum);
    free(itmp);
}


/*
 * computemnlists: set max list sizes for list1..list4. In 2D these
 * are hardcoded constants; the routine does not examine the tree.
 */
void FNAME(computemnlists)(const fint *nlevels, const fint *nboxes,
                           const fint *itree, const fint *ltree,
                           const fint *iptr, const double *centers,
                           const double *boxsize, const fint *iper,
                           fint *mnlist1, fint *mnlist2,
                           fint *mnlist3, fint *mnlist4)
{
    (void)nlevels;
    (void)nboxes;
    (void)itree;
    (void)ltree;
    (void)iptr;
    (void)centers;
    (void)boxsize;
    (void)iper;

    *mnlist1 = 13;
    *mnlist2 = 27;
    *mnlist3 = 20;
    *mnlist4 = 5;
}


/*
 * computelists: compute FMM lists 1, 2, 3, 4 from a packed itree.
 * The packed itree contains (at offsets into iptr):
 *   iptr(1)  laddr(2, 0:nlevels)    - level box range table
 *   iptr(2)  ilevel(nboxes)
 *   iptr(3)  iparent(nboxes)
 *   iptr(4)  nchild(nboxes)
 *   iptr(5)  ichild(4, nboxes)
 *   iptr(6)  nnbors(nboxes)
 *   iptr(7)  nbors(9, nboxes)
 *
 * Note the Fortran uses `itree(2*ilev+1)` and `itree(2*ilev+2)` to
 * read laddr(1, ilev) and laddr(2, ilev) directly, taking advantage
 * of the fact that iptr(1) == 1 and laddr's (1,0) element lives at
 * itree(1) (so laddr(1, ilev) is at itree(2*ilev+1)). We preserve
 * that exactly.
 *
 * Note: the Fortran source has `cc goto 1120` (commented out) plus a
 * `1120 continue` label at the bottom of the inner loop. We keep the
 * continue label as a no-op to match structure.
 */
void FNAME(computelists)(const fint *nlevels, const fint *nboxes,
                         const fint *itree, const fint *ltree,
                         const fint *iptr, const double *centers,
                         const double *boxsize, const fint *iper,
                         const fint *mnlist1, fint *nlist1, fint *list1,
                         const fint *mnlist2, fint *nlist2, fint *list2,
                         const fint *mnlist3, fint *nlist3, fint *list3,
                         const fint *mnlist4, fint *nlist4, fint *list4)
{
    fint nlevels_v = *nlevels;
    fint nboxes_v = *nboxes;
    fint mnlist1_v = *mnlist1;
    fint mnlist2_v = *mnlist2;
    fint mnlist3_v = *mnlist3;
    fint mnlist4_v = *mnlist4;
    fint ilev, ibox, jbox, kbox, dad;
    fint i, j, ifirstbox, ilastbox;
    double distest, xdis, ydis;

    (void)ltree;
    (void)iper;

    for (i = 1; i <= nboxes_v; i++) {
        nlist1[i - 1] = 0;
        nlist2[i - 1] = 0;
        nlist3[i - 1] = 0;
        nlist4[i - 1] = 0;
    }

    /*
     * Setting parameters for level = 0
     *
     * Fortran: if (itree(iptr(4)).eq.0) ...   (i.e. itree at offset
     * iptr(4), which is the first element of nchild, i.e. nchild(1))
     */
    if (itree[iptr[3] - 1] == 0) {
        nlist1[0] = 1;
        list1[FA2(1, 1, mnlist1_v)] = 1;
    } else {
        nlist1[0] = 0;
    }
    nlist2[0] = 0;
    nlist3[0] = 0;
    nlist4[0] = 0;

    for (ilev = 1; ilev <= nlevels_v; ilev++) {
        /*
         * Fortran: ifirstbox = itree(2*ilev+1); ilastbox = itree(2*ilev+2)
         * This reads laddr(1,ilev) and laddr(2,ilev) from the packed itree.
         * The 1-based Fortran index k is k in C index k-1.
         */
        ifirstbox = itree[2 * ilev + 1 - 1];
        ilastbox = itree[2 * ilev + 2 - 1];

        for (ibox = ifirstbox; ibox <= ilastbox; ibox++) {
            /* dad = itree(iptr(3) + ibox - 1) */
            dad = itree[iptr[2] + ibox - 1 - 1];

            /* Loop over colleagues of dad; list 2 candidates are
             * children of dad's colleagues that are far from ibox. */
            /* Fortran: do i = 1, itree(iptr(6)+dad-1) */
            {
                fint nnb_dad = itree[iptr[5] + dad - 1 - 1];
                for (i = 1; i <= nnb_dad; i++) {
                    /* jbox = itree(iptr(7) + 9*(dad-1) + i - 1) */
                    jbox = itree[iptr[6] + 9 * (dad - 1) + i - 1 - 1];
                    {
                        fint nch_jbox = itree[iptr[3] + jbox - 1 - 1];
                        for (j = 1; j <= nch_jbox; j++) {
                            /* kbox = itree(iptr(5) + 4*(jbox-1) + j - 1) */
                            kbox = itree[iptr[4] + 4 * (jbox - 1) + j - 1 - 1];

                            if ((fabs(centers[FA2(1, kbox, 2)] -
                                      centers[FA2(1, ibox, 2)]) >=
                                 1.05 * boxsize[ilev]) ||
                                (fabs(centers[FA2(2, kbox, 2)] -
                                      centers[FA2(2, ibox, 2)]) >=
                                 1.05 * boxsize[ilev])) {
                                nlist2[ibox - 1] = nlist2[ibox - 1] + 1;
                                list2[FA2(nlist2[ibox - 1], ibox, mnlist2_v)] = kbox;
                            }
                        }
                    }
                }
            }

            /* cc goto 1120 (commented out in Fortran source) */

            /* Compute list 1 and list 3 of ibox if ibox is childless.
             * Fortran: if (itree(iptr(4)+ibox-1).eq.0) then */
            if (itree[iptr[3] + ibox - 1 - 1] == 0) {
                /* Loop over all colleagues of ibox.
                 * Fortran: do i = 1, itree(iptr(6)+ibox-1) */
                fint nnb_ibox = itree[iptr[5] + ibox - 1 - 1];
                for (i = 1; i <= nnb_ibox; i++) {
                    /* jbox = itree(iptr(7) + 9*(ibox-1) + i - 1) */
                    jbox = itree[iptr[6] + 9 * (ibox - 1) + i - 1 - 1];

                    /* If jbox is childless, it goes in list 1.
                     * Fortran: if (itree(iptr(4)+jbox-1).eq.0) then */
                    if (itree[iptr[3] + jbox - 1 - 1] == 0) {
                        nlist1[ibox - 1] = nlist1[ibox - 1] + 1;
                        list1[FA2(nlist1[ibox - 1], ibox, mnlist1_v)] = jbox;
                    } else {
                        distest = 1.05 * (boxsize[ilev] + boxsize[ilev + 1]) / 2.0;
                        /* Loop over children of colleague box jbox.
                         * Fortran: do j = 1, itree(iptr(4)+jbox-1) */
                        {
                            fint nch_jbox2 = itree[iptr[3] + jbox - 1 - 1];
                            for (j = 1; j <= nch_jbox2; j++) {
                                /* kbox = itree(iptr(5)+4*(jbox-1)+j-1) */
                                kbox = itree[iptr[4] + 4 * (jbox - 1) + j - 1 - 1];
                                xdis = fabs(centers[FA2(1, kbox, 2)] -
                                            centers[FA2(1, ibox, 2)]);
                                ydis = fabs(centers[FA2(2, kbox, 2)] -
                                            centers[FA2(2, ibox, 2)]);
                                if (xdis < distest && ydis < distest) {
                                    nlist1[ibox - 1] = nlist1[ibox - 1] + 1;
                                    list1[FA2(nlist1[ibox - 1], ibox, mnlist1_v)] = kbox;

                                    nlist1[kbox - 1] = nlist1[kbox - 1] + 1;
                                    list1[FA2(nlist1[kbox - 1], kbox, mnlist1_v)] = ibox;
                                } else {
                                    nlist3[ibox - 1] = nlist3[ibox - 1] + 1;
                                    list3[FA2(nlist3[ibox - 1], ibox, mnlist3_v)] = kbox;

                                    nlist4[kbox - 1] = nlist4[kbox - 1] + 1;
                                    list4[FA2(nlist4[kbox - 1], kbox, mnlist4_v)] = ibox;
                                }
                            }
                        }
                    }
                }
            }

            /* Fortran label 1120 continue (target of a commented-out
             * goto; preserved as a no-op for structural fidelity). */
            /* label_1120: (unused) */
            ;
        }
    }
}
