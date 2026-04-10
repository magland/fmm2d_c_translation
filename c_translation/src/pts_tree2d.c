/*
 * pts_tree2d.c - C translation of src/common/pts_tree2d.f
 *
 * Routines for building a level-restricted quad-tree from a set of
 * source / target points in 2D. Translated 1:1 from the Fortran
 * reference: same control flow, same operation order, same allocations,
 * same loop bounds. OpenMP directives (C$OMP, C$OMP$) in the source
 * are comments and are ignored here (sequential C).
 *
 * Cross-file dependency on tree_routs2d: tree_copy, tree_refine_boxes,
 * computecoll, tree_refine_boxes_flag, updateflags. These are called
 * via bare Fortran symbol names (extern declarations below) so the
 * diff tests isolate this file: when linked alongside libfmm2d.a,
 * these symbols resolve to the gfortran-compiled reference; when the
 * full library is built with -DFMM2D_DROP_IN they resolve to the C
 * port's drop-in symbols.
 */

#include "pts_tree2d.h"

/* Cross-file calls: use BARE Fortran symbol names. */
extern void tree_copy_(const fint *nb,
                       const double *centers, const fint *ilevel,
                       const fint *iparent, const fint *nchild,
                       const fint *ichild,
                       double *centers2, fint *ilevel2,
                       fint *iparent2, fint *nchild2,
                       fint *ichild2);

extern void tree_refine_boxes_(const fint *irefinebox, const fint *nboxes,
                               const fint *ifirstbox, const fint *nbloc,
                               double *centers, const double *bs,
                               fint *nbctr, const fint *nlctr,
                               fint *ilevel, fint *iparent,
                               fint *nchild, fint *ichild);

extern void computecoll_(const fint *nlevels, const fint *nboxes,
                         const fint *laddr, const double *boxsize,
                         const double *centers, const fint *iparent,
                         const fint *nchild, const fint *ichild,
                         const fint *iper,
                         fint *nnbors, fint *nbors);

extern void tree_refine_boxes_flag_(fint *iflag, const fint *nboxes,
                                    const fint *ifirstbox, const fint *nbloc,
                                    double *centers, const double *bs,
                                    fint *nbctr, const fint *nlctr,
                                    fint *ilevel, fint *iparent,
                                    fint *nchild, fint *ichild);

extern void updateflags_(const fint *curlev, const fint *nboxes,
                         const fint *nlevels, const fint *laddr,
                         const fint *nchild, const fint *ichild,
                         const fint *nnbors, const fint *nbors,
                         const double *centers, const double *boxsize,
                         fint *iflag);

/*
 * laddr(k, ilev): k is 1-based, ilev is 0-based, leading dim 2.
 * Used for laddr(2, 0:nlmax), laddrtail(2, 0:nlmax), tladdr(2, 0:nlmax).
 */
#define LADDR(k, ilev) ((ilev) * 2 + ((k) - 1))


/* ---------------------------------------------------------------- */
/* pts_tree_mem - get memory requirements for the tree              */
/* ---------------------------------------------------------------- */
void FNAME(pts_tree_mem)(const double *src, const fint *ns,
                         const double *targ, const fint *nt,
                         const fint *idivflag, const fint *ndiv,
                         const fint *nlmin, const fint *nlmax,
                         const fint *ifunif, const fint *iper,
                         fint *nlevels, fint *nboxes, fint *ltree)
{
    fint ns_v = *ns;
    fint nt_v = *nt;
    fint idivflag_v = *idivflag;
    fint ndiv_v = *ndiv;
    fint nlmin_v = *nlmin;
    fint nlmax_v = *nlmax;
    fint ifunif_v = *ifunif;
    fint iper_v = *iper;

    fint nbmax;
    fint i, j;
    fint nbloc, nbctr, nbadd, irefine, ilev, ifirstbox, ilastbox;
    fint nbtot;
    /* nn is set by one of three branches keyed on idivflag (0, 1, or 2);
       initializing to 0 silences a -Wmaybe-uninitialized warning and
       gives deterministic behavior for any out-of-range idivflag. */
    fint ibox, nn = 0, nss, ntt;
    double sizey;
    double xmin, xmax, ymin, ymax;
    double xmin2, xmax2, ymin2, ymax2;

    fint *laddr = NULL;
    fint *ilevel = NULL;
    fint *iparent = NULL;
    fint *nchild = NULL;
    fint *ichild = NULL;
    double *centers = NULL;
    fint *nbors = NULL;
    fint *nnbors = NULL;
    fint *isrc = NULL;
    fint *itarg = NULL;
    fint *isrcse = NULL;
    fint *itargse = NULL;
    fint *ilevel2 = NULL;
    fint *iparent2 = NULL;
    fint *nchild2 = NULL;
    fint *ichild2 = NULL;
    fint *isrcse2 = NULL;
    fint *itargse2 = NULL;
    double *centers2 = NULL;
    double *boxsize = NULL;
    fint *irefinebox = NULL;

    (void)iper_v;

    nbmax = 100000;

    boxsize = (double *)malloc((nlmax_v + 1) * sizeof(double));

    laddr = (fint *)malloc(2 * (nlmax_v + 1) * sizeof(fint));
    ilevel = (fint *)malloc(nbmax * sizeof(fint));
    iparent = (fint *)malloc(nbmax * sizeof(fint));
    nchild = (fint *)malloc(nbmax * sizeof(fint));
    ichild = (fint *)malloc(4 * nbmax * sizeof(fint));

    centers = (double *)malloc(2 * nbmax * sizeof(double));
    isrcse = (fint *)malloc(2 * nbmax * sizeof(fint));
    itargse = (fint *)malloc(2 * nbmax * sizeof(fint));
    isrc = (fint *)malloc(ns_v * sizeof(fint));
    itarg = (fint *)malloc(nt_v * sizeof(fint));

    /* step 1: find enclosing box */
    xmin = src[FA2(1, 1, 2)];
    xmax = src[FA2(1, 1, 2)];
    ymin = src[FA2(2, 1, 2)];
    ymax = src[FA2(2, 1, 2)];

    for (i = 1; i <= ns_v; i++) {
        if (src[FA2(1, i, 2)] < xmin) xmin = src[FA2(1, i, 2)];
        if (src[FA2(1, i, 2)] > xmax) xmax = src[FA2(1, i, 2)];
        if (src[FA2(2, i, 2)] < ymin) ymin = src[FA2(2, i, 2)];
        if (src[FA2(2, i, 2)] > ymax) ymax = src[FA2(2, i, 2)];
        isrc[i - 1] = i;
    }

    xmin2 = xmin;
    xmax2 = xmax;
    ymin2 = ymin;
    ymax2 = ymax;
    for (i = 1; i <= nt_v; i++) {
        if (targ[FA2(1, i, 2)] < xmin2) xmin2 = targ[FA2(1, i, 2)];
        if (targ[FA2(1, i, 2)] > xmax2) xmax2 = targ[FA2(1, i, 2)];
        if (targ[FA2(2, i, 2)] < ymin2) ymin2 = targ[FA2(2, i, 2)];
        if (targ[FA2(2, i, 2)] > ymax2) ymax2 = targ[FA2(2, i, 2)];
        itarg[i - 1] = i;
    }

    if (xmax2 > xmax) xmax = xmax2;
    if (xmin2 < xmin) xmin = xmin2;
    if (ymax2 > ymax) ymax = ymax2;
    if (ymin2 < ymin) ymin = ymin2;

    boxsize[0] = (xmax - xmin);
    sizey = (ymax - ymin);
    if (sizey > boxsize[0]) boxsize[0] = sizey;

    /* set tree info for level 0 */
    laddr[LADDR(1, 0)] = 1;
    laddr[LADDR(2, 0)] = 1;
    ilevel[0] = 0;
    iparent[0] = -1;
    nchild[0] = 0;
    for (i = 1; i <= 4; i++) {
        ichild[FA2(i, 1, 4)] = -1;
    }

    centers[FA2(1, 1, 2)] = (xmin + xmax) / 2;
    centers[FA2(2, 1, 2)] = (ymin + ymax) / 2;

    isrcse[FA2(1, 1, 2)] = 1;
    isrcse[FA2(2, 1, 2)] = ns_v;

    itargse[FA2(1, 1, 2)] = 1;
    itargse[FA2(2, 1, 2)] = nt_v;

    nbctr = 1;

    ilev = 0; /* In case the loop never executes (nlmax_v == 0). */
    for (ilev = 0; ilev <= nlmax_v - 1; ilev++) {
        irefine = 0;

        ifirstbox = laddr[LADDR(1, ilev)];
        ilastbox = laddr[LADDR(2, ilev)];

        nbloc = ilastbox - ifirstbox + 1;

        irefinebox = (fint *)malloc(nbloc * sizeof(fint));

        /* determine which boxes need to be refined */
        irefine = 0;

        if (ilev >= nlmin_v) {
            for (i = 1; i <= nbloc; i++) {
                irefinebox[i - 1] = 0;
                ibox = ifirstbox + i - 1;
                nss = isrcse[FA2(2, ibox, 2)] - isrcse[FA2(1, ibox, 2)] + 1;
                ntt = itargse[FA2(2, ibox, 2)] - itargse[FA2(1, ibox, 2)] + 1;
                if (idivflag_v == 0) nn = nss;
                if (idivflag_v == 1) nn = ntt;
                if (idivflag_v == 2) nn = (ntt > nss) ? ntt : nss;

                if (nn > ndiv_v) irefinebox[i - 1] = 1;
            }
            irefine = 0;
            for (i = 1; i <= nbloc; i++) {
                if (irefinebox[i - 1] > irefine) irefine = irefinebox[i - 1];
            }
            if (ifunif_v == 1) {
                for (i = 1; i <= nbloc; i++) {
                    irefinebox[i - 1] = irefine;
                }
            }
        } else {
            irefine = 1;
            for (i = 1; i <= nbloc; i++) {
                irefinebox[i - 1] = 1;
            }
        }

        nbadd = 0;
        for (i = 1; i <= nbloc; i++) {
            if (irefinebox[i - 1] == 1) nbadd = nbadd + 4;
        }

        /* figure out if current allocation of boxes is sufficient */
        nbtot = nbctr + nbadd;

        /* if current memory is not sufficient reallocate */
        if (nbtot > nbmax) {
            centers2 = (double *)malloc(2 * nbmax * sizeof(double));
            ilevel2 = (fint *)malloc(nbmax * sizeof(fint));
            iparent2 = (fint *)malloc(nbmax * sizeof(fint));
            nchild2 = (fint *)malloc(nbmax * sizeof(fint));
            ichild2 = (fint *)malloc(4 * nbmax * sizeof(fint));
            isrcse2 = (fint *)malloc(2 * nbmax * sizeof(fint));
            itargse2 = (fint *)malloc(2 * nbmax * sizeof(fint));

            tree_copy_(&nbctr, centers, ilevel, iparent, nchild,
                       ichild, centers2, ilevel2, iparent2,
                       nchild2, ichild2);

            for (i = 1; i <= nbctr; i++) {
                isrcse2[FA2(1, i, 2)] = isrcse[FA2(1, i, 2)];
                isrcse2[FA2(2, i, 2)] = isrcse[FA2(2, i, 2)];
                itargse2[FA2(1, i, 2)] = itargse[FA2(1, i, 2)];
                itargse2[FA2(2, i, 2)] = itargse[FA2(2, i, 2)];
            }

            free(centers);
            free(ilevel);
            free(iparent);
            free(nchild);
            free(ichild);
            free(isrcse);
            free(itargse);

            nbmax = nbtot;
            centers = (double *)malloc(2 * nbmax * sizeof(double));
            ilevel = (fint *)malloc(nbmax * sizeof(fint));
            iparent = (fint *)malloc(nbmax * sizeof(fint));
            nchild = (fint *)malloc(nbmax * sizeof(fint));
            ichild = (fint *)malloc(4 * nbmax * sizeof(fint));
            isrcse = (fint *)malloc(2 * nbmax * sizeof(fint));
            itargse = (fint *)malloc(2 * nbmax * sizeof(fint));

            tree_copy_(&nbctr, centers2, ilevel2, iparent2,
                       nchild2, ichild2, centers, ilevel, iparent,
                       nchild, ichild);

            for (i = 1; i <= nbctr; i++) {
                isrcse[FA2(1, i, 2)] = isrcse2[FA2(1, i, 2)];
                isrcse[FA2(2, i, 2)] = isrcse2[FA2(2, i, 2)];
                itargse[FA2(1, i, 2)] = itargse2[FA2(1, i, 2)];
                itargse[FA2(2, i, 2)] = itargse2[FA2(2, i, 2)];
            }

            free(centers2);
            free(ilevel2);
            free(iparent2);
            free(nchild2);
            free(ichild2);
            free(isrcse2);
            free(itargse2);
            centers2 = NULL;
            ilevel2 = NULL;
            iparent2 = NULL;
            nchild2 = NULL;
            ichild2 = NULL;
            isrcse2 = NULL;
            itargse2 = NULL;
        }

        if (irefine == 1) {
            fint ilevp1 = ilev + 1;
            boxsize[ilev + 1] = boxsize[ilev] / 2;
            laddr[LADDR(1, ilev + 1)] = nbctr + 1;

            tree_refine_boxes_(irefinebox, &nbmax,
                               &ifirstbox, &nbloc, centers, &boxsize[ilev + 1],
                               &nbctr, &ilevp1,
                               ilevel, iparent, nchild, ichild);

            for (i = 1; i <= nbloc; i++) {
                ibox = ifirstbox + i - 1;
                if (irefinebox[i - 1] == 1) {
                    FNAME(sort_pts_to_children)(&ibox, &nbmax, centers, ichild,
                                                src, &ns_v, isrc, isrcse);
                    FNAME(sort_pts_to_children)(&ibox, &nbmax, centers, ichild,
                                                targ, &nt_v, itarg, itargse);
                }
            }
            laddr[LADDR(2, ilev + 1)] = nbctr;
        } else {
            free(irefinebox);
            irefinebox = NULL;
            break;
        }

        free(irefinebox);
        irefinebox = NULL;
    }

    *nboxes = nbctr;
    *nlevels = ilev;

    if (*nlevels >= 2) {
        nbtot = 8 * (*nboxes);
        if (nbtot > nbmax) {
            centers2 = (double *)malloc(2 * nbmax * sizeof(double));
            ilevel2 = (fint *)malloc(nbmax * sizeof(fint));
            iparent2 = (fint *)malloc(nbmax * sizeof(fint));
            nchild2 = (fint *)malloc(nbmax * sizeof(fint));
            ichild2 = (fint *)malloc(4 * nbmax * sizeof(fint));
            isrcse2 = (fint *)malloc(2 * nbmax * sizeof(fint));
            itargse2 = (fint *)malloc(2 * nbmax * sizeof(fint));

            tree_copy_(&nbctr, centers, ilevel, iparent, nchild,
                       ichild, centers2, ilevel2, iparent2,
                       nchild2, ichild2);

            for (i = 1; i <= nbctr; i++) {
                isrcse2[FA2(1, i, 2)] = isrcse[FA2(1, i, 2)];
                isrcse2[FA2(2, i, 2)] = isrcse[FA2(2, i, 2)];
                itargse2[FA2(1, i, 2)] = itargse[FA2(1, i, 2)];
                itargse2[FA2(2, i, 2)] = itargse[FA2(2, i, 2)];
            }

            free(centers);
            free(ilevel);
            free(iparent);
            free(nchild);
            free(ichild);
            free(isrcse);
            free(itargse);

            nbmax = nbtot;
            centers = (double *)malloc(2 * nbmax * sizeof(double));
            ilevel = (fint *)malloc(nbmax * sizeof(fint));
            iparent = (fint *)malloc(nbmax * sizeof(fint));
            nchild = (fint *)malloc(nbmax * sizeof(fint));
            ichild = (fint *)malloc(4 * nbmax * sizeof(fint));
            isrcse = (fint *)malloc(2 * nbmax * sizeof(fint));
            itargse = (fint *)malloc(2 * nbmax * sizeof(fint));

            tree_copy_(&nbctr, centers2, ilevel2, iparent2,
                       nchild2, ichild2, centers, ilevel, iparent,
                       nchild, ichild);
            for (i = 1; i <= nbctr; i++) {
                isrcse[FA2(1, i, 2)] = isrcse2[FA2(1, i, 2)];
                isrcse[FA2(2, i, 2)] = isrcse2[FA2(2, i, 2)];
                itargse[FA2(1, i, 2)] = itargse2[FA2(1, i, 2)];
                itargse[FA2(2, i, 2)] = itargse2[FA2(2, i, 2)];
            }

            free(centers2);
            free(ilevel2);
            free(iparent2);
            free(nchild2);
            free(ichild2);
            free(isrcse2);
            free(itargse2);
            centers2 = NULL;
            ilevel2 = NULL;
            iparent2 = NULL;
            nchild2 = NULL;
            ichild2 = NULL;
            isrcse2 = NULL;
            itargse2 = NULL;
        }

        nnbors = (fint *)malloc(nbmax * sizeof(fint));
        nbors = (fint *)malloc(9 * nbmax * sizeof(fint));

        for (i = 1; i <= *nboxes; i++) {
            nnbors[i - 1] = 0;
            for (j = 1; j <= 9; j++) {
                nbors[FA2(j, i, 9)] = -1;
            }
        }

        computecoll_(nlevels, nboxes, laddr, boxsize, centers,
                     iparent, nchild, ichild, iper, nnbors, nbors);

        if (*nlevels >= 2) {
            FNAME(pts_tree_fix_lr)(centers, nlevels, nboxes, boxsize, &nbmax,
                                   nlmax, iper, laddr, ilevel, iparent,
                                   nchild, ichild, nnbors, nbors);
        }
    }

    *ltree = 17 * (*nboxes) + 2 * (*nlevels + 1);

    /* Free everything */
    if (irefinebox) free(irefinebox);
    free(boxsize);
    free(laddr);
    free(ilevel);
    free(iparent);
    free(nchild);
    free(ichild);
    free(centers);
    free(isrcse);
    free(itargse);
    free(isrc);
    free(itarg);
    if (nnbors) free(nnbors);
    if (nbors) free(nbors);

    return;
}


/* ---------------------------------------------------------------- */
/* pts_tree_build - build the tree                                  */
/* ---------------------------------------------------------------- */
void FNAME(pts_tree_build)(const double *src, const fint *ns,
                           const double *targ, const fint *nt,
                           const fint *idivflag, const fint *ndiv,
                           const fint *nlmin, const fint *nlmax,
                           const fint *ifunif, const fint *iper,
                           fint *nlevels, const fint *nboxes,
                           const fint *ltree, fint *itree, fint *iptr,
                           double *centers, double *boxsize)
{
    fint ns_v = *ns;
    fint nt_v = *nt;
    fint idivflag_v = *idivflag;
    fint ndiv_v = *ndiv;
    fint nlmin_v = *nlmin;
    fint nlmax_v = *nlmax;
    fint ifunif_v = *ifunif;
    fint nboxes_v = *nboxes;

    fint i, ilev, irefine, j;
    fint ifirstbox, ilastbox, nbctr, nbloc;
    fint nboxes0;
    /* See note on nn in pts_tree_mem above. */
    fint ibox, nn = 0, nss, ntt;

    double xmin, xmax, ymin, ymax, sizey;

    fint *irefinebox = NULL;
    fint *isrc = NULL;
    fint *itarg = NULL;
    fint *isrcse = NULL;
    fint *itargse = NULL;

    (void)nlmax_v;

    iptr[0] = 1;
    iptr[1] = 2 * (*nlevels + 1) + 1;
    iptr[2] = iptr[1] + nboxes_v;
    iptr[3] = iptr[2] + nboxes_v;
    iptr[4] = iptr[3] + nboxes_v;
    iptr[5] = iptr[4] + 4 * nboxes_v;
    iptr[6] = iptr[5] + nboxes_v;
    iptr[7] = iptr[6] + 9 * nboxes_v;

    xmin = src[FA2(1, 1, 2)];
    xmax = src[FA2(1, 1, 2)];
    ymin = src[FA2(2, 1, 2)];
    ymax = src[FA2(2, 1, 2)];

    for (i = 1; i <= ns_v; i++) {
        if (src[FA2(1, i, 2)] < xmin) xmin = src[FA2(1, i, 2)];
        if (src[FA2(1, i, 2)] > xmax) xmax = src[FA2(1, i, 2)];
        if (src[FA2(2, i, 2)] < ymin) ymin = src[FA2(2, i, 2)];
        if (src[FA2(2, i, 2)] > ymax) ymax = src[FA2(2, i, 2)];
    }

    for (i = 1; i <= nt_v; i++) {
        if (targ[FA2(1, i, 2)] < xmin) xmin = targ[FA2(1, i, 2)];
        if (targ[FA2(1, i, 2)] > xmax) xmax = targ[FA2(1, i, 2)];
        if (targ[FA2(2, i, 2)] < ymin) ymin = targ[FA2(2, i, 2)];
        if (targ[FA2(2, i, 2)] > ymax) ymax = targ[FA2(2, i, 2)];
    }
    boxsize[0] = xmax - xmin;
    sizey = ymax - ymin;
    if (sizey > boxsize[0]) boxsize[0] = sizey;

    centers[FA2(1, 1, 2)] = (xmin + xmax) / 2;
    centers[FA2(2, 1, 2)] = (ymin + ymax) / 2;

    isrc = (fint *)malloc(ns_v * sizeof(fint));
    itarg = (fint *)malloc(nt_v * sizeof(fint));
    isrcse = (fint *)malloc(2 * nboxes_v * sizeof(fint));
    itargse = (fint *)malloc(2 * nboxes_v * sizeof(fint));

    /* set tree info for level 0 */
    itree[0] = 1;
    itree[1] = 1;
    itree[iptr[1] - 1] = 0;          /* itree(iptr(2)) */
    itree[iptr[2] - 1] = -1;         /* itree(iptr(3)) */
    itree[iptr[3] - 1] = 0;          /* itree(iptr(4)) */
    for (i = 1; i <= 4; i++) {
        itree[iptr[4] + i - 1 - 1] = -1;   /* itree(iptr(5)+i-1) */
    }

    isrcse[FA2(1, 1, 2)] = 1;
    isrcse[FA2(2, 1, 2)] = ns_v;
    itargse[FA2(1, 1, 2)] = 1;
    itargse[FA2(2, 1, 2)] = nt_v;

    for (i = 1; i <= ns_v; i++) {
        isrc[i - 1] = i;
    }

    for (i = 1; i <= nt_v; i++) {
        itarg[i - 1] = i;
    }

    /* Reset nlevels, nboxes */
    nbctr = 1;

    ilev = 0; /* In case the loop never executes. */
    for (ilev = 0; ilev <= *nlevels - 1; ilev++) {
        irefine = 0;

        ifirstbox = itree[2 * ilev];       /* itree(2*ilev+1) */
        ilastbox = itree[2 * ilev + 1];    /* itree(2*ilev+2) */

        nbloc = ilastbox - ifirstbox + 1;
        irefinebox = (fint *)malloc(nbloc * sizeof(fint));

        /* determine which boxes need to be refined */
        if (ilev >= nlmin_v) {
            for (i = 1; i <= nbloc; i++) {
                irefinebox[i - 1] = 0;
                ibox = ifirstbox + i - 1;
                nss = isrcse[FA2(2, ibox, 2)] - isrcse[FA2(1, ibox, 2)] + 1;
                ntt = itargse[FA2(2, ibox, 2)] - itargse[FA2(1, ibox, 2)] + 1;

                if (idivflag_v == 0) nn = nss;
                if (idivflag_v == 1) nn = ntt;
                if (idivflag_v == 2) nn = (ntt > nss) ? ntt : nss;

                if (nn > ndiv_v) irefinebox[i - 1] = 1;
            }
            irefine = 0;
            for (i = 1; i <= nbloc; i++) {
                if (irefinebox[i - 1] > irefine) irefine = irefinebox[i - 1];
            }
            if (ifunif_v == 1) {
                for (i = 1; i <= nbloc; i++) {
                    irefinebox[i - 1] = irefine;
                }
            }
        } else {
            for (i = 1; i <= nbloc; i++) {
                irefinebox[i - 1] = 1;
            }
        }

        irefine = 0;
        for (i = 1; i <= nbloc; i++) {
            if (irefinebox[i - 1] > irefine) irefine = irefinebox[i - 1];
        }

        if (irefine == 1) {
            fint ilevp1 = ilev + 1;
            boxsize[ilev + 1] = boxsize[ilev] / 2;
            itree[2 * ilev + 2] = nbctr + 1;       /* itree(2*ilev+3) */

            tree_refine_boxes_(irefinebox, nboxes, &ifirstbox, &nbloc,
                               centers, &boxsize[ilev + 1], &nbctr, &ilevp1,
                               &itree[iptr[1] - 1],
                               &itree[iptr[2] - 1],
                               &itree[iptr[3] - 1],
                               &itree[iptr[4] - 1]);

            /* re sort points in refined boxes */
            for (i = 1; i <= nbloc; i++) {
                ibox = ifirstbox + i - 1;
                if (irefinebox[i - 1] == 1) {
                    FNAME(sort_pts_to_children)(&ibox, nboxes, centers,
                                                &itree[iptr[4] - 1],
                                                src, &ns_v, isrc, isrcse);
                    FNAME(sort_pts_to_children)(&ibox, nboxes, centers,
                                                &itree[iptr[4] - 1],
                                                targ, &nt_v, itarg, itargse);
                }
            }

            itree[2 * ilev + 3] = nbctr;           /* itree(2*ilev+4) */
        } else {
            free(irefinebox);
            irefinebox = NULL;
            break;
        }

        free(irefinebox);
        irefinebox = NULL;
    }

    nboxes0 = nbctr;
    *nlevels = ilev;

    for (i = 1; i <= nboxes0; i++) {
        itree[iptr[5] + i - 1 - 1] = 0;              /* itree(iptr(6)+i-1) */
        for (j = 1; j <= 9; j++) {
            itree[iptr[6] + 9 * (i - 1) + j - 1 - 1] = -1; /* itree(iptr(7)+9*(i-1)+j-1) */
        }
    }

    computecoll_(nlevels, &nboxes0, &itree[iptr[0] - 1], boxsize, centers,
                 &itree[iptr[2] - 1], &itree[iptr[3] - 1], &itree[iptr[4] - 1],
                 iper, &itree[iptr[5] - 1], &itree[iptr[6] - 1]);

    if (*nlevels >= 2) {
        FNAME(pts_tree_fix_lr)(centers, nlevels,
                               &nboxes0, boxsize, nboxes, nlevels, iper,
                               &itree[iptr[0] - 1],
                               &itree[iptr[1] - 1],
                               &itree[iptr[2] - 1],
                               &itree[iptr[3] - 1],
                               &itree[iptr[4] - 1],
                               &itree[iptr[5] - 1],
                               &itree[iptr[6] - 1]);
    }

    if (irefinebox) free(irefinebox);
    free(isrc);
    free(itarg);
    free(isrcse);
    free(itargse);

    (void)ltree;

    return;
}


/* ---------------------------------------------------------------- */
/* sort_pts_to_children - partition points of box ibox into children */
/* ---------------------------------------------------------------- */
void FNAME(sort_pts_to_children)(const fint *ibox, const fint *nboxes,
                                 const double *centers, const fint *ichild,
                                 const double *src, const fint *ns,
                                 fint *isrc, fint *isrcse)
{
    fint ibox_v = *ibox;
    fint nboxes_v = *nboxes;
    fint ns_v = *ns;

    fint i12, i34, npts, iss, i, jbox, istart;
    fint nsc[4];
    fint *isrctmp = NULL;

    (void)nboxes_v;
    (void)ns_v;

    /*
     * Allocate temporary array to figure out which child you
     * belong to?  1,2,3,4? counter ns1,ns2,ns3,ns4
     * The box nomenclature is as follows
     *     3   4
     *     1   2
     */

    i12 = isrcse[FA2(1, ibox_v, 2)] - 1;
    i34 = 0;
    npts = isrcse[FA2(2, ibox_v, 2)] - isrcse[FA2(1, ibox_v, 2)] + 1;
    isrctmp = (fint *)malloc(npts * sizeof(fint));

    for (iss = isrcse[FA2(1, ibox_v, 2)]; iss <= isrcse[FA2(2, ibox_v, 2)]; iss++) {
        if (src[FA2(2, isrc[iss - 1], 2)] - centers[FA2(2, ibox_v, 2)] < 0) {
            i12 = i12 + 1;
            isrc[i12 - 1] = isrc[iss - 1];
        } else {
            i34 = i34 + 1;
            isrctmp[i34 - 1] = isrc[iss - 1];
        }
    }

    /* reorder sources */
    for (i = 1; i <= i34; i++) {
        isrc[i12 + i - 1] = isrctmp[i - 1];
    }

    /* now sort into boxes 1,2 */
    nsc[0] = 0;
    nsc[1] = 0;
    nsc[2] = 0;
    nsc[3] = 0;

    /* sort into boxes 1,2 */
    for (iss = isrcse[FA2(1, ibox_v, 2)]; iss <= i12; iss++) {
        if (src[FA2(1, isrc[iss - 1], 2)] - centers[FA2(1, ibox_v, 2)] < 0) {
            isrc[isrcse[FA2(1, ibox_v, 2)] + nsc[0] - 1] = isrc[iss - 1];
            nsc[0] = nsc[0] + 1;
        } else {
            nsc[1] = nsc[1] + 1;
            isrctmp[nsc[1] - 1] = isrc[iss - 1];
        }
    }

    /* reorder sources so that sources in 2 are at end of array */
    for (i = 1; i <= nsc[1]; i++) {
        isrc[isrcse[FA2(1, ibox_v, 2)] + nsc[0] + i - 1 - 1] = isrctmp[i - 1];
    }

    /* sort into boxes 3,4 */
    for (iss = i12 + 1; iss <= isrcse[FA2(2, ibox_v, 2)]; iss++) {
        if (src[FA2(1, isrc[iss - 1], 2)] - centers[FA2(1, ibox_v, 2)] < 0) {
            isrc[i12 + 1 + nsc[2] - 1] = isrc[iss - 1];
            nsc[2] = nsc[2] + 1;
        } else {
            nsc[3] = nsc[3] + 1;
            isrctmp[nsc[3] - 1] = isrc[iss - 1];
        }
    }

    for (i = 1; i <= nsc[3]; i++) {
        isrc[i12 + nsc[2] + i - 1] = isrctmp[i - 1];
    }

    istart = isrcse[FA2(1, ibox_v, 2)];
    for (i = 1; i <= 4; i++) {
        jbox = ichild[FA2(i, ibox_v, 4)];
        isrcse[FA2(1, jbox, 2)] = istart;
        isrcse[FA2(2, jbox, 2)] = istart + nsc[i - 1] - 1;
        istart = istart + nsc[i - 1];
    }

    free(isrctmp);
    return;
}


/* ---------------------------------------------------------------- */
/* pts_tree_fix_lr - convert adaptive tree to level-restricted tree */
/* ---------------------------------------------------------------- */
void FNAME(pts_tree_fix_lr)(double *centers, const fint *nlevels,
                            fint *nboxes, const double *boxsize,
                            const fint *nbmax, const fint *nlmax,
                            const fint *iper, fint *laddr, fint *ilevel,
                            fint *iparent, fint *nchild, fint *ichild,
                            fint *nnbors, fint *nbors)
{
    fint nlevels_v = *nlevels;
    fint nbmax_v = *nbmax;
    fint nlmax_v = *nlmax;

    fint i, j, ibox, jbox, kbox, ilev, idad, igranddad;
    fint nbloc, ict;
    double xdis, ydis, distest;

    fint *iflag = NULL;
    fint *laddrtail = NULL;

    laddrtail = (fint *)malloc(2 * (nlmax_v + 1) * sizeof(fint));
    iflag = (fint *)malloc(nbmax_v * sizeof(fint));

    /* Initialize flag array */
    for (i = 1; i <= *nboxes; i++) {
        iflag[i - 1] = 0;
    }

    /*
     * Flag boxes that violate level restriction by "1"
     * For such boxes, we set iflag(i) = 1
     */
    for (ilev = nlevels_v; ilev >= 2; ilev--) {
        /* distance to test if two boxes separated by two levels are touching */
        distest = 1.05 * (boxsize[ilev - 1] + boxsize[ilev - 2]) / 2.0;
        for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            idad = iparent[ibox - 1];
            igranddad = iparent[idad - 1];

            /* Loop over colleagues of granddad */
            for (i = 1; i <= nnbors[igranddad - 1]; i++) {
                jbox = nbors[FA2(i, igranddad, 9)];
                /* Check if the colleague of grandad is a leaf node */
                if (nchild[jbox - 1] == 0 && iflag[jbox - 1] == 0) {
                    xdis = centers[FA2(1, jbox, 2)] - centers[FA2(1, idad, 2)];
                    ydis = centers[FA2(2, jbox, 2)] - centers[FA2(2, idad, 2)];
                    ict = 0;
                    if (fabs(xdis) <= distest) ict = ict + 1;
                    if (fabs(ydis) <= distest) ict = ict + 1;
                    if (ict == 2) {
                        iflag[jbox - 1] = 1;
                    }
                }
            }
        }
    }

    /*
     * Find all boxes that need to be given a flag+
     */
    for (ilev = nlevels_v; ilev >= 1; ilev--) {
        /* distance to test if two boxes separated by one level are touching */
        distest = 1.05 * (boxsize[ilev] + boxsize[ilev - 1]) / 2.0;
        for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            if (iflag[ibox - 1] == 1 || iflag[ibox - 1] == 2) {
                idad = iparent[ibox - 1];
                /* Loop over dad's colleagues */
                for (i = 1; i <= nnbors[idad - 1]; i++) {
                    jbox = nbors[FA2(i, idad, 9)];
                    if (nchild[jbox - 1] == 0 && iflag[jbox - 1] == 0) {
                        xdis = centers[FA2(1, jbox, 2)] - centers[FA2(1, ibox, 2)];
                        ydis = centers[FA2(2, jbox, 2)] - centers[FA2(2, ibox, 2)];
                        ict = 0;
                        if (fabs(xdis) <= distest) ict = ict + 1;
                        if (fabs(ydis) <= distest) ict = ict + 1;
                        if (ict == 2) {
                            iflag[jbox - 1] = 2;
                        }
                    }
                }
            }
        }
    }

    /*
     * Subdivide all flag and flag+ boxes.
     */
    for (ilev = 0; ilev <= nlevels_v; ilev++) {
        laddrtail[LADDR(1, ilev)] = 0;
        laddrtail[LADDR(2, ilev)] = -1;
    }

    for (ilev = 1; ilev <= nlevels_v - 2; ilev++) {
        laddrtail[LADDR(1, ilev + 1)] = *nboxes + 1;

        nbloc = laddr[LADDR(2, ilev)] - laddr[LADDR(1, ilev)] + 1;
        tree_refine_boxes_flag_(iflag, &nbmax_v, &laddr[LADDR(1, ilev)],
                                &nbloc, centers, &boxsize[ilev + 1], nboxes,
                                &ilev, ilevel, iparent, nchild, ichild);

        laddrtail[LADDR(2, ilev + 1)] = *nboxes;
    }

    /* Reorganize the tree to get it back in the standard format */
    FNAME(pts_tree_reorg)(nboxes, centers, nlevels, laddr,
                          laddrtail, ilevel, iparent, nchild, ichild,
                          iflag);

    /* Compute colleague information again */
    for (i = 1; i <= *nboxes; i++) {
        nnbors[i - 1] = 0;
        for (j = 1; j <= 9; j++) {
            nbors[FA2(j, i, 9)] = -1;
        }
    }
    computecoll_(nlevels, nboxes, laddr, boxsize,
                 centers, iparent, nchild, ichild, iper, nnbors, nbors);

    /*
     * Reset the flags array to remove all the flag and flag+ cases.
     */
    for (ibox = 1; ibox <= *nboxes; ibox++) {
        if (iflag[ibox - 1] != 3) iflag[ibox - 1] = 0;
    }

    for (ilev = 0; ilev <= nlevels_v; ilev++) {
        laddrtail[LADDR(1, ilev)] = 0;
        laddrtail[LADDR(2, ilev)] = -1;
    }

    for (ilev = 2; ilev <= nlevels_v - 2; ilev++) {
        /* Step 1: Determine which flag++ boxes need further division */
        updateflags_(&ilev, nboxes, nlevels, laddr, nchild, ichild,
                     nnbors, nbors, centers, boxsize, iflag);

        updateflags_(&ilev, nboxes, nlevels, laddrtail, nchild, ichild,
                     nnbors, nbors, centers, boxsize, iflag);

        /* Step 2: Subdivide all the boxes that need subdivision */
        laddrtail[LADDR(1, ilev + 1)] = *nboxes + 1;

        nbloc = laddr[LADDR(2, ilev)] - laddr[LADDR(1, ilev)] + 1;
        tree_refine_boxes_flag_(iflag, &nbmax_v, &laddr[LADDR(1, ilev)],
                                &nbloc, centers, &boxsize[ilev + 1], nboxes,
                                &ilev, ilevel, iparent, nchild, ichild);

        nbloc = laddrtail[LADDR(2, ilev)] - laddrtail[LADDR(1, ilev)] + 1;
        tree_refine_boxes_flag_(iflag, &nbmax_v, &laddrtail[LADDR(1, ilev)],
                                &nbloc, centers, &boxsize[ilev + 1], nboxes,
                                &ilev, ilevel, iparent, nchild, ichild);

        laddrtail[LADDR(2, ilev + 1)] = *nboxes;

        /* Step 3: Update colleague information for newly created boxes */
        for (ibox = laddrtail[LADDR(1, ilev + 1)];
             ibox <= laddrtail[LADDR(2, ilev + 1)]; ibox++) {
            nnbors[ibox - 1] = 0;
            idad = iparent[ibox - 1];
            for (i = 1; i <= nnbors[idad - 1]; i++) {
                jbox = nbors[FA2(i, idad, 9)];
                for (j = 1; j <= 4; j++) {
                    kbox = ichild[FA2(j, jbox, 4)];
                    if (kbox > 0) {
                        if ((fabs(centers[FA2(1, kbox, 2)] - centers[FA2(1, ibox, 2)])
                                 <= 1.05 * boxsize[ilev + 1]) &&
                            (fabs(centers[FA2(2, kbox, 2)] - centers[FA2(2, ibox, 2)])
                                 <= 1.05 * boxsize[ilev + 1])) {
                            nnbors[ibox - 1] = nnbors[ibox - 1] + 1;
                            nbors[FA2(nnbors[ibox - 1], ibox, 9)] = kbox;
                        }
                    }
                }
            }
        }
    }

    /* Reorganize tree once again and we are all done */
    FNAME(pts_tree_reorg)(nboxes, centers, nlevels, laddr,
                          laddrtail, ilevel, iparent, nchild, ichild,
                          iflag);

    /* Compute colleague information again */
    for (i = 1; i <= *nboxes; i++) {
        nnbors[i - 1] = 0;
        for (j = 1; j <= 9; j++) {
            nbors[FA2(j, i, 9)] = -1;
        }
    }

    computecoll_(nlevels, nboxes, laddr, boxsize,
                 centers, iparent, nchild, ichild, iper, nnbors, nbors);

    free(iflag);
    free(laddrtail);
    return;
}


/* ---------------------------------------------------------------- */
/* pts_tree_reorg - reorganize tree arrays after fix_lr             */
/* ---------------------------------------------------------------- */
void FNAME(pts_tree_reorg)(const fint *nboxes, double *centers,
                           const fint *nlevels, fint *laddr,
                           const fint *laddrtail, fint *ilevel,
                           fint *iparent, fint *nchild, fint *ichild,
                           fint *iflag)
{
    fint nboxes_v = *nboxes;
    fint nlevels_v = *nlevels;

    fint i, ibox, ilev, curbox, nblev;

    fint *tilevel = NULL;
    fint *tiparent = NULL;
    fint *tnchild = NULL;
    fint *tichild = NULL;
    fint *tiflag = NULL;
    fint *iboxtocurbox = NULL;
    fint *ilevptr = NULL;
    fint *ilevptr2 = NULL;
    double *tcenters = NULL;
    fint *tladdr = NULL;

    tladdr = (fint *)malloc(2 * (nlevels_v + 1) * sizeof(fint));

    tilevel = (fint *)malloc(nboxes_v * sizeof(fint));
    tiparent = (fint *)malloc(nboxes_v * sizeof(fint));
    tnchild = (fint *)malloc(nboxes_v * sizeof(fint));
    tichild = (fint *)malloc(4 * nboxes_v * sizeof(fint));
    tiflag = (fint *)malloc(nboxes_v * sizeof(fint));
    iboxtocurbox = (fint *)malloc(nboxes_v * sizeof(fint));
    tcenters = (double *)malloc(2 * nboxes_v * sizeof(double));

    for (ilev = 0; ilev <= nlevels_v; ilev++) {
        tladdr[LADDR(1, ilev)] = laddr[LADDR(1, ilev)];
        tladdr[LADDR(2, ilev)] = laddr[LADDR(2, ilev)];
    }
    tree_copy_(nboxes, centers, ilevel, iparent, nchild,
               ichild, tcenters, tilevel, tiparent,
               tnchild, tichild);

    for (ibox = 1; ibox <= nboxes_v; ibox++) {
        tiflag[ibox - 1] = iflag[ibox - 1];
    }

    /* Rearrange old arrays now */
    for (ilev = 0; ilev <= 1; ilev++) {
        for (ibox = laddr[LADDR(1, ilev)]; ibox <= laddr[LADDR(2, ilev)]; ibox++) {
            iboxtocurbox[ibox - 1] = ibox;
        }
    }

    ilevptr = (fint *)malloc((nlevels_v + 1) * sizeof(fint));
    ilevptr2 = (fint *)malloc(nlevels_v * sizeof(fint));

    /* ilevptr(2) = laddr(1,2) -- ilevptr is 1-indexed in Fortran, so C offset 1 */
    ilevptr[2 - 1] = laddr[LADDR(1, 2)];

    for (ilev = 2; ilev <= nlevels_v; ilev++) {
        nblev = laddr[LADDR(2, ilev)] - laddr[LADDR(1, ilev)] + 1;
        ilevptr2[ilev - 1] = ilevptr[ilev - 1] + nblev;
        nblev = laddrtail[LADDR(2, ilev)] - laddrtail[LADDR(1, ilev)] + 1;
        ilevptr[ilev + 1 - 1] = ilevptr2[ilev - 1] + nblev;
    }

    curbox = laddr[LADDR(1, 2)];
    for (ilev = 2; ilev <= nlevels_v; ilev++) {
        laddr[LADDR(1, ilev)] = curbox;
        for (ibox = tladdr[LADDR(1, ilev)]; ibox <= tladdr[LADDR(2, ilev)]; ibox++) {
            ilevel[curbox - 1] = tilevel[ibox - 1];
            nchild[curbox - 1] = tnchild[ibox - 1];
            centers[FA2(1, curbox, 2)] = tcenters[FA2(1, ibox, 2)];
            centers[FA2(2, curbox, 2)] = tcenters[FA2(2, ibox, 2)];
            iflag[curbox - 1] = tiflag[ibox - 1];
            iboxtocurbox[ibox - 1] = curbox;

            curbox = curbox + 1;
        }
        for (ibox = laddrtail[LADDR(1, ilev)]; ibox <= laddrtail[LADDR(2, ilev)]; ibox++) {
            ilevel[curbox - 1] = tilevel[ibox - 1];
            centers[FA2(1, curbox, 2)] = tcenters[FA2(1, ibox, 2)];
            centers[FA2(2, curbox, 2)] = tcenters[FA2(2, ibox, 2)];
            nchild[curbox - 1] = tnchild[ibox - 1];
            iflag[curbox - 1] = tiflag[ibox - 1];
            iboxtocurbox[ibox - 1] = curbox;

            curbox = curbox + 1;
        }
        laddr[LADDR(2, ilev)] = curbox - 1;
    }

    /* Handle parent/children using iboxtocurbox mapping */
    for (ibox = 1; ibox <= nboxes_v; ibox++) {
        if (tiparent[ibox - 1] == -1) iparent[iboxtocurbox[ibox - 1] - 1] = -1;
        if (tiparent[ibox - 1] > 0)
            iparent[iboxtocurbox[ibox - 1] - 1] = iboxtocurbox[tiparent[ibox - 1] - 1];
        for (i = 1; i <= 4; i++) {
            if (tichild[FA2(i, ibox, 4)] == -1)
                ichild[FA2(i, iboxtocurbox[ibox - 1], 4)] = -1;
            if (tichild[FA2(i, ibox, 4)] > 0)
                ichild[FA2(i, iboxtocurbox[ibox - 1], 4)] =
                    iboxtocurbox[tichild[FA2(i, ibox, 4)] - 1];
        }
    }

    free(tilevel);
    free(tiparent);
    free(tnchild);
    free(tichild);
    free(tiflag);
    free(iboxtocurbox);
    free(tcenters);
    free(ilevptr);
    free(ilevptr2);
    free(tladdr);
    return;
}


/* ---------------------------------------------------------------- */
/* pts_tree_sort - sort points into leaf boxes of an existing tree  */
/* ---------------------------------------------------------------- */
void FNAME(pts_tree_sort)(const fint *n, const double *xys,
                          const fint *itree, const fint *ltree,
                          const fint *nboxes, const fint *nlevels,
                          const fint *iptr, const double *centers,
                          fint *ixy, fint *ixyse)
{
    fint n_v = *n;
    fint nlevels_v = *nlevels;

    fint i, ilev, ibox;

    (void)ltree;

    for (i = 1; i <= n_v; i++) {
        ixy[i - 1] = i;
    }

    ixyse[FA2(1, 1, 2)] = 1;
    ixyse[FA2(2, 1, 2)] = n_v;

    for (ilev = 0; ilev <= nlevels_v - 1; ilev++) {
        for (ibox = itree[2 * ilev]; ibox <= itree[2 * ilev + 1]; ibox++) {
            /* itree(iptr(4)+ibox-1) > 0 */
            if (itree[iptr[3] + ibox - 1 - 1] > 0) {
                FNAME(sort_pts_to_children)(&ibox, nboxes, centers,
                                            &itree[iptr[4] - 1],
                                            xys, &n_v, ixy, ixyse);
            }
        }
    }

    return;
}
