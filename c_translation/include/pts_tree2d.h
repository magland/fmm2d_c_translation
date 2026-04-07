/*
 * pts_tree2d.h - C translation of src/common/pts_tree2d.f
 *
 * Routines for building a level-restricted quad-tree from a set of
 * source / target points in 2D, used by the FMM. All 6 Fortran
 * subroutines of pts_tree2d.f are translated here:
 *
 *   pts_tree_mem           - size the tree (output nlevels, nboxes, ltree)
 *   pts_tree_build         - build the tree (output itree, iptr, centers, boxsize)
 *   sort_pts_to_children   - partition the points of one box into its 4 children
 *   pts_tree_fix_lr        - convert adaptive tree to level-restricted tree
 *   pts_tree_reorg         - reorganize tree arrays after fix_lr
 *   pts_tree_sort          - sort points into the leaf boxes of an existing tree
 *
 * The cross-file calls (tree_copy, tree_refine_boxes, computecoll,
 * tree_refine_boxes_flag, updateflags) are issued with bare Fortran
 * symbol names so the diff test can isolate this file.
 */

#ifndef FMM2D_PTS_TREE2D_H
#define FMM2D_PTS_TREE2D_H

#include "fmm2d_c.h"

void FNAME(pts_tree_mem)(const double *src, const fint *ns,
                         const double *targ, const fint *nt,
                         const fint *idivflag, const fint *ndiv,
                         const fint *nlmin, const fint *nlmax,
                         const fint *ifunif, const fint *iper,
                         fint *nlevels, fint *nboxes, fint *ltree);

void FNAME(pts_tree_build)(const double *src, const fint *ns,
                           const double *targ, const fint *nt,
                           const fint *idivflag, const fint *ndiv,
                           const fint *nlmin, const fint *nlmax,
                           const fint *ifunif, const fint *iper,
                           fint *nlevels, const fint *nboxes,
                           const fint *ltree, fint *itree, fint *iptr,
                           double *centers, double *boxsize);

void FNAME(sort_pts_to_children)(const fint *ibox, const fint *nboxes,
                                 const double *centers, const fint *ichild,
                                 const double *src, const fint *ns,
                                 fint *isrc, fint *isrcse);

void FNAME(pts_tree_fix_lr)(double *centers, const fint *nlevels,
                            fint *nboxes, const double *boxsize,
                            const fint *nbmax, const fint *nlmax,
                            const fint *iper, fint *laddr, fint *ilevel,
                            fint *iparent, fint *nchild, fint *ichild,
                            fint *nnbors, fint *nbors);

void FNAME(pts_tree_reorg)(const fint *nboxes, double *centers,
                           const fint *nlevels, fint *laddr,
                           const fint *laddrtail, fint *ilevel,
                           fint *iparent, fint *nchild, fint *ichild,
                           fint *iflag);

void FNAME(pts_tree_sort)(const fint *n, const double *xys,
                          const fint *itree, const fint *ltree,
                          const fint *nboxes, const fint *nlevels,
                          const fint *iptr, const double *centers,
                          fint *ixy, fint *ixyse);

#endif /* FMM2D_PTS_TREE2D_H */
