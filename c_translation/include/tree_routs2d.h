/*
 * tree_routs2d.h - C translation of src/common/tree_routs2d.f
 *
 * Common routines for generating and processing a level-restricted
 * quad tree in 2D for the FMM. These manage integer-index bookkeeping
 * (parent / child / neighbor / list arrays) on top of an explicit
 * tree representation.
 *
 * print_tree (file I/O) is intentionally not translated.
 */

#ifndef FMM2D_TREE_ROUTS2D_H
#define FMM2D_TREE_ROUTS2D_H

#include "fmm2d_c.h"

void FNAME(tree_refine_boxes)(const fint *irefinebox, const fint *nboxes,
                              const fint *ifirstbox, const fint *nbloc,
                              double *centers, const double *bs,
                              fint *nbctr, const fint *nlctr,
                              fint *ilevel, fint *iparent,
                              fint *nchild, fint *ichild);

void FNAME(tree_copy)(const fint *nb,
                      const double *centers, const fint *ilevel,
                      const fint *iparent, const fint *nchild,
                      const fint *ichild,
                      double *centers2, fint *ilevel2,
                      fint *iparent2, fint *nchild2,
                      fint *ichild2);

void FNAME(computecoll)(const fint *nlevels, const fint *nboxes,
                        const fint *laddr, const double *boxsize,
                        const double *centers, const fint *iparent,
                        const fint *nchild, const fint *ichild,
                        const fint *iper,
                        fint *nnbors, fint *nbors);

void FNAME(updateflags)(const fint *curlev, const fint *nboxes,
                        const fint *nlevels, const fint *laddr,
                        const fint *nchild, const fint *ichild,
                        const fint *nnbors, const fint *nbors,
                        const double *centers, const double *boxsize,
                        fint *iflag);

void FNAME(tree_refine_boxes_flag)(fint *iflag, const fint *nboxes,
                                   const fint *ifirstbox, const fint *nbloc,
                                   double *centers, const double *bs,
                                   fint *nbctr, const fint *nlctr,
                                   fint *ilevel, fint *iparent,
                                   fint *nchild, fint *ichild);

void FNAME(computemnlists)(const fint *nlevels, const fint *nboxes,
                           const fint *itree, const fint *ltree,
                           const fint *iptr, const double *centers,
                           const double *boxsize, const fint *iper,
                           fint *mnlist1, fint *mnlist2,
                           fint *mnlist3, fint *mnlist4);

void FNAME(computelists)(const fint *nlevels, const fint *nboxes,
                         const fint *itree, const fint *ltree,
                         const fint *iptr, const double *centers,
                         const double *boxsize, const fint *iper,
                         const fint *mnlist1, fint *nlist1, fint *list1,
                         const fint *mnlist2, fint *nlist2, fint *list2,
                         const fint *mnlist3, fint *nlist3, fint *list3,
                         const fint *mnlist4, fint *nlist4, fint *list4);

#endif /* FMM2D_TREE_ROUTS2D_H */
