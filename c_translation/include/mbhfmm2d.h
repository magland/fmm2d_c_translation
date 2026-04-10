#ifndef MBHFMM2D_H
#define MBHFMM2D_H

#include "fmm2d_c.h"

void FNAME(mbhfmm2d)(const fint *nd, const double *eps, const double *beta,
    const fint *ns, const double *sources,
    const fint *ifcharge, const double *charge,
    const fint *ifdipole, const double *dipstr, const double *dipvec,
    const fint *ifquadpole, const double *quadstr, const double *quadvec,
    const fint *ifoctpole, const double *octstr, const double *octvec,
    const fint *iper, const fint *ifpgh,
    double *pot, double *grad, double *hess,
    const fint *nt, const double *targ, const fint *ifpghtarg,
    double *pottarg, double *gradtarg, double *hesstarg,
    fint *ier);

void FNAME(mbhfmm2dmain)(const fint *nd, const double *eps,
    const double *beta, const fint *nsource, const double *source,
    const fint *ntermsmps_p, const fcomplex *mbhmps, const fcomplex *ymps,
    const fint *ntarget, const double *target,
    const fint *iaddr, double *rmlexp,
    const double *carray, const fint *ldc_p,
    const fint *itree, const fint *ltree,
    const fint *iptr, const fint *ndiv, const fint *nlevels,
    const fint *nboxes, const fint *iper,
    const double *boxsize, const double *rscales,
    const double *centers, const fint *laddr,
    const fint *isrcse, const fint *itargse,
    const fint *nterms,
    const fint *ifpgh, double *pot, double *grad, double *hess,
    const fint *ifpghtarg, double *pottarg, double *gradtarg,
    double *hesstarg, fint *ier);

void FNAME(mbhfmm2dmps_direct_vec)(const fint *nd,
    const fint *istart, const fint *iend,
    const fint *jstart, const fint *jend,
    const double *beta, const double *source,
    const fcomplex *mbhmps, const fcomplex *ymps,
    const fint *ntermsmps,
    const double *targ, const fint *ifpgh,
    double *pot, double *grad, double *hess,
    const double *thresh);

void FNAME(mbh2dmpalloc)(const fint *nd, const fint *laddr, fint *iaddr,
    const fint *nlevels, fint *lmptot, const fint *nterms);

void FNAME(mbh2dmpzero)(const fint *nd, fcomplex *mpole, const fint *nterms);

#endif /* MBHFMM2D_H */
