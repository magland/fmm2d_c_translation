#ifndef MBHROUTS2D_H
#define MBHROUTS2D_H

#include "fmm2d_c.h"

void FNAME(mbh2d_rk)(double *pow, double *dpow,
    const double *r, const double *beta, const double *rscale,
    const fint *nterms);
void FNAME(mbh2d_rksc)(double *pow, double *dpow,
    const double *r, const double *beta, const double *rscale,
    const fint *nterms);
void FNAME(mbh2d_rmk)(double *pow, double *dpow,
    const double *r, const double *beta, const double *rscale,
    const fint *nterms);
void FNAME(mbh2d_init_carray)(double *carray, const fint *ldc);

void FNAME(mbh2dconvtomp_vec)(
    const fint *nd, const double *beta, const fint *ns,
    const fint *ifcharge, const double *charge,
    const fint *ifdipole, const double *dipstr, const double *dipvec,
    const fint *ifquad, const double *quadstr, const double *quadvec,
    const fint *ifoct, const double *octstr, const double *octvec,
    const fint *nterms, fcomplex *mbhmpole, fcomplex *ympole);

void FNAME(mbh2dmpevalp_vec)(const fint *nd, const double *beta,
    const double *rscale, const double *center,
    const fcomplex *mbhmpole, const fcomplex *ympole,
    const fint *nterms, const double *ztarg, const fint *ntarg,
    double *pot);
void FNAME(mbh2dmpevalg_vec)(const fint *nd, const double *beta,
    const double *rscale, const double *center,
    const fcomplex *mbhmpole, const fcomplex *ympole,
    const fint *nterms, const double *ztarg, const fint *ntarg,
    double *pot, double *grad);
void FNAME(mbh2dmpevalh_vec)(const fint *nd, const double *beta,
    const double *rscale, const double *center,
    const fcomplex *mbhmpole, const fcomplex *ympole,
    const fint *nterms, const double *ztarg, const fint *ntarg,
    double *pot, double *grad, double *hess);

void FNAME(mbh2dtaevalp_vec)(const fint *nd, const double *beta,
    const double *rscale, const double *center,
    const fcomplex *mbhloc, const fcomplex *lloc,
    const fint *nterms, const double *ztarg, const fint *ntarg,
    double *pot);
void FNAME(mbh2dtaevalg_vec)(const fint *nd, const double *beta,
    const double *rscale, const double *center,
    const fcomplex *mbhloc, const fcomplex *lloc,
    const fint *nterms, const double *ztarg, const fint *ntarg,
    double *pot, double *grad);
void FNAME(mbh2dtaevalh_vec)(const fint *nd, const double *beta,
    const double *rscale, const double *center,
    const fcomplex *mbhloc, const fcomplex *lloc,
    const fint *nterms, const double *ztarg, const fint *ntarg,
    double *pot, double *grad, double *hess);

void FNAME(mbh2dformmpmp_vec)(const fint *nd, const double *beta,
    const double *rscale, const double *source, const fint *ns,
    fcomplex *mbhmpolesrc, fcomplex *ympolesrc,
    const fint *ntermsrc, const double *center,
    const fint *nterms, fcomplex *mbhmpole, fcomplex *ympole);
void FNAME(mbh2dformtamp_vec)(const fint *nd, const double *beta,
    const double *rscale, const double *source, const fint *ns,
    fcomplex *mbhmpolesrc, fcomplex *ympolesrc,
    const fint *ntermsrc, const double *center,
    const fint *nterms, fcomplex *mbhloc, fcomplex *lloc);

void FNAME(mbh2dmpmp_vec)(const fint *nd, const double *beta,
    const double *rscale1, const double *center1,
    const fcomplex *mbhmpole1, const fcomplex *ympole1, const fint *nterms1,
    const double *rscale2, const double *center2,
    fcomplex *mbhmpole2, fcomplex *ympole2, const fint *nterms2);
void FNAME(mbh2dmploc_vec)(const fint *nd, const double *beta,
    const double *rscale1, const double *center1,
    const fcomplex *mbhmpole, const fcomplex *ympole, const fint *nterms1,
    const double *rscale2, const double *center2,
    fcomplex *mbhloc, fcomplex *lloc, const fint *nterms2);
void FNAME(mbh2dlocloc_vec)(const fint *nd, const double *beta,
    const double *rscale1, const double *center1,
    const fcomplex *mbhloc1, const fcomplex *lloc1, const fint *nterms1,
    const double *rscale2, const double *center2,
    fcomplex *mbhloc2, fcomplex *lloc2, const fint *nterms2,
    const double *carray, const fint *ldc);

#endif /* MBHROUTS2D_H */
