#ifndef MBHKERNELS2D_H
#define MBHKERNELS2D_H

#include "fmm2d_c.h"

void FNAME(mbh2d_directcp_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *targ,
    const fint *nt,
    double *pot,
    const double *thresh);

void FNAME(mbh2d_directdp_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *dipstr,
    const double *dipvec,
    const double *targ,
    const fint *nt,
    double *pot,
    const double *thresh);

void FNAME(mbh2d_directcdp_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *dipstr,
    const double *dipvec,
    const double *targ,
    const fint *nt,
    double *pot,
    const double *thresh);

void FNAME(mbh2d_directqp_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *quadstr,
    const double *quadvec,
    const double *targ,
    const fint *nt,
    double *pot,
    const double *thresh);

void FNAME(mbh2d_directcqp_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *quadstr,
    const double *quadvec,
    const double *targ,
    const fint *nt,
    double *pot,
    const double *thresh);

void FNAME(mbh2d_directdqp_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *dipstr,
    const double *dipvec,
    const double *quadstr,
    const double *quadvec,
    const double *targ,
    const fint *nt,
    double *pot,
    const double *thresh);

void FNAME(mbh2d_directcdqp_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *dipstr,
    const double *dipvec,
    const double *quadstr,
    const double *quadvec,
    const double *targ,
    const fint *nt,
    double *pot,
    const double *thresh);

void FNAME(mbh2d_directop_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    const double *thresh);

void FNAME(mbh2d_directcop_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    const double *thresh);

void FNAME(mbh2d_directdop_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *dipstr,
    const double *dipvec,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    const double *thresh);

void FNAME(mbh2d_directcdop_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *dipstr,
    const double *dipvec,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    const double *thresh);

void FNAME(mbh2d_directqop_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *quadstr,
    const double *quadvec,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    const double *thresh);

void FNAME(mbh2d_directcqop_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *quadstr,
    const double *quadvec,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    const double *thresh);

void FNAME(mbh2d_directdqop_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *dipstr,
    const double *dipvec,
    const double *quadstr,
    const double *quadvec,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    const double *thresh);

void FNAME(mbh2d_directcdqop_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *dipstr,
    const double *dipvec,
    const double *quadstr,
    const double *quadvec,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    const double *thresh);

void FNAME(mbh2d_directcg_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    const double *thresh);

void FNAME(mbh2d_directdg_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *dipstr,
    const double *dipvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    const double *thresh);

void FNAME(mbh2d_directcdg_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *dipstr,
    const double *dipvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    const double *thresh);

void FNAME(mbh2d_directqg_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *quadstr,
    const double *quadvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    const double *thresh);

void FNAME(mbh2d_directcqg_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *quadstr,
    const double *quadvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    const double *thresh);

void FNAME(mbh2d_directdqg_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *dipstr,
    const double *dipvec,
    const double *quadstr,
    const double *quadvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    const double *thresh);

void FNAME(mbh2d_directcdqg_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *dipstr,
    const double *dipvec,
    const double *quadstr,
    const double *quadvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    const double *thresh);

void FNAME(mbh2d_directog_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    const double *thresh);

void FNAME(mbh2d_directcog_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    const double *thresh);

void FNAME(mbh2d_directdog_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *dipstr,
    const double *dipvec,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    const double *thresh);

void FNAME(mbh2d_directcdog_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *dipstr,
    const double *dipvec,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    const double *thresh);

void FNAME(mbh2d_directqog_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *quadstr,
    const double *quadvec,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    const double *thresh);

void FNAME(mbh2d_directcqog_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *quadstr,
    const double *quadvec,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    const double *thresh);

void FNAME(mbh2d_directdqog_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *dipstr,
    const double *dipvec,
    const double *quadstr,
    const double *quadvec,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    const double *thresh);

void FNAME(mbh2d_directcdqog_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *dipstr,
    const double *dipvec,
    const double *quadstr,
    const double *quadvec,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    const double *thresh);

void FNAME(mbh2d_directch_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    double *hess,
    const double *thresh);

void FNAME(mbh2d_directdh_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *dipstr,
    const double *dipvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    double *hess,
    const double *thresh);

void FNAME(mbh2d_directcdh_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *dipstr,
    const double *dipvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    double *hess,
    const double *thresh);

void FNAME(mbh2d_directqh_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *quadstr,
    const double *quadvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    double *hess,
    const double *thresh);

void FNAME(mbh2d_directcqh_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *quadstr,
    const double *quadvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    double *hess,
    const double *thresh);

void FNAME(mbh2d_directdqh_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *dipstr,
    const double *dipvec,
    const double *quadstr,
    const double *quadvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    double *hess,
    const double *thresh);

void FNAME(mbh2d_directcdqh_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *dipstr,
    const double *dipvec,
    const double *quadstr,
    const double *quadvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    double *hess,
    const double *thresh);

void FNAME(mbh2d_directoh_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    double *hess,
    const double *thresh);

void FNAME(mbh2d_directcoh_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    double *hess,
    const double *thresh);

void FNAME(mbh2d_directdoh_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *dipstr,
    const double *dipvec,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    double *hess,
    const double *thresh);

void FNAME(mbh2d_directcdoh_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *dipstr,
    const double *dipvec,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    double *hess,
    const double *thresh);

void FNAME(mbh2d_directqoh_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *quadstr,
    const double *quadvec,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    double *hess,
    const double *thresh);

void FNAME(mbh2d_directcqoh_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *quadstr,
    const double *quadvec,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    double *hess,
    const double *thresh);

void FNAME(mbh2d_directdqoh_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *dipstr,
    const double *dipvec,
    const double *quadstr,
    const double *quadvec,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    double *hess,
    const double *thresh);

void FNAME(mbh2d_directcdqoh_vec)(
    const fint *nd,
    const double *beta,
    const double *source,
    const fint *ns,
    const double *charge,
    const double *dipstr,
    const double *dipvec,
    const double *quadstr,
    const double *quadvec,
    const double *octstr,
    const double *octvec,
    const double *targ,
    const fint *nt,
    double *pot,
    double *grad,
    double *hess,
    const double *thresh);

void FNAME(mbh2d_directmpsp_vec)(const fint *nd, const double *beta,
    const double *source, const fint *ns, const fcomplex *mbhmps,
    const fcomplex *ymps, const fint *ntermsmps, const double *targ,
    const fint *nt, double *pot, const double *thresh);

void FNAME(mbh2d_directmpsg_vec)(const fint *nd, const double *beta,
    const double *source, const fint *ns, const fcomplex *mbhmps,
    const fcomplex *ymps, const fint *ntermsmps, const double *targ,
    const fint *nt, double *pot, double *grad, const double *thresh);

void FNAME(mbh2d_directmpsh_vec)(const fint *nd, const double *beta,
    const double *source, const fint *ns, const fcomplex *mbhmps,
    const fcomplex *ymps, const fint *ntermsmps, const double *targ,
    const fint *nt, double *pot, double *grad, double *hess,
    const double *thresh);

#endif /* MBHKERNELS2D_H */
