/*
 * dfft_threadsafe.h - C translation of src/common/dfft_threadsafe.f
 *
 * Complex FFT routines from FFTPACK (thread-safe variant).
 * Only the complex-FFT subset (zffti/zfftf/zfftb and their
 * helpers) is translated here — the real FFT and trig-transform
 * routines are not needed by fmm2d.
 */

#ifndef FMM2D_DFFT_THREADSAFE_H
#define FMM2D_DFFT_THREADSAFE_H

#include "fmm2d_c.h"

/* Top-level entry points */
void FNAME(zffti)(const fint *n, double *wsave);
void FNAME(zfftf)(const fint *n, double *c, double *wsave);
void FNAME(zfftb)(const fint *n, double *c, double *wsave);

/* Initialization helper */
void FNAME(zffti1)(const fint *n, double *wa, fint *ifac);

/* Forward/backward dispatch helpers */
void FNAME(zfftf1)(const fint *n, double *c, double *ch,
                   double *wa, fint *ifac);
void FNAME(zfftb1)(const fint *n, double *c, double *ch,
                   double *wa, fint *ifac);

/* Backward butterfly passes */
void FNAME(dpassb)(fint *nac, const fint *ido, const fint *ip,
                   const fint *l1, const fint *idl1,
                   double *cc, double *c1, double *c2,
                   double *ch, double *ch2, const double *wa);
void FNAME(dpassb2)(const fint *ido, const fint *l1,
                    const double *cc, double *ch, const double *wa1);
void FNAME(dpassb3)(const fint *ido, const fint *l1,
                    const double *cc, double *ch,
                    const double *wa1, const double *wa2);
void FNAME(dpassb4)(const fint *ido, const fint *l1,
                    const double *cc, double *ch,
                    const double *wa1, const double *wa2,
                    const double *wa3);
void FNAME(dpassb5)(const fint *ido, const fint *l1,
                    const double *cc, double *ch,
                    const double *wa1, const double *wa2,
                    const double *wa3, const double *wa4);

/* Forward butterfly passes */
void FNAME(dpassf)(fint *nac, const fint *ido, const fint *ip,
                   const fint *l1, const fint *idl1,
                   double *cc, double *c1, double *c2,
                   double *ch, double *ch2, const double *wa);
void FNAME(dpassf2)(const fint *ido, const fint *l1,
                    const double *cc, double *ch, const double *wa1);
void FNAME(dpassf3)(const fint *ido, const fint *l1,
                    const double *cc, double *ch,
                    const double *wa1, const double *wa2);
void FNAME(dpassf4)(const fint *ido, const fint *l1,
                    const double *cc, double *ch,
                    const double *wa1, const double *wa2,
                    const double *wa3);
void FNAME(dpassf5)(const fint *ido, const fint *l1,
                    const double *cc, double *ch,
                    const double *wa1, const double *wa2,
                    const double *wa3, const double *wa4);

#endif
