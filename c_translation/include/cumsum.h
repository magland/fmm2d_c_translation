/*
 * cumsum.h - C translation of src/common/cumsum.f
 */

#ifndef FMM2D_CUMSUM_H
#define FMM2D_CUMSUM_H

#include "fmm2d_c.h"

/*
 * void cumsum(n, a, b)
 *
 * Computes the cumulative (prefix) sum of an integer array a of length n
 * and stores the result in b. This is the wrapper that decides between
 * the serial code (cumsum1) and the OpenMP parallel code (cumsum_para).
 * In this C translation OpenMP pragmas are ignored, so cumsum_para runs
 * sequentially but produces the same numerical result.
 */
void FNAME(cumsum)(const fint *n, const fint *a, fint *b);

/*
 * void cumsum1(n, a, b)
 *
 * Serial prefix sum: b(i) = sum_{j<=i} a(j) for i = 1,...,n.
 */
void FNAME(cumsum1)(const fint *n, const fint *a, fint *b);

/*
 * void cumsum_para(n, a, b, nd, d)
 *
 * OpenMP parallel prefix sum (sequentialized in this C translation).
 * d is a work array of length nd; nd should be at least the number of
 * threads to be used.
 */
void FNAME(cumsum_para)(const fint *n, const fint *a, fint *b,
                        const fint *nd, fint *d);

#endif /* FMM2D_CUMSUM_H */
