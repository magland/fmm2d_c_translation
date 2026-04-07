/*
 * cumsum.c - C translation of src/common/cumsum.f
 *
 * Original Fortran routines:
 *
 *   subroutine cumsum(n,a,b)        - top-level wrapper
 *   subroutine cumsum1(n,a,b)       - serial prefix sum
 *   subroutine cumsum_para(n,a,b,nd,d) - openmp parallel prefix sum
 *
 * The OpenMP pragmas (c$OMP, c$) in the Fortran source are comments;
 * this translation ignores them and emits sequential C. The serialized
 * version of cumsum_para still produces the correct cumulative sum
 * because, with a single thread, the per-thread offset accumulation
 * adds nothing.
 */

#include "cumsum.h"

void FNAME(cumsum)(const fint *n, const fint *a, fint *b)
{
    /* parameter (lend = 200) */
    const fint lend = 200;
    fint d[200];
    fint *d2;
    fint maxth;
    fint nd_arg;

    /* data nsmall / 10000 / */
    fint nsmall = 10000;

    /* if small problem, don't do parallel */
    if (*n < nsmall) goto L300;

    /* get upper bound of number of threads on hand */
    maxth = 1;
    /* c$    maxth = omp_get_max_threads()  -- ignored */

    /* no benefit to parallelize if only 2 threads */
    if (maxth <= 2) goto L300;

    /* if not a ton of processors, use d on stack */
    if (maxth <= lend) goto L200;

    /* if tons of processors, allocate d2 */
    d2 = (fint *)malloc(maxth * sizeof(fint));
    FNAME(cumsum_para)(n, a, b, &maxth, d2);
    free(d2);
    return;

L200:
    nd_arg = lend;
    FNAME(cumsum_para)(n, a, b, &nd_arg, d);
    return;

L300:
    FNAME(cumsum1)(n, a, b);
    return;
}

void FNAME(cumsum1)(const fint *n, const fint *a, fint *b)
{
    fint i, isum;

    isum = 0;

    for (i = 1; i <= *n; i++) {
        isum = isum + a[i - 1];
        b[i - 1] = isum;
    }

    return;
}

void FNAME(cumsum_para)(const fint *n, const fint *a, fint *b,
                        const fint *nd, fint *d)
{
    fint i, isum, offset;
    fint id;

    (void)nd;

    /* c$OMP PARALLEL DEFAULT(none) SHARED(a,b,n,d) ... -- ignored */

    id = 1;
    /* c$    id = omp_get_thread_num()+1  -- ignored */

    /* compute cumulative sums of portions (parallel) */

    isum = 0;
    /* c$OMP DO SCHEDULE(static) */
    for (i = 1; i <= *n; i++) {
        isum = isum + a[i - 1];
        b[i - 1] = isum;
    }
    /* c$OMP END DO no wait */
    d[id - 1] = isum;

    /* c$OMP BARRIER */

    /* accumulate the ends of the partial sums (each thread) */

    offset = 0;
    for (i = 1; i <= id - 1; i++) {
        offset = offset + d[i - 1];
    }

    /* c$OMP DO SCHEDULE(static) */
    for (i = 1; i <= *n; i++) {
        b[i - 1] = b[i - 1] + offset;
    }
    /* c$OMP END DO no wait */

    /* c$OMP END PARALLEL */

    return;
}
