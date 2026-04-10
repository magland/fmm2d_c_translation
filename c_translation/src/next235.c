/*
 * next235.c - C translation of src/common/next235.f
 *
 * Returns the smallest integer >= base that is a product of
 * powers of 2, 3, and 5 (with at least one factor of 2).
 */

#include "next235.h"

fint FNAME(next235)(const double *base_p)
{
    double base_v = *base_p;
    fint result, numdiv;

    result = 2 * (fint)(base_v / 2.0 + 0.9999);
    if (result <= 0) result = 2;

    for (;;) {
        numdiv = result;
        while (numdiv / 2 * 2 == numdiv) {
            numdiv = numdiv / 2;
        }
        while (numdiv / 3 * 3 == numdiv) {
            numdiv = numdiv / 3;
        }
        while (numdiv / 5 * 5 == numdiv) {
            numdiv = numdiv / 5;
        }
        if (numdiv == 1) return result;
        result = result + 2;
    }
}
