/* arpack-ng\UTIL\iset.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/* ----------------------------------------------------------------------- */

/*     Only work with increment equal to 1 right now. */

int __iset_(int *n, int *value, int *array, int *inc)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    int i;

    /* Parameter adjustments */
    --array;

    /* Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
        array[i] = *value;
    }

    return 0;
} /* iset_ */

