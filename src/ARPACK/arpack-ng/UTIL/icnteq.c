/* arpack-ng\UTIL\icnteq.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/* ----------------------------------------------------------------------- */

/*     Count the number of elements equal to a specified integer value. */

int __icnteq_(int *n, int *array, int *value)
{
    /* System generated locals */
    int ret_val, i__1;

    /* Local variables */
    int i, k;

    /* Parameter adjustments */
    --array;

    /* Function Body */
    k = 0;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
        if (array[i] == *value)
        {
            ++k;
        }
    }
    ret_val = k;

    return ret_val;
} /* icnteq_ */

