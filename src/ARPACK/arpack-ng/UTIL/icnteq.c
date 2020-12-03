/* D:\Projekte\ARPACK\arpack-ng\UTIL\icnteq.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/* ----------------------------------------------------------------------- */

/*     Count the number of elements equal to a specified integer value. */

int32_t icnteq_(int32_t *n, int32_t *array, int32_t *value)
{
    /* System generated locals */
    int32_t ret_val, i__1;

    /* Local variables */
    int32_t i__, k;

    /* Parameter adjustments */
    --array;

    /* Function Body */
    k = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (array[i__] == *value) {
	    ++k;
	}
/* L10: */
    }
    ret_val = k;

    return ret_val;
} /* icnteq_ */

