/* arpack-ng\UTIL\iset.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/* ----------------------------------------------------------------------- */

/*     Only work with increment equal to 1 right now. */

int iset_(int32_t *n, int32_t *value, int32_t *array, 
	int32_t *inc)
{
    /* System generated locals */
    int32_t i__1;

    /* Local variables */
    int32_t i;

    /* Parameter adjustments */
    --array;

    /* Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	array[i] = *value;
    }

    return 0;
} /* iset_ */

