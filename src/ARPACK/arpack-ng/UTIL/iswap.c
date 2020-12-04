/* arpack-ng\UTIL\iswap.f -- translated by f2c (version 20100827). */

#include "arpack.h"

int iswap_(int32_t *n, int32_t *sx, int32_t *incx, int32_t *sy, int32_t *incy)
{
    /* System generated locals */
    int32_t i__1;

    /* Local variables */
    int32_t i, m, ix, iy, mp1, stemp;

    /*     interchanges two vectors. */
    /*     uses unrolled loops for increments equal to 1. */
    /*     jack dongarra, linpack, 3/11/78. */

    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0)
    {
        return 0;
    }
    if (*incx == 1 && *incy == 1)
    {
        goto L20;
    }

    /*       code for unequal increments or equal increments not equal */
    /*         to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0)
    {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0)
    {
        iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
        stemp = sx[ix];
        sx[ix] = sy[iy];
        sy[iy] = stemp;
        ix += *incx;
        iy += *incy;
    }
    return 0;

    /*       code for both increments equal to 1 */

    /*       clean-up loop */

L20:
    m = *n % 3;
    if (m == 0)
    {
        goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i)
    {
        stemp = sx[i];
        sx[i] = sy[i];
        sy[i] = stemp;
    }
    if (*n < 3)
    {
        return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 3)
    {
        stemp = sx[i];
        sx[i] = sy[i];
        sy[i] = stemp;
        stemp = sx[i + 1];
        sx[i + 1] = sy[i + 1];
        sy[i + 1] = stemp;
        stemp = sx[i + 2];
        sx[i + 2] = sy[i + 2];
        sy[i + 2] = stemp;
    }
    return 0;
} /* iswap_ */

