/* D:\Projekte\ARPACK\arpack-ng\UTIL\icopy.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/* -------------------------------------------------------------------- */
/* \Documentation */

/* \Name: ICOPY */

/* \Description: */
/*     ICOPY copies an integer vector lx to an integer vector ly. */

/* \Usage: */
/*     call icopy ( n, lx, inc, ly, incy ) */

/* \Arguments: */
/*    n        integer (input) */
/*             On entry, n is the number of elements of lx to be */
/*             copied to ly. */

/*    lx       integer array (input) */
/*             On entry, lx is the integer vector to be copied. */

/*    incx     integer (input) */
/*             On entry, incx is the increment between elements of lx. */

/*    ly       integer array (input) */
/*             On exit, ly is the integer vector that contains the */
/*             copy of lx. */

/*    incy     integer (input) */
/*             On entry, incy is the increment between elements of ly. */

/* \Enddoc */

/* -------------------------------------------------------------------- */

/* Subroutine */ int icopy_(integer *n, integer *lx, integer *incx, integer *
	ly, integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, ix, iy;


/*     ---------------------------- */
/*     Specifications for arguments */
/*     ---------------------------- */

/*     ---------------------------------- */
/*     Specifications for local variables */
/*     ---------------------------------- */

/*     -------------------------- */
/*     First executable statement */
/*     -------------------------- */
    /* Parameter adjustments */
    --ly;
    --lx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/* .....code for unequal increments or equal increments */
/*     not equal to 1 */
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ly[iy] = lx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/* .....code for both increments equal to 1 */

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ly[i__] = lx[i__];
/* L30: */
    }
    return 0;
} /* icopy_ */

