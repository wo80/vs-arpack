/* D:\Projekte\csparse-interop\vs-arpack\src\ARPACK\arpack-ng\UTIL\iset.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"


/* ----------------------------------------------------------------------- */

/*     Only work with increment equal to 1 right now. */

/* Subroutine */ int iset_(integer *n, integer *value, integer *array, 
	integer *inc)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;



    /* Parameter adjustments */
    --array;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	array[i__] = *value;
/* L10: */
    }

    return 0;
} /* iset_ */

