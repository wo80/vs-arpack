/* D:\Projekte\ARPACK\arpack-ng\SRC\dsortr.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \BeginDoc
 *
 * \Name: dsortr
 *
 * \Description:
 *  Sort the array X1 in the order specified by WHICH and optionally
 *  applies the permutation to the array X2.
 *
 * \Usage:
 *  call dsortr
 *     ( WHICH, APPLY, N, X1, X2 )
 *
 * \Arguments
 *  WHICH   Character*2.  (Input)
 *          'LM' -> X1 is sorted into increasing order of magnitude.
 *          'SM' -> X1 is sorted into decreasing order of magnitude.
 *          'LA' -> X1 is sorted into increasing order of algebraic.
 *          'SA' -> X1 is sorted into decreasing order of algebraic.
 *
 *  APPLY   Logical.  (Input)
 *          APPLY = .TRUE.  -> apply the sorted order to X2.
 *          APPLY = .FALSE. -> do not apply the sorted order to X2.
 *
 *  N       Integer.  (INPUT)
 *          Size of the arrays.
 *
 *  X1      Double precision array of length N.  (INPUT/OUTPUT)
 *          The array to be sorted.
 *
 *  X2      Double precision array of length N.  (INPUT/OUTPUT)
 *          Only referenced if APPLY = .TRUE.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Author
 *     Danny Sorensen               Phuong Vu
 *     Richard Lehoucq              CRPC / Rice University
 *     Dept. of Computational &     Houston, Texas
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * \Revision history:
 *     12/16/93: Version ' 2.1'.
 *               Adapted from the sort routine in LANSO.
 *
 * \SCCS Information: @(#)
 * FILE: sortr.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2
 *
 * \EndLib
 */

int dsortr_(char *which, bool *apply, int32_t *n, double *x1, double *x2)
{
    /* System generated locals */
    int32_t i__1;
    double d__1, d__2;

    /* Builtin functions */
    

    /* Local variables */
    int32_t i__, j, igap;
    double temp;

     /* ---------------- */
     /* Scalar Arguments */
     /* ---------------- */

     /* --------------------- */
     /* Executable Statements */
     /* --------------------- */

    igap = *n / 2;

    if (strcmp(which, "SA") == 0) {

/*        X1 is sorted into decreasing order of algebraic. */

L10:
	if (igap == 0) {
	    goto L9000;
	}
	i__1 = *n - 1;
	for (i__ = igap; i__ <= i__1; ++i__) {
	    j = i__ - igap;
L20:

	    if (j < 0) {
		goto L30;
	    }

	    if (x1[j] < x1[j + igap]) {
		temp = x1[j];
		x1[j] = x1[j + igap];
		x1[j + igap] = temp;
		if (*apply) {
		    temp = x2[j];
		    x2[j] = x2[j + igap];
		    x2[j + igap] = temp;
		}
	    } else {
		goto L30;
	    }
	    j -= igap;
	    goto L20;
L30:
	    ;
	}
	igap /= 2;
	goto L10;

    } else if (strcmp(which, "SM") == 0) {

/*        X1 is sorted into decreasing order of magnitude. */

L40:
	if (igap == 0) {
	    goto L9000;
	}
	i__1 = *n - 1;
	for (i__ = igap; i__ <= i__1; ++i__) {
	    j = i__ - igap;
L50:

	    if (j < 0) {
		goto L60;
	    }

	    if ((d__1 = x1[j], abs(d__1)) < (d__2 = x1[j + igap], abs(d__2))) 
		    {
		temp = x1[j];
		x1[j] = x1[j + igap];
		x1[j + igap] = temp;
		if (*apply) {
		    temp = x2[j];
		    x2[j] = x2[j + igap];
		    x2[j + igap] = temp;
		}
	    } else {
		goto L60;
	    }
	    j -= igap;
	    goto L50;
L60:
	    ;
	}
	igap /= 2;
	goto L40;

    } else if (strcmp(which, "LA") == 0) {

/*        X1 is sorted into increasing order of algebraic. */

L70:
	if (igap == 0) {
	    goto L9000;
	}
	i__1 = *n - 1;
	for (i__ = igap; i__ <= i__1; ++i__) {
	    j = i__ - igap;
L80:

	    if (j < 0) {
		goto L90;
	    }

	    if (x1[j] > x1[j + igap]) {
		temp = x1[j];
		x1[j] = x1[j + igap];
		x1[j + igap] = temp;
		if (*apply) {
		    temp = x2[j];
		    x2[j] = x2[j + igap];
		    x2[j + igap] = temp;
		}
	    } else {
		goto L90;
	    }
	    j -= igap;
	    goto L80;
L90:
	    ;
	}
	igap /= 2;
	goto L70;

    } else if (strcmp(which, "LM") == 0) {

/*        X1 is sorted into increasing order of magnitude. */

L100:
	if (igap == 0) {
	    goto L9000;
	}
	i__1 = *n - 1;
	for (i__ = igap; i__ <= i__1; ++i__) {
	    j = i__ - igap;
L110:

	    if (j < 0) {
		goto L120;
	    }

	    if ((d__1 = x1[j], abs(d__1)) > (d__2 = x1[j + igap], abs(d__2))) 
		    {
		temp = x1[j];
		x1[j] = x1[j + igap];
		x1[j + igap] = temp;
		if (*apply) {
		    temp = x2[j];
		    x2[j] = x2[j + igap];
		    x2[j + igap] = temp;
		}
	    } else {
		goto L120;
	    }
	    j -= igap;
	    goto L110;
L120:
	    ;
	}
	igap /= 2;
	goto L100;
    }

L9000:
    return 0;

     /* ------------- */
     /* End of dsortr */
     /* ------------- */

} /* dsortr_ */

