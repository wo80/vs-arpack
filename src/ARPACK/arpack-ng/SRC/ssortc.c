/* D:\Projekte\ARPACK\arpack-ng\SRC\ssortc.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \BeginDoc
 *
 * \Name: ssortc
 *
 * \Description:
 *  Sorts the complex array in XREAL and XIMAG into the order
 *  specified by WHICH and optionally applies the permutation to the
 *  float array Y. It is assumed that if an element of XIMAG is
 *  nonzero, then its negative is also an element. In other words,
 *  both members of a complex conjugate pair are to be sorted and the
 *  pairs are kept adjacent to each other.
 *
 * \Usage:
 *  call ssortc
 *     ( WHICH, APPLY, N, XREAL, XIMAG, Y )
 *
 * \Arguments
 *  WHICH   Character*2.  (Input)
 *          'LM' -> sort XREAL,XIMAG into increasing order of magnitude.
 *          'SM' -> sort XREAL,XIMAG into decreasing order of magnitude.
 *          'LR' -> sort XREAL into increasing order of algebraic.
 *          'SR' -> sort XREAL into decreasing order of algebraic.
 *          'LI' -> sort XIMAG into increasing order of magnitude.
 *          'SI' -> sort XIMAG into decreasing order of magnitude.
 *          NOTE: If an element of XIMAG is non-zero, then its negative
 *                is also an element.
 *
 *  APPLY   Logical.  (Input)
 *          APPLY = .TRUE.  -> apply the sorted order to array Y.
 *          APPLY = .FALSE. -> do not apply the sorted order to array Y.
 *
 *  N       Integer.  (INPUT)
 *          Size of the arrays.
 *
 *  XREAL,  Real array of length N.  (INPUT/OUTPUT)
 *  XIMAG   Real and imaginary part of the array to be sorted.
 *
 *  Y       Real array of length N.  (INPUT/OUTPUT)
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
 *     xx/xx/92: Version ' 2.1'
 *               Adapted from the sort routine in LANSO.
 *
 * \SCCS Information: @(#)
 * FILE: sortc.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
 *
 * \EndLib
 */

int ssortc_(char *which, bool *apply, int32_t *n, float *xfloat, float *ximag, float *y)
{
    /* System generated locals */
    int32_t i__1;
    float r__1, r__2;

    /* Builtin functions */
    
    /* Local variables */
    int32_t i, j, igap;
    float temp, temp1, temp2;

     /* ---------------- */
     /* Scalar Arguments */
     /* ---------------- */

     /* --------------------- */
     /* Executable Statements */
     /* --------------------- */

    igap = *n / 2;

    if (strcmp(which, "LM") == 0) {

        /* ---------------------------------------------------- */
        /* Sort XREAL,XIMAG into increasing order of magnitude. */
        /* ---------------------------------------------------- */

L10:
	if (igap == 0) {
	    goto L9000;
	}

	i__1 = *n - 1;
	for (i = igap; i <= i__1; ++i) {
	    j = i - igap;
L20:

	    if (j < 0) {
		goto L30;
	    }

	    temp1 = slapy2_(&xfloat[j], &ximag[j]);
	    temp2 = slapy2_(&xfloat[j + igap], &ximag[j + igap]);

	    if (temp1 > temp2) {
		temp = xfloat[j];
		xfloat[j] = xfloat[j + igap];
		xfloat[j + igap] = temp;

		temp = ximag[j];
		ximag[j] = ximag[j + igap];
		ximag[j + igap] = temp;

		if (*apply) {
		    temp = y[j];
		    y[j] = y[j + igap];
		    y[j + igap] = temp;
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

        /* ---------------------------------------------------- */
        /* Sort XREAL,XIMAG into decreasing order of magnitude. */
        /* ---------------------------------------------------- */

L40:
	if (igap == 0) {
	    goto L9000;
	}

	i__1 = *n - 1;
	for (i = igap; i <= i__1; ++i) {
	    j = i - igap;
L50:

	    if (j < 0) {
		goto L60;
	    }

	    temp1 = slapy2_(&xfloat[j], &ximag[j]);
	    temp2 = slapy2_(&xfloat[j + igap], &ximag[j + igap]);

	    if (temp1 < temp2) {
		temp = xfloat[j];
		xfloat[j] = xfloat[j + igap];
		xfloat[j + igap] = temp;

		temp = ximag[j];
		ximag[j] = ximag[j + igap];
		ximag[j + igap] = temp;

		if (*apply) {
		    temp = y[j];
		    y[j] = y[j + igap];
		    y[j + igap] = temp;
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

    } else if (strcmp(which, "LR") == 0) {

        /* ---------------------------------------------- */
        /* Sort XREAL into increasing order of algebraic. */
        /* ---------------------------------------------- */

L70:
	if (igap == 0) {
	    goto L9000;
	}

	i__1 = *n - 1;
	for (i = igap; i <= i__1; ++i) {
	    j = i - igap;
L80:

	    if (j < 0) {
		goto L90;
	    }

	    if (xfloat[j] > xfloat[j + igap]) {
		temp = xfloat[j];
		xfloat[j] = xfloat[j + igap];
		xfloat[j + igap] = temp;

		temp = ximag[j];
		ximag[j] = ximag[j + igap];
		ximag[j + igap] = temp;

		if (*apply) {
		    temp = y[j];
		    y[j] = y[j + igap];
		    y[j + igap] = temp;
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

    } else if (strcmp(which, "SR") == 0) {

        /* ---------------------------------------------- */
        /* Sort XREAL into decreasing order of algebraic. */
        /* ---------------------------------------------- */

L100:
	if (igap == 0) {
	    goto L9000;
	}
	i__1 = *n - 1;
	for (i = igap; i <= i__1; ++i) {
	    j = i - igap;
L110:

	    if (j < 0) {
		goto L120;
	    }

	    if (xfloat[j] < xfloat[j + igap]) {
		temp = xfloat[j];
		xfloat[j] = xfloat[j + igap];
		xfloat[j + igap] = temp;

		temp = ximag[j];
		ximag[j] = ximag[j + igap];
		ximag[j + igap] = temp;

		if (*apply) {
		    temp = y[j];
		    y[j] = y[j + igap];
		    y[j + igap] = temp;
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

    } else if (strcmp(which, "LI") == 0) {

        /* ---------------------------------------------- */
        /* Sort XIMAG into increasing order of magnitude. */
        /* ---------------------------------------------- */

L130:
	if (igap == 0) {
	    goto L9000;
	}
	i__1 = *n - 1;
	for (i = igap; i <= i__1; ++i) {
	    j = i - igap;
L140:

	    if (j < 0) {
		goto L150;
	    }

	    if ((r__1 = ximag[j], dabs(r__1)) > (r__2 = ximag[j + igap], dabs(
		    r__2))) {
		temp = xfloat[j];
		xfloat[j] = xfloat[j + igap];
		xfloat[j + igap] = temp;

		temp = ximag[j];
		ximag[j] = ximag[j + igap];
		ximag[j + igap] = temp;

		if (*apply) {
		    temp = y[j];
		    y[j] = y[j + igap];
		    y[j + igap] = temp;
		}
	    } else {
		goto L150;
	    }
	    j -= igap;
	    goto L140;
L150:
	    ;
	}
	igap /= 2;
	goto L130;

    } else if (strcmp(which, "SI") == 0) {

        /* ---------------------------------------------- */
        /* Sort XIMAG into decreasing order of magnitude. */
        /* ---------------------------------------------- */

L160:
	if (igap == 0) {
	    goto L9000;
	}
	i__1 = *n - 1;
	for (i = igap; i <= i__1; ++i) {
	    j = i - igap;
L170:

	    if (j < 0) {
		goto L180;
	    }

	    if ((r__1 = ximag[j], dabs(r__1)) < (r__2 = ximag[j + igap], dabs(
		    r__2))) {
		temp = xfloat[j];
		xfloat[j] = xfloat[j + igap];
		xfloat[j + igap] = temp;

		temp = ximag[j];
		ximag[j] = ximag[j + igap];
		ximag[j + igap] = temp;

		if (*apply) {
		    temp = y[j];
		    y[j] = y[j + igap];
		    y[j + igap] = temp;
		}
	    } else {
		goto L180;
	    }
	    j -= igap;
	    goto L170;
L180:
	    ;
	}
	igap /= 2;
	goto L160;
    }

L9000:
    return 0;

     /* ------------- */
     /* End of ssortc */
     /* ------------- */

} /* ssortc_ */

