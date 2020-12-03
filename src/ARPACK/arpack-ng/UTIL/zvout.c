/* arpack-ng\UTIL\zvout.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/* ----------------------------------------------------------------------- */

/* \SCCS Information: @(#) */
/* FILE: zvout.f   SID: 2.1   DATE OF SID: 11/16/95   RELEASE: 2 */

/* ----------------------------------------------------------------------- */
/*  Routine:    ZVOUT */

/*  Purpose:    Complex*16 vector output routine. */

/*  Usage:      CALL ZVOUT (LOUT, N, CX, IDIGIT, IFMT) */

/*  Arguments */
/*     N      - Length of array CX.  (Input) */
/*     CX     - Complex*16 array to be printed.  (Input) */
/*     IFMT   - Format to be used in printing array CX.  (Input) */
/*     IDIGIT - Print up to IABS(IDIGIT) decimal digits per number.  (In) */
/*              If IDIGIT .LT. 0, printing is done with 72 columns. */
/*              If IDIGIT .GT. 0, printing is done with 132 columns. */

/* ----------------------------------------------------------------------- */

int zvout_(int32_t *n, zomplex *cx, 
	int32_t *idigit, char *ifmt)
{

    /* System generated locals */
    int32_t i__1, i__2, i__3;

    /* Local variables */
    int32_t i, k1, k2, lll;
    char line[80];
    int32_t ndigit;

    /* Parameter adjustments */
    --cx;

    /* Function Body */
/* Computing MIN */
    i__1 = strlen(ifmt);
    lll = min(i__1,80);
    i__1 = lll;
    for (i = 1; i <= i__1; ++i) {
	*&line[i - 1] = '-';
    }

    for (i = lll + 1; i <= 80; ++i) {
	*&line[i - 1] = ' ';
    }

    printf("\n %s\n %s", ifmt, line);

    if (*n <= 0) {
	return 0;
    }
    ndigit = *idigit;
    if (*idigit == 0) {
	ndigit = 4;
    }

/* ======================================================================= */
/*             CODE FOR OUTPUT USING 72 COLUMNS FORMAT */
/* ======================================================================= */

    if (*idigit < 0) {
	ndigit = -(*idigit);
	if (ndigit <= 4) {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 2) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 1;
		k2 = min(i__2,i__3);
	    printf("\n  %4d -   %4d:  ", k1, k2);
	    for (i = k1; i <= k2; ++i) {
		printf("(  %10.3f   %10.3f)  ", cx[i].r, cx[i].i);
	    }
	    }
	} else if (ndigit <= 6) {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 2) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 1;
		k2 = min(i__2,i__3);
	    printf("\n  %4d -   %4d:  ", k1, k2);
	    for (i = k1; i <= k2; ++i) {
		printf("(  %12.5f   %12.5f)  ", cx[i].r, cx[i].i);
	    }
	    }
	} else if (ndigit <= 8) {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 2) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 1;
		k2 = min(i__2,i__3);
	    printf("\n  %4d -   %4d:  ", k1, k2);
	    for (i = k1; i <= k2; ++i) {
		printf("(  %14.7f   %14.7f)  ", cx[i].r, cx[i].i);
	    }
	    }
	} else {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; ++k1) {
		printf("  %4d:  (  %20.13f   %20.13f)  ", k1, cx[i].r, cx[i].i);
	    }
	}

/* ======================================================================= */
/*             CODE FOR OUTPUT USING 132 COLUMNS FORMAT */
/* ======================================================================= */

    } else {
	if (ndigit <= 4) {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 4) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 3;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d:  ", k1, k2);
		for (i = k1; i <= k2; ++i) {
		printf("(  %10.3f   %10.3f)  ", cx[i].r, cx[i].i);
		}
	    }
	} else if (ndigit <= 6) {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 3) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 2;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d:  ", k1, k2);
		for (i = k1; i <= k2; ++i) {
		printf("(  %12.5f   %12.5f)  ", cx[i].r, cx[i].i);
		}
	    }
	} else if (ndigit <= 8) {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 3) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 2;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d:  ", k1, k2);
		for (i = k1; i <= k2; ++i) {
		printf("(  %14.7f   %14.7f)  ", cx[i].r, cx[i].i);
		}
	    }
	} else {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 2) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 1;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d:  ", k1, k2);
		for (i = k1; i <= k2; ++i) {
		printf("(  %20.13f   %20.13f)  ", cx[i].r, cx[i].i);
		}
	    }
	}
    }
    printf("\n");
    return 0;

} /* zvout_ */

