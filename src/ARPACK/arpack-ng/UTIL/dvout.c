/* arpack-ng\UTIL\dvout.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/* ----------------------------------------------------------------------- */
/*  Routine:    DVOUT */

/*  Purpose:    Real vector output routine. */

/*  Usage:      CALL DVOUT (LOUT, N, SX, IDIGIT, IFMT) */

/*  Arguments */
/*     N      - Length of array SX.  (Input) */
/*     SX     - Real array to be printed.  (Input) */
/*     IFMT   - Format to be used in printing array SX.  (Input) */
/*     IDIGIT - Print up to IABS(IDIGIT) decimal digits per number.  (In) */
/*              If IDIGIT .LT. 0, printing is done with 72 columns. */
/*              If IDIGIT .GT. 0, printing is done with 132 columns. */

/* ----------------------------------------------------------------------- */

int dvout_(int32_t *n, double *sx, 
	int32_t *idigit, char *ifmt)
{

    /* System generated locals */
    int32_t i__1, i__2, i__3;

    /* Local variables */
    int32_t i, k1, k2, lll;
    char line[80];
    int32_t ndigit;

/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --sx;

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
	    for (k1 = 1; k1 <= i__1; k1 += 5) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 4;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d:", k1, k2);
		for (i = k1; i <= k2; ++i) {
		    printf(" %12.3f", sx[i]);
		}
	    }
	} else if (ndigit <= 6) {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 4) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 3;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d:", k1, k2);
		for (i = k1; i <= k2; ++i) {
		    printf(" %14.5f", sx[i]);
		}
	    }
	} else if (ndigit <= 10) {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 3) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 2;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d:", k1, k2);
		for (i = k1; i <= k2; ++i) {
		    printf(" %18.9f", sx[i]);
		}
	    }
	} else {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 2) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 1;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d:", k1, k2);
		for (i = k1; i <= k2; ++i) {
		    printf(" %24.13f", sx[i]);
		}
	    }
	}

/* ======================================================================= */
/*             CODE FOR OUTPUT USING 132 COLUMNS FORMAT */
/* ======================================================================= */

    } else {
	if (ndigit <= 4) {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 10) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 9;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d:", k1, k2);
		for (i = k1; i <= k2; ++i) {
		    printf(" %12.3f", sx[i]);
		}
	    }
	} else if (ndigit <= 6) {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 8) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 7;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d:", k1, k2);
		for (i = k1; i <= k2; ++i) {
		    printf(" %14.5f", sx[i]);
		}
	    }
	} else if (ndigit <= 10) {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 6) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 5;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d:", k1, k2);
		for (i = k1; i <= k2; ++i) {
		    printf(" %18.9f", sx[i]);
		}
	    }
	} else {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 5) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 4;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d:", k1, k2);
		for (i = k1; i <= k2; ++i) {
		    printf(" %24.13f", sx[i]);
		}
	    }
	}
    }
    printf("\n");
    return 0;
} /* dvout_ */

