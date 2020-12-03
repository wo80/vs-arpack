/* arpack-ng\UTIL\ivout.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/* ----------------------------------------------------------------------- */
/*  Routine:    IVOUT */

/*  Purpose:    Integer vector output routine. */

/*  Usage:      CALL IVOUT (LOUT, N, IX, IDIGIT, IFMT) */

/*  Arguments */
/*     N      - Length of array IX. (Input) */
/*     IX     - Integer array to be printed. (Input) */
/*     IFMT   - Format to be used in printing array IX. (Input) */
/*     IDIGIT - Print up to ABS(IDIGIT) decimal digits / number. (Input) */
/*              If IDIGIT .LT. 0, printing is done with 72 columns. */
/*              If IDIGIT .GT. 0, printing is done with 132 columns. */

/* ----------------------------------------------------------------------- */

int ivout_(int32_t *n, int32_t *ix, int32_t *
	idigit, char *ifmt)
{

    /* System generated locals */
    int32_t i__1, i__2, i__3;

    /* Local variables */
    int32_t i, k1, k2, lll;
    char line[80];
    int32_t ndigit;

    /* Parameter adjustments */
    --ix;

    /* Function Body */
/* Computing MIN */
    i__1 = strlen(ifmt);
    lll = min(i__1,80);
    i__1 = lll;
    for (i = 1; i <= i__1; ++i) {
	*&line[i - 1] = '-';
/* L1: */
    }

    for (i = lll + 1; i <= 80; ++i) {
	*&line[i - 1] = ' ';
/* L2: */
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
	    for (k1 = 1; k1 <= i__1; k1 += 10) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 9;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d: ", k1, k2);
		for (i = k1; i <= k2; ++i) {
		    printf("   %5d", ix[i]);
		}
	    }

	} else if (ndigit <= 6) {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 7) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 6;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d: ", k1, k2);
		for (i = k1; i <= k2; ++i) {
		    printf("  %7d", ix[i]);
		}
	    }

	} else if (ndigit <= 10) {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 5) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 4;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d: ", k1, k2);
		for (i = k1; i <= k2; ++i) {
		    printf("  %11d", ix[i]);
		}
	    }

	} else {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 3) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 2;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d: ", k1, k2);
		for (i = k1; i <= k2; ++i) {
		    printf("  %15d", ix[i]);
		}
	    }
	}

/* ======================================================================= */
/*             CODE FOR OUTPUT USING 132 COLUMNS FORMAT */
/* ======================================================================= */

    } else {

	if (ndigit <= 4) {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 20) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 19;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d: ", k1, k2);
		for (i = k1; i <= k2; ++i) {
		    printf("   %5d", ix[i]);
		}
	    }

	} else if (ndigit <= 6) {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 15) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 14;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d: ", k1, k2);
		for (i = k1; i <= k2; ++i) {
		    printf("  %7d", ix[i]);
		}
	    }

	} else if (ndigit <= 10) {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 10) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 9;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d: ", k1, k2);
		for (i = k1; i <= k2; ++i) {
		    printf("  %11d", ix[i]);
		}
	    }

	} else {
	    i__1 = *n;
	    for (k1 = 1; k1 <= i__1; k1 += 7) {
/* Computing MIN */
		i__2 = *n, i__3 = k1 + 6;
		k2 = min(i__2,i__3);
		printf("\n  %4d -   %4d: ", k1, k2);
		for (i = k1; i <= k2; ++i) {
		    printf("  %15d", ix[i]);
		}
	    }
	}
    }
    printf("\n");

    return 0;
} /* ivout_ */

