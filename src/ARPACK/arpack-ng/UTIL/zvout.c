/* arpack-ng\UTIL\zvout.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \Name: ZVOUT
 *        Complex*16 vector output routine.
 *
 *     N      - Length of array CX.  (Input) */
/*     CX     - Complex*16 array to be printed.  (Input)
 *     IFMT   - Format to be used in printing array CX.  (Input)
 *     IDIGIT - Print up to IABS(IDIGIT) decimal digits per number.  (In)
 *              If IDIGIT .LT. 0, printing is done with 72 columns.
 *              If IDIGIT .GT. 0, printing is done with 132 columns.
 */
int zvout_(int *n, zomplex *cx, int *idigit, char *ifmt)
{

    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    int i, k1, k2, lll;
    char line[80];
    int ndigit;
    int len = *n;

    /* Parameter adjustments */
    --cx;

    /* Function Body */
    /* Computing MIN */
    i__1 = strlen(ifmt);
    lll = min(i__1,79);
    i__1 = lll;
    for (i = 1; i <= i__1; ++i)
    {
        line[i - 1] = '-';
    }
    line[lll] = '\0';

    printf("\n %s\n %s", ifmt, line);

    if (*n <= 0)
    {
        return 0;
    }
    ndigit = *idigit;
    if (ndigit == 0)
    {
        ndigit = 4;
    }
    else if (ndigit < 0)
    {
        /* 132 COLUMNS FORMAT WAS REMOVED */
        ndigit = -ndigit;
    }

    /* =============================================================== */
    /*         CODE FOR OUTPUT USING 72 COLUMNS FORMAT                 */
    /* =============================================================== */

    if (ndigit <= 4)
    {
        for (k1 = 1; k1 <= len; k1 += 2)
        {
            /* Computing MIN */
            i__3 = k1 + 1;
            k2 = min(len,i__3);
            printf("\n  %4d - %4d:  ", k1, k2);
            for (i = k1; i <= k2; ++i)
            {
                printf("(  %10.3e   %10.3e)  ", cx[i].r, cx[i].i);
            }
        }
    }
    else if (ndigit <= 6)
    {
        for (k1 = 1; k1 <= len; k1 += 2)
        {
            /* Computing MIN */
            i__3 = k1 + 1;
            k2 = min(len,i__3);
            printf("\n  %4d - %4d:  ", k1, k2);
            for (i = k1; i <= k2; ++i)
            {
                printf("(  %12.5e   %12.5e)  ", cx[i].r, cx[i].i);
            }
        }
    }
    else if (ndigit <= 8)
    {
        for (k1 = 1; k1 <= len; k1 += 2)
        {
            /* Computing MIN */
            i__3 = k1 + 1;
            k2 = min(len,i__3);
            printf("\n  %4d - %4d:  ", k1, k2);
            for (i = k1; i <= k2; ++i)
            {
                printf("(  %14.7e   %14.7e)  ", cx[i].r, cx[i].i);
            }
        }
    }
    else
    {
        for (k1 = 1; k1 <= len; ++k1)
        {
            printf("  %4d:  (  %20.13e   %20.13e)  ", k1, cx[i].r, cx[i].i);
        }
    }

    printf("\n");
    return 0;

} /* zvout_ */

