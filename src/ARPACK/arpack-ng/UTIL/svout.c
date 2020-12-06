/* arpack-ng\UTIL\svout.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \Name: SVOUT
 *        Real vector output routine.
 *
 *     N      - Length of array SX.  (Input)
 *     SX     - Real array to be printed.  (Input)
 *     IFMT   - Format to be used in printing array SX.  (Input)
 *     IDIGIT - Print up to IABS(IDIGIT) decimal digits per number.  (In)
 *              If IDIGIT .LT. 0, printing is done with 72 columns.
 *              If IDIGIT .GT. 0, printing is done with 132 columns.
 */
int svout_(int *n, float *sx, int *idigit, char *ifmt)
{

    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    int i, k1, k2, lll;
    char line[80];
    int ndigit;
    int len = *n;

    /* Parameter adjustments */
    --sx;

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
        for (k1 = 1; k1 <= len; k1 += 5)
        {
            /* Computing MIN */
            i__2 = *n, i__3 = k1 + 4;
            k2 = min(i__2,i__3);
            printf("\n  %4d - %4d:", k1, k2);
            for (i = k1; i <= k2; ++i)
            {
                printf(" %12.3e", sx[i]);
            }
        }
    }
    else if (ndigit <= 6)
    {
        for (k1 = 1; k1 <= len; k1 += 4)
        {
            /* Computing MIN */
            i__3 = k1 + 3;
            k2 = min(len,i__3);
            printf("\n  %4d - %4d:", k1, k2);
            for (i = k1; i <= k2; ++i)
            {
                printf(" %14.5e", sx[i]);
            }
        }
    }
    else if (ndigit <= 10)
    {
        for (k1 = 1; k1 <= len; k1 += 3)
        {
            /* Computing MIN */
            i__3 = k1 + 2;
            k2 = min(len,i__3);
            printf("\n  %4d - %4d:", k1, k2);
            for (i = k1; i <= k2; ++i)
            {
                printf(" %18.9e", sx[i]);
            }
        }
    }
    else
    {
        for (k1 = 1; k1 <= len; k1 += 2)
        {
            /* Computing MIN */
            i__3 = k1 + 1;
            k2 = min(len,i__3);
            printf("\n  %4d - %4d:", k1, k2);
            for (i = k1; i <= k2; ++i)
            {
                printf(" %24.13e", sx[i]);
            }
        }
    }

    printf("\n");
    return 0;
} /* svout_ */

