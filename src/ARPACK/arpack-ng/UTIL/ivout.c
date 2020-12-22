/* arpack-ng\UTIL\ivout.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \Name: IVOUT
 *        Integer vector output routine.
 *
 *     N      - Length of array IX. (Input)
 *     IX     - Integer array to be printed. (Input)
 *     IFMT   - Format to be used in printing array IX. (Input)
 *     IDIGIT - Print up to ABS(IDIGIT) decimal digits / number. (Input)
 *              If IDIGIT .LT. 0, printing is done with 72 columns.
 *              If IDIGIT .GT. 0, printing is done with 132 columns.
 */
int ivout_(int *n, int *ix, int *idigit, char *ifmt)
{
    /* Local variables */
    int i, k1, k2, len, m;
    char line[80];
    int ndigit;

    /* Function Body */
    len = strlen(ifmt);
    len = min(len,80);
    for (i = 0; i < len; ++i)
    {
        line[i] = '-';
    }
    i = min(i,79);
    line[i] = '\0';

    printf("\n %s\n %s", ifmt, line);

    len = *n;
    if (len <= 0)
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

    /* ======================================================================= */
    /*             CODE FOR OUTPUT USING 72 COLUMNS FORMAT */
    /* ======================================================================= */


    if (ndigit <= 4)
    {
        for (k1 = 0; k1 < len; k1 += 10)
        {
            m = k1 + 10;
            k2 = min(len,m);
            printf("\n  %4d - %4d: ", k1, k2);
            for (i = k1; i < k2; ++i)
            {
                printf("   %5d", ix[i]);
            }
        }
    }
    else if (ndigit <= 6)
    {
        for (k1 = 0; k1 < len; k1 += 7)
        {
            m = k1 + 7;
            k2 = min(len,m);
            printf("\n  %4d - %4d: ", k1, k2);
            for (i = k1; i < k2; ++i)
            {
                printf("  %7d", ix[i]);
            }
        }
    }
    else if (ndigit <= 10)
    {
        for (k1 = 0; k1 < len; k1 += 5)
        {
            m = k1 + 5;
            k2 = min(len,m);
            printf("\n  %4d - %4d: ", k1, k2);
            for (i = k1; i < k2; ++i)
            {
                printf("  %11d", ix[i]);
            }
        }
    }
    else
    {
        for (k1 = 0; k1 < len; k1 += 3)
        {
            m = k1 + 3;
            k2 = min(len,m);
            printf("\n  %4d - %4d: ", k1, k2);
            for (i = k1; i < k2; ++i)
            {
                printf("  %15d", ix[i]);
            }
        }
    }

    printf("\n");

    return 0;
} /* ivout_ */

