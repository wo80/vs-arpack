/* arpack-ng\UTIL\ivout.f -- translated by f2c (version 20100827). */

#include "arpack_internal.h"

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
int ivout_(const int n, const int *ix, const int idigit, const char *ifmt)
{
    /* Local variables */
    int i, j, k, len, m;
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

    if (n <= 0)
    {
        return 0;
    }
    ndigit = idigit;
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
        for (j = 0; j < n; j += 10)
        {
            m = j + 10;
            k = min(n,m);
            printf("\n  %4d - %4d: ", j, k);
            for (i = j; i < k; ++i)
            {
                printf("   %5d", ix[i]);
            }
        }
    }
    else if (ndigit <= 6)
    {
        for (j = 0; j < n; j += 7)
        {
            m = j + 7;
            k = min(n,m);
            printf("\n  %4d - %4d: ", j, k);
            for (i = j; i < k; ++i)
            {
                printf("  %7d", ix[i]);
            }
        }
    }
    else if (ndigit <= 10)
    {
        for (j = 0; j < n; j += 5)
        {
            m = j + 5;
            k = min(n,m);
            printf("\n  %4d - %4d: ", j, k);
            for (i = j; i < k; ++i)
            {
                printf("  %11d", ix[i]);
            }
        }
    }
    else
    {
        for (j = 0; j < n; j += 3)
        {
            m = j + 3;
            k = min(n,m);
            printf("\n  %4d - %4d: ", j, k);
            for (i = j; i < k; ++i)
            {
                printf("  %15d", ix[i]);
            }
        }
    }

    printf("\n");

    return 0;
} /* ivout_ */

