/* arpack-ng\UTIL\zvout.f -- translated by f2c (version 20100827). */

#include "arpack_internal.h"

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
int zvout_(const int n, const a_dcomplex *cx, const int idigit, const char *ifmt)
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

    /* =============================================================== */
    /*         CODE FOR OUTPUT USING 72 COLUMNS FORMAT                 */
    /* =============================================================== */

    if (ndigit <= 4)
    {
        for (j = 0; j < n; j += 2)
        {
            m = j + 2;
            k = min(n,m);
            printf("\n  %4d - %4d:  ", j, k);
            for (i = j; i < k; ++i)
            {
                printf("(  %10.3e   %10.3e)  ", cx[i].r, cx[i].i);
            }
        }
    }
    else if (ndigit <= 6)
    {
        for (j = 0; j < n; j += 2)
        {
            m = j + 2;
            k = min(n,m);
            printf("\n  %4d - %4d:  ", j, k);
            for (i = j; i < k; ++i)
            {
                printf("(  %12.5e   %12.5e)  ", cx[i].r, cx[i].i);
            }
        }
    }
    else if (ndigit <= 8)
    {
        for (j = 0; j < n; j += 2)
        {
            m = j + 2;
            k = min(n,m);
            printf("\n  %4d - %4d:  ", j, k);
            for (i = j; i < k; ++i)
            {
                printf("(  %14.7e   %14.7e)  ", cx[i].r, cx[i].i);
            }
        }
    }
    else
    {
        for (j = 0; j < n; ++j)
        {
            printf("\n  %4d:  (  %20.13e   %20.13e)  ", j, cx[i].r, cx[i].i);
        }
    }

    printf("\n");
    return 0;

} /* zvout_ */

