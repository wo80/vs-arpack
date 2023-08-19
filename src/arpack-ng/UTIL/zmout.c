/* arpack-ng\UTIL\zmout.f -- translated by f2c (version 20100827). */

#include "arpack_internal.h"

/**
 * \Name: ZMOUT
 *        Complex*16 matrix output routine.
 *
 *     M      - Number of rows of A.  (Input)
 *     N      - Number of columns of A.  (Input)
 *     A      - Complex*16 M by N matrix to be printed.  (Input)
 *     LDA    - Leading dimension of A exactly as specified in the
 *              dimension statement of the calling program.  (Input)
 *     IFMT   - Format to be used in printing matrix A.  (Input)
 *     IDIGIT - Print up to IABS(IDIGIT) decimal digits per number.  (In)
 *              If IDIGIT .LT. 0, printing is done with 72 columns.
 *              If IDIGIT .GT. 0, printing is done with 132 columns.
 */
int zmout_(int *m, int *n, zomplex *a, int *lda, int *idigit, char *ifmt)
{
    /* System generated locals */
    zomplex d__1;

    /* Local variables */
    int i, j, k1, k2, len, p;
    char line[80];
    int ndigit;
    int cols = *n;
    int rows = *m;

    int a_dim = *lda;

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

    if (*m <= 0 || *n <= 0 || *lda <= 0)
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
        for (k1 = 0; k1 < cols; k1 += 2)
        {
            p = k1 + 2;
            k2 = min(cols, p);
            printf("\n           ");
            for (i = k1; i < k2; ++i)
            {
                printf("          Col %4d         ", i);
            }
            for (i = 0; i < rows; ++i)
            {
                printf("\n  Row %4d:  ", i);
                for (j = k1; j < k2; ++j)
                {
                    d__1 = a[i + j * a_dim];
                    printf("(  %10.3e   %10.3e)  ", d__1.r, d__1.i);
                }
            }
        }
    }
    else if (ndigit <= 6)
    {
        for (k1 = 0; k1 < cols; k1 += 2)
        {
            p = k1 + 2;
            k2 = min(cols,p);
            printf("\n          ");
            for (i = k1; i < k2; ++i)
            {
                printf("            Col %4d           ", i);
            }
            for (i = 0; i < rows; ++i)
            {
                printf("\n  Row %4d:  ", i);
                for (j = k1; j < k2; ++j)
                {
                    d__1 = a[i + j * a_dim];
                    printf("(  %12.5e   %12.5e)  ", d__1.r, d__1.i);
                }
            }
        }
    }
    else if (ndigit <= 8)
    {
        for (k1 = 0; k1 < cols; k1 += 2)
        {
            p = k1 + 2;
            k2 = min(cols,p);
            printf("\n          ");
            for (i = k1; i < k2; ++i)
            {
                printf("              Col %4d             ", i);
            }
            for (i = 0; i < rows; ++i)
            {
                printf("\n  Row %4d:  ", i);
                for (j = k1; j < k2; ++j)
                {
                    d__1 = a[i + j * a_dim];
                    printf("(  %14.7e   %14.7e)  ", d__1.r, d__1.i);
                }
            }
        }
    }
    else
    {
        for (k1 = 0; k1 < cols; ++k1)
        {
            printf("                               Col %4d                  ", i);
            for (i = 0; i < rows; ++i)
            {
                d__1 = a[i + k1 * a_dim];
                printf("\n  Row %4d:  (  %20.13e   %20.13e)", i, d__1.r, d__1.i);
            }
        }
    }

    printf("\n");

    return 0;
} /* zmout_ */

