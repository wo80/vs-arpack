/* arpack-ng\UTIL\smout.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \Name: SMOUT
 *        Real matrix output routine.
 *
 *     M      - Number of rows of A.  (Input)
 *     N      - Number of columns of A.  (Input)
 *     A      - Real M by N matrix to be printed.  (Input)
 *     LDA    - Leading dimension of A exactly as specified in the
 *              dimension statement of the calling program.  (Input)
 *     IFMT   - Format to be used in printing matrix A.  (Input)
 *     IDIGIT - Print up to IABS(IDIGIT) decimal digits per number.  (In)
 *              If IDIGIT .LT. 0, printing is done with 72 columns.
 *              If IDIGIT .GT. 0, printing is done with 132 columns.
 */
int smout_(int *m, int *n, float *a, int *lda, int *idigit, char *ifmt)
{
    /* Initialized data */

    /* System generated locals */
    int a_dim, a_offset, i__1, i__3;

    /* Local variables */
    int i, j, k1, k2, lll;
    char line[80];
    int ndigit;
    int cols = *n;
    int rows = *m;

    /* Parameter adjustments */
    a_dim = *lda;
    a_offset = 1 + a_dim;
    a -= a_offset;

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
        for (k1 = 1; k1 <= cols; k1 += 5)
        {
            /* Computing MIN */
            i__3 = k1 + 4;
            k2 = min(cols,i__3);
            printf("\n          ");
            for (i = k1; i <= k2; ++i)
            {
                printf("     Col %4d ", i);
            }
            for (i = 1; i <= rows; ++i)
            {
                printf("\n  Row %4d:  ", i);
                for (j = k1; j <= k2; ++j)
                {
                    printf(" %12.3e", a[i + j * a_dim]);
                }
            }
        }
    }
    else if (ndigit <= 6)
    {
        for (k1 = 1; k1 <= cols; k1 += 4)
        {
            /* Computing MIN */
            i__3 = k1 + 3;
            k2 = min(cols,i__3);
            printf("\n          ");
            for (i = k1; i <= k2; ++i)
            {
                printf("      Col %4d  ", i);
            }
            for (i = 1; i <= rows; ++i)
            {
                printf("\n  Row %4d:  ", i);
                for (j = k1; j <= k2; ++j)
                {
                    printf(" %14.5e", a[i + j * a_dim]);
                }
            }
        }
    }
    else if (ndigit <= 10)
    {
        for (k1 = 1; k1 <= cols; k1 += 3)
        {
            /* Computing MIN */
            i__3 = k1 + 2;
            k2 = min(cols,i__3);
            printf("\n          ");
            for (i = k1; i <= k2; ++i)
            {
                printf("        Col %4d    ", i);
            }
            for (i = 1; i <= rows; ++i)
            {
                printf("\n  Row %4d:  ", i);
                for (j = k1; j <= k2; ++j)
                {
                    printf(" %18.9e", a[i + j * a_dim]);
                }
            }
        }
    }
    else
    {
        for (k1 = 1; k1 <= cols; k1 += 2)
        {
            /* Computing MIN */
            i__3 = k1 + 1;
            k2 = min(cols,i__3);
            printf("\n          ");
            for (i = k1; i <= k2; ++i)
            {
                printf("          Col %4d      ", i);
            }
            for (i = 1; i <= rows; ++i)
            {
                printf("\n  Row %4d:  ", i);
                for (j = k1; j <= k2; ++j)
                {
                    printf(" %22.13e", a[i + j * a_dim]);
                }
            }
        }
    }

    printf("\n");

    return 0;
} /* smout_ */

