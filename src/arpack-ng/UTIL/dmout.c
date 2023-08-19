/* arpack-ng\UTIL\dmout.f -- translated by f2c (version 20100827). */

#include "arpack_internal.h"

/**
 * \Name: DMOUT
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
int dmout_(int *m, int *n, double *a, int *lda, int *idigit, char *ifmt)
{
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
        for (k1 = 0; k1 < cols; k1 += 5)
        {
            p = k1 + 5;
            k2 = min(cols,p);
            printf("\n          ");
            for (i = k1; i < k2; ++i)
            {
                printf("     Col %4d ", i);
            }
            for (i = 0; i < rows; ++i)
            {
                printf("\n  Row %4d:  ", i);
                for (j = k1; j < k2; ++j)
                {
                    printf(" %12.3e", a[i + j * a_dim]);
                }
            }
        }
    }
    else if (ndigit <= 6)
    {
        for (k1 = 0; k1 < cols; k1 += 4)
        {
            p = k1 + 4;
            k2 = min(cols,p);
            printf("\n          ");
            for (i = k1; i < k2; ++i)
            {
                printf("      Col %4d  ", i);
            }
            for (i = 0; i < rows; ++i)
            {
                printf("\n  Row %4d:  ", i);
                for (j = k1; j < k2; ++j)
                {
                    printf(" %14.5e", a[i + j * a_dim]);
                }
            }
        }
    }
    else if (ndigit <= 10)
    {
        for (k1 = 0; k1 < cols; k1 += 3)
        {
            p = k1 + 3;
            k2 = min(cols,p);
            printf("\n          ");
            for (i = k1; i < k2; ++i)
            {
                printf("        Col %4d    ", i);
            }
            for (i = 0; i < rows; ++i)
            {
                printf("\n  Row %4d:  ", i);
                for (j = k1; j < k2; ++j)
                {
                    printf(" %18.9e", a[i + j * a_dim]);
                }
            }
        }
    }
    else
    {
        for (k1 = 0; k1 < cols; k1 += 2)
        {
            p = k1 + 2;
            k2 = min(cols,p);
            printf("\n          ");
            for (i = k1; i < k2; ++i)
            {
                printf("          Col %4d      ", i);
            }
            for (i = 0; i < rows; ++i)
            {
                printf("\n  Row %4d:  ", i);
                for (j = k1; j < k2; ++j)
                {
                    printf(" %22.13e", a[i + j * a_dim]);
                }
            }
        }
    }

    printf("\n");

    return 0;
} /* dmout_ */

