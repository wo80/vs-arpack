/* arpack-ng\SRC\ssesrt.f -- translated by f2c (version 20100827). */

#include "arpack_internal.h"

/**
 * \BeginDoc
 *
 * \Name: ssesrt
 *
 * \Description:
 *  Sort the array X in the order specified by WHICH and optionally
 *  apply the permutation to the columns of the matrix A.
 *
 * \Usage:
 *  call ssesrt
 *     ( WHICH, APPLY, N, X, NA, A, LDA)
 *
 * \Arguments
 *  WHICH   Character*2.  (Input)
 *          'LM' -> X is sorted into increasing order of magnitude.
 *          'SM' -> X is sorted into decreasing order of magnitude.
 *          'LA' -> X is sorted into increasing order of algebraic.
 *          'SA' -> X is sorted into decreasing order of algebraic.
 *
 *  APPLY   Logical.  (Input)
 *          APPLY = .TRUE.  -> apply the sorted order to A.
 *          APPLY = .FALSE. -> do not apply the sorted order to A.
 *
 *  N       Integer.  (INPUT)
 *          Dimension of the array X.
 *
 *  X      Real array of length N.  (INPUT/OUTPUT)
 *          The array to be sorted.
 *
 *  NA      Integer.  (INPUT)
 *          Number of rows of the matrix A.
 *
 *  A      Real array of length NA by N.  (INPUT/OUTPUT)
 *
 *  LDA     Integer.  (INPUT)
 *          Leading dimension of A.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Routines
 *     sswap  Level 1 BLAS that swaps the contents of two vectors.
 *
 * \Authors
 *     Danny Sorensen               Phuong Vu
 *     Richard Lehoucq              CRPC / Rice University
 *     Dept. of Computational &     Houston, Texas
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * \Revision history:
 *     12/15/93: Version ' 2.1'.
 *               Adapted from the sort routine in LANSO and
 *               the ARPACK code ssortr
 *
 * \SCCS Information: @(#)
 * FILE: sesrt.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2
 *
 * \EndLib
 */
int ssesrt_(char *which, bool *apply, int *n, float *x,int *na, float *a,
            int *lda)
{
    /* System generated locals */
    int a_dim, a_offset, i__1;
    float r__1, r__2;

    /* Local variables */
    int i, j, igap;
    float temp;

    /* Parameter adjustments */
    a_dim = *lda;
    a_offset = 1 + a_dim * 0;
    a -= a_offset;

    /* Function Body */
    igap = *n / 2;

    if (strcmp(which, "SA") == 0)
    {
        /*        X is sorted into decreasing order of algebraic. */

L10:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i = igap; i <= i__1; ++i)
        {
            j = i - igap;
L20:

            if (j < 0)
            {
                goto L30;
            }

            if (x[j] < x[j + igap])
            {
                temp = x[j];
                x[j] = x[j + igap];
                x[j + igap] = temp;
                if (*apply)
                {
                    sswap_(na, &a[j * a_dim + 1], &c__1, &a[(j + igap) * a_dim + 1], &c__1);
                }
            }
            else
            {
                goto L30;
            }
            j -= igap;
            goto L20;
L30:
            ;
        }
        igap /= 2;
        goto L10;
    }
    else if (strcmp(which, "SM") == 0)
    {
        /*        X is sorted into decreasing order of magnitude. */

L40:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i = igap; i <= i__1; ++i)
        {
            j = i - igap;
L50:

            if (j < 0)
            {
                goto L60;
            }

            if ((r__1 = x[j], dabs(r__1)) < (r__2 = x[j + igap], dabs(r__2)))
            {
                temp = x[j];
                x[j] = x[j + igap];
                x[j + igap] = temp;
                if (*apply)
                {
                    sswap_(na, &a[j * a_dim + 1], &c__1, &a[(j + igap) * a_dim + 1], &c__1);
                }
            }
            else
            {
                goto L60;
            }
            j -= igap;
            goto L50;
L60:
            ;
        }
        igap /= 2;
        goto L40;
    }
    else if (strcmp(which, "LA") == 0)
    {
        /*        X is sorted into increasing order of algebraic. */

L70:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i = igap; i <= i__1; ++i)
        {
            j = i - igap;
L80:

            if (j < 0)
            {
                goto L90;
            }

            if (x[j] > x[j + igap])
            {
                temp = x[j];
                x[j] = x[j + igap];
                x[j + igap] = temp;
                if (*apply)
                {
                    sswap_(na, &a[j * a_dim + 1], &c__1, &a[(j + igap) * a_dim + 1], &c__1);
                }
            }
            else
            {
                goto L90;
            }
            j -= igap;
            goto L80;
L90:
            ;
        }
        igap /= 2;
        goto L70;
    }
    else if (strcmp(which, "LM") == 0)
    {
        /*        X is sorted into increasing order of magnitude. */

L100:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i = igap; i <= i__1; ++i)
        {
            j = i - igap;
L110:

            if (j < 0)
            {
                goto L120;
            }

            if ((r__1 = x[j], dabs(r__1)) > (r__2 = x[j + igap], dabs(r__2)))
            {
                temp = x[j];
                x[j] = x[j + igap];
                x[j + igap] = temp;
                if (*apply)
                {
                    sswap_(na, &a[j * a_dim + 1], &c__1, &a[(j + igap) * a_dim + 1], &c__1);
                }
            }
            else
            {
                goto L120;
            }
            j -= igap;
            goto L110;
L120:
            ;
        }
        igap /= 2;
        goto L100;
    }

L9000:
    return 0;

    /* ------------- */
    /* End of ssesrt */
    /* ------------- */

} /* ssesrt_ */

