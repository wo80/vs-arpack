/* arpack-ng\SRC\zsortc.f -- translated by f2c (version 20100827). */

#include "arpack_internal.h"
/**
 * \BeginDoc
 *
 * \Name: zsortc
 *
 * \Description:
 *  Sorts the Complex*16 array in X into the order
 *  specified by WHICH and optionally applies the permutation to the
 *  Double precision  array Y.
 *
 * \Usage:
 *  call zsortc
 *     ( WHICH, APPLY, N, X, Y )
 *
 * \Arguments
 *  WHICH   Character*2.  (Input)
 *          'LM' -> sort X into increasing order of magnitude.
 *          'SM' -> sort X into decreasing order of magnitude.
 *          'LR' -> sort X with real(X) in increasing algebraic order
 *          'SR' -> sort X with real(X) in decreasing algebraic order
 *          'LI' -> sort X with imag(X) in increasing algebraic order
 *          'SI' -> sort X with imag(X) in decreasing algebraic order
 *
 *  APPLY   Logical.  (Input)
 *          APPLY = .TRUE.  -> apply the sorted order to array Y.
 *          APPLY = .FALSE. -> do not apply the sorted order to array Y.
 *
 *  N       Integer.  (INPUT)
 *          Size of the arrays.
 *
 *  X       Complex*16 array of length N.  (INPUT/OUTPUT)
 *          This is the array to be sorted.
 *
 *  Y       Complex*16 array of length N.  (INPUT/OUTPUT)
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Routines called:
 *     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *
 * \Author
 *     Danny Sorensen               Phuong Vu
 *     Richard Lehoucq              CRPC / Rice University
 *     Dept. of Computational &     Houston, Texas
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 *     Adapted from the sort routine in LANSO.
 *
 * \SCCS Information: @(#)
 * FILE: sortc.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2
 *
 * \EndLib
 */
int zsortc_(char *which, logical *apply, int *n, a_dcomplex *x, a_dcomplex *y)
{
    /* System generated locals */
    int i__1, k;
    double d__1, d__2;

    /* Local variables */
    int i, j, igap;
    a_dcomplex temp;
    double temp1, temp2;


    igap = *n / 2;

    if (strcmp(which, "LM") == 0)
    {
        /* ------------------------------------------ */
        /* Sort X into increasing order of magnitude. */
        /* ------------------------------------------ */

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

            d__1 = x[j].r;
            d__2 = x[j].i;
            temp1 = dlapy2_(&d__1, &d__2);
            k = j + igap;
            d__1 = x[k].r;
            d__2 = x[k].i;
            temp2 = dlapy2_(&d__1, &d__2);

            if (temp1 > temp2)
            {
                temp.r = x[j].r, temp.i = x[j].i;
                x[j].r = x[k].r, x[j].i = x[k].i;
                x[k].r = temp.r, x[k].i = temp.i;

                if (*apply)
                {
                    temp.r = y[j].r, temp.i = y[j].i;
                    y[j].r = y[k].r, y[j].i = y[k].i;
                    y[k].r = temp.r, y[k].i = temp.i;
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
        /* ------------------------------------------ */
        /* Sort X into decreasing order of magnitude. */
        /* ------------------------------------------ */

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

            d__1 = x[j].r;
            d__2 = x[j].i;
            temp1 = dlapy2_(&d__1, &d__2);
            k = j + igap;
            d__1 = x[k].r;
            d__2 = x[k].i;
            temp2 = dlapy2_(&d__1, &d__2);

            if (temp1 < temp2)
            {
                temp.r = x[j].r, temp.i = x[j].i;
                x[j].r = x[k].r, x[j].i = x[k].i;
                x[k].r = temp.r, x[k].i = temp.i;

                if (*apply)
                {
                    temp.r = y[j].r, temp.i = y[j].i;
                    y[j].r = y[k].r, y[j].i = y[k].i;
                    y[k].r = temp.r, y[k].i = temp.i;
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
    else if (strcmp(which, "LR") == 0)
    {
        /* ---------------------------------------------- */
        /* Sort XREAL into increasing order of algebraic. */
        /* ---------------------------------------------- */

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

            k = j + igap;
            if (x[j].r > x[k].r)
            {
                temp.r = x[j].r, temp.i = x[j].i;
                x[j].r = x[k].r, x[j].i = x[k].i;
                x[k].r = temp.r, x[k].i = temp.i;

                if (*apply)
                {
                    temp.r = y[j].r, temp.i = y[j].i;
                    y[j].r = y[k].r, y[j].i = y[k].i;
                    y[k].r = temp.r, y[k].i = temp.i;
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
    else if (strcmp(which, "SR") == 0)
    {
        /* ---------------------------------------------- */
        /* Sort XREAL into decreasing order of algebraic. */
        /* ---------------------------------------------- */

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

            k = j + igap;
            if (x[j].r < x[k].r)
            {
                temp.r = x[j].r, temp.i = x[j].i;
                x[j].r = x[k].r, x[j].i = x[k].i;
                x[k].r = temp.r, x[k].i = temp.i;

                if (*apply)
                {
                    temp.r = y[j].r, temp.i = y[j].i;
                    y[j].r = y[k].r, y[j].i = y[k].i;
                    y[k].r = temp.r, y[k].i = temp.i;
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
    else if (strcmp(which, "LI") == 0)
    {
        /* ------------------------------------------ */
        /* Sort XIMAG into increasing algebraic order */
        /* ------------------------------------------ */

L130:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i = igap; i <= i__1; ++i)
        {
            j = i - igap;
L140:

            if (j < 0)
            {
                goto L150;
            }

            k = j + igap;
            if (x[j].i > x[k].i)
            {
                temp.r = x[j].r, temp.i = x[j].i;
                x[j].r = x[k].r, x[j].i = x[k].i;
                x[k].r = temp.r, x[k].i = temp.i;

                if (*apply)
                {
                    temp.r = y[j].r, temp.i = y[j].i;
                    y[j].r = y[k].r, y[j].i = y[k].i;
                    y[k].r = temp.r, y[k].i = temp.i;
                }
            }
            else
            {
                goto L150;
            }
            j -= igap;
            goto L140;
L150:
            ;
        }
        igap /= 2;
        goto L130;
    }
    else if (strcmp(which, "SI") == 0)
    {
        /* ------------------------------------------- */
        /* Sort XIMAG into decreasing algebraic order  */
        /* ------------------------------------------- */

L160:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i = igap; i <= i__1; ++i)
        {
            j = i - igap;
L170:

            if (j < 0)
            {
                goto L180;
            }

            k = j + igap;
            if (x[j].i < x[k].i)
            {
                temp.r = x[j].r, temp.i = x[j].i;
                x[j].r = x[k].r, x[j].i = x[k].i;
                x[k].r = temp.r, x[k].i = temp.i;

                if (*apply)
                {
                    temp.r = y[j].r, temp.i = y[j].i;
                    y[j].r = y[k].r, y[j].i = y[k].i;
                    y[k].r = temp.r, y[k].i = temp.i;
                }
            }
            else
            {
                goto L180;
            }
            j -= igap;
            goto L170;
L180:
            ;
        }
        igap /= 2;
        goto L160;
    }

L9000:
    return 0;

    /* ------------- */
    /* End of zsortc */
    /* ------------- */

} /* zsortc_ */

