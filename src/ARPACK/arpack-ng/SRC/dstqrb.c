/* arpack-ng\SRC\dstqrb.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \BeginDoc
 *
 * \Name: dstqrb
 *
 * \Description:
 *  Computes all eigenvalues and the last component of the eigenvectors
 *  of a symmetric tridiagonal matrix using the implicit QL or QR method.
 *
 *  This is mostly a modification of the LAPACK routine dsteqr.
 *  See Remarks.
 *
 * \Usage:
 *  call dstqrb
 *     ( N, D, E, Z, WORK, INFO )
 *
 * \Arguments
 *  N       Integer.  (INPUT)
 *          The number of rows and columns in the matrix.  N >= 0.
 *
 *  D       Double precision array, dimension (N).  (INPUT/OUTPUT)
 *          On entry, D contains the diagonal elements of the
 *          tridiagonal matrix.
 *          On exit, D contains the eigenvalues, in ascending order.
 *          If an error exit is made, the eigenvalues are correct
 *          for indices 1,2,...,INFO-1, but they are unordered and
 *          may not be the smallest eigenvalues of the matrix.
 *
 *  E       Double precision array, dimension (N-1).  (INPUT/OUTPUT)
 *          On entry, E contains the subdiagonal elements of the
 *          tridiagonal matrix in positions 1 through N-1.
 *          On exit, E has been destroyed.
 *
 *  Z       Double precision array, dimension (N).  (OUTPUT)
 *          On exit, Z contains the last row of the orthonormal
 *          eigenvector matrix of the symmetric tridiagonal matrix.
 *          If an error exit is made, Z contains the last row of the
 *          eigenvector matrix associated with the stored eigenvalues.
 *
 *  WORK    Double precision array, dimension (max(1,2*N-2)).  (WORKSPACE)
 *          Workspace used in accumulating the transformation for
 *          computing the last components of the eigenvectors.
 *
 *  INFO    Integer.  (OUTPUT)
 *          = 0:  normal return.
 *          < 0:  if INFO = -i, the i-th argument had an illegal value.
 *          > 0:  if INFO = +i, the i-th eigenvalue has not converged
 *                              after a total of  30*N  iterations.
 *
 * \Remarks
 *  1. None.
 *
 * -----------------------------------------------------------------------
 * \EndDoc
 *
 * \BeginLib
 *
 * \Local variables:
 *     xxxxxx  real
 *
 * \Routines called:
 *     daxpy   Level 1 BLAS that computes a vector triad.
 *     dcopy   Level 1 BLAS that copies one vector to another.
 *     dswap   Level 1 BLAS that swaps the contents of two vectors.
 *     lsame   LAPACK character comparison routine.
 *     dlae2   LAPACK routine that computes the eigenvalues of a 2-by-2
 *             symmetric matrix.
 *     dlaev2  LAPACK routine that eigendecomposition of a 2-by-2 symmetric
 *             matrix.
 *     dlamch  LAPACK routine that determines machine constants.
 *     dlanst  LAPACK routine that computes the norm of a matrix.
 *     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     dlartg  LAPACK Givens rotation construction routine.
 *     dlascl  LAPACK routine for careful scaling of a matrix.
 *     dlaset  LAPACK matrix initialization routine.
 *     dlasr   LAPACK routine that applies an orthogonal transformation to
 *             a matrix.
 *     dlasrt  LAPACK sorting routine.
 *     dsteqr  LAPACK routine that computes eigenvalues and eigenvectors
 *             of a symmetric tridiagonal matrix.
 *     xerbla  LAPACK error handler routine.
 *
 * \Authors
 *     Danny Sorensen               Phuong Vu
 *     Richard Lehoucq              CRPC / Rice University
 *     Dept. of Computational &     Houston, Texas
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * \SCCS Information: @(#)
 * FILE: stqrb.F   SID: 2.5   DATE OF SID: 8/27/96   RELEASE: 2
 *
 * \Remarks
 *     1. Starting with version 2.5, this routine is a modified version
 *        of LAPACK version 2.0 subroutine SSTEQR. No lines are deleted,
 *        only commeted out and new lines inserted.
 *        All lines commented out have "c$$$" at the beginning.
 *        Note that the LAPACK version 1.0 subroutine SSTEQR contained
 *        bugs.
 *
 * \EndLib
 */
int dstqrb_(int *n, double *d, double *e,
            double *z, double *work, int *info)
{
    /* System generated locals */
    int i__1, i__2;
    double d__1, d__2;

    /* Builtin functions */
    double sqrt(double), d_sign(double *, double *);

    /* Local variables */
    double b, c, f, g;
    int i, j, k, l, m;
    double p, r, s;
    int l1, ii, mm, lm1, mm1, nm1;
    double rt1, rt2, eps;
    int lsv;
    double tst, eps2;
    int lend, jtot;
    double anorm;
    int lendm1, lendp1;
    int iscale;
    double safmin;
    double safmax;
    int lendsv, nmaxit, icompz;
    double ssfmax, ssfmin;

    /* Parameter adjustments */
    --work;
    --z;
    --e;
    --d;

    /* Function Body */
    *info = 0;

    /* $$$      IF( LSAME( COMPZ, 'N' ) ) THEN */
    /* $$$         ICOMPZ = 0 */
    /* $$$      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN */
    /* $$$         ICOMPZ = 1 */
    /* $$$      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN */
    /* $$$         ICOMPZ = 2 */
    /* $$$      ELSE */
    /* $$$         ICOMPZ = -1 */
    /* $$$      END IF */
    /* $$$      IF( ICOMPZ.LT.0 ) THEN */
    /* $$$         INFO = -1 */
    /* $$$      ELSE IF( N.LT.0 ) THEN */
    /* $$$         INFO = -2 */
    /* $$$      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1, */
    /* $$$     $         N ) ) ) THEN */
    /* $$$         INFO = -6 */
    /* $$$      END IF */
    /* $$$      IF( INFO.NE.0 ) THEN */
    /* $$$         CALL XERBLA( 'SSTEQR', -INFO ) */
    /* $$$         RETURN */
    /* $$$      END IF */

    /*    *** New starting with version 2.5 *** */

    icompz = 2;
    /*    ************************************* */

    /*     quick return if possible */

    if (*n == 0)
    {
        return 0;
    }

    if (*n == 1)
    {
        if (icompz == 2)
        {
            z[1] = 1.0;
        }
        return 0;
    }

    /*     determine the unit roundoff and over/underflow thresholds. */

    eps = dlamch_("E");
    /* Computing 2nd power */
    d__1 = eps;
    eps2 = d__1 * d__1;
    safmin = dlamch_("S");
    safmax = 1.0 / safmin;
    ssfmax = sqrt(safmax) / 3.0;
    ssfmin = sqrt(safmin) / eps2;

    /*     compute the eigenvalues and eigenvectors of the tridiagonal */
    /*     matrix. */

    /* $$      if( icompz.eq.2 ) */
    /* $$$     $   call dlaset( 'full', n, n, zero, one, z, ldz ) */

    /*     *** New starting with version 2.5 *** */

    if (icompz == 2)
    {
        i__1 = *n - 1;
        for (j = 1; j <= i__1; ++j)
        {
            z[j] = 0.0;

        }
        z[*n] = 1.0;
    }
    /*     ************************************* */

    nmaxit = *n * 30;
    jtot = 0;

    /*     determine where the matrix splits and choose ql or qr iteration */
    /*     for each block, according to whether top or bottom diagonal */
    /*     element is smaller. */

    l1 = 1;
    nm1 = *n - 1;

L10:
    if (l1 > *n)
    {
        goto L160;
    }
    if (l1 > 1)
    {
        e[l1 - 1] = 0.0;
    }
    if (l1 <= nm1)
    {
        i__1 = nm1;
        for (m = l1; m <= i__1; ++m)
        {
            tst = (d__1 = e[m], abs(d__1));
            if (tst == 0.0)
            {
                goto L30;
            }
            if (tst <= sqrt((d__1 = d[m], abs(d__1))) * sqrt((d__2 = d[m
                    + 1], abs(d__2))) * eps)
            {
                e[m] = 0.0;
                goto L30;
            }

        }
    }
    m = *n;

L30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l)
    {
        goto L10;
    }

    /*     scale submatrix in rows and columns l to lend */

    i__1 = lend - l + 1;
    anorm = dlanst_("i", &i__1, &d[l], &e[l]);
    iscale = 0;
    if (anorm == 0.0)
    {
        goto L10;
    }
    if (anorm > ssfmax)
    {
        iscale = 1;
        i__1 = lend - l + 1;
        dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d[l], n, info);
        i__1 = lend - l;
        dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, info);
    }
    else if (anorm < ssfmin)
    {
        iscale = 2;
        i__1 = lend - l + 1;
        dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d[l], n, info);
        i__1 = lend - l;
        dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, info);
    }

    /*     choose between ql and qr iteration */

    if ((d__1 = d[lend], abs(d__1)) < (d__2 = d[l], abs(d__2)))
    {
        lend = lsv;
        l = lendsv;
    }

    if (lend > l)
    {
        /*        ql iteration */

        /*        look for small subdiagonal element. */

L40:
        if (l != lend)
        {
            lendm1 = lend - 1;
            i__1 = lendm1;
            for (m = l; m <= i__1; ++m)
            {
                /* Computing 2nd power */
                d__2 = (d__1 = e[m], abs(d__1));
                tst = d__2 * d__2;
                if (tst <= eps2 * (d__1 = d[m], abs(d__1)) * (d__2 = d[m
                        + 1], abs(d__2)) + safmin)
                {
                    goto L60;
                }
            }
        }

        m = lend;

L60:
        if (m < lend)
        {
            e[m] = 0.0;
        }
        p = d[l];
        if (m == l)
        {
            goto L80;
        }

        /*        if remaining matrix is 2-by-2, use dlae2 or dlaev2 */
        /*        to compute its eigensystem. */

        if (m == l + 1)
        {
            if (icompz > 0)
            {
                dlaev2_(&d[l], &e[l], &d[l + 1], &rt1, &rt2, &c, &s);
                work[l] = c;
                work[*n - 1 + l] = s;
                /* $$$               call dlasr( 'r', 'v', 'b', n, 2, work( l ), */
                /* $$$     $                     work( n-1+l ), z( 1, l ), ldz ) */

                /*              *** New starting with version 2.5 *** */

                tst = z[l + 1];
                z[l + 1] = c * tst - s * z[l];
                z[l] = s * tst + c * z[l];
                /*              ************************************* */
            }
            else
            {
                dlae2_(&d[l], &e[l], &d[l + 1], &rt1, &rt2);
            }
            d[l] = rt1;
            d[l + 1] = rt2;
            e[l] = 0.0;
            l += 2;
            if (l <= lend)
            {
                goto L40;
            }
            goto L140;
        }

        if (jtot == nmaxit)
        {
            goto L140;
        }
        ++jtot;

        /*        form shift. */

        g = (d[l + 1] - p) / (e[l] * 2.);
        r = dlapy2_(&g, &d_one);
        g = d[m] - p + e[l] / (g + d_sign(&r, &g));

        s = 1.0;
        c = 1.0;
        p = 0.0;

        /*        inner loop */

        mm1 = m - 1;
        i__1 = l;
        for (i = mm1; i >= i__1; --i)
        {
            f = s * e[i];
            b = c * e[i];
            dlartg_(&g, &f, &c, &s, &r);
            if (i != m - 1)
            {
                e[i + 1] = r;
            }
            g = d[i + 1] - p;
            r = (d[i] - g) * s + c * 2. * b;
            p = s * r;
            d[i + 1] = g + p;
            g = c * r - b;

            /*           if eigenvectors are desired, then save rotations. */

            if (icompz > 0)
            {
                work[i] = c;
                work[*n - 1 + i] = -s;
            }

        }

        /*        if eigenvectors are desired, then apply saved rotations. */

        if (icompz > 0)
        {
            mm = m - l + 1;
            /* $$$            call dlasr( 'r', 'v', 'b', n, mm, work( l ), work( n-1+l ), */
            /* $$$     $                  z( 1, l ), ldz ) */

            /*             *** New starting with version 2.5 *** */

            dlasr_("r", "v", "b", &c__1, &mm, &work[l], &work[*n - 1 + l], &z[l], &c__1);
            /*             ************************************* */
        }

        d[l] -= p;
        e[l] = g;
        goto L40;

        /*        eigenvalue found. */

L80:
        d[l] = p;

        ++l;
        if (l <= lend)
        {
            goto L40;
        }
        goto L140;
    }
    else
    {
        /*        qr iteration */

        /*        look for small superdiagonal element. */

L90:
        if (l != lend)
        {
            lendp1 = lend + 1;
            i__1 = lendp1;
            for (m = l; m >= i__1; --m)
            {
                /* Computing 2nd power */
                d__2 = (d__1 = e[m - 1], abs(d__1));
                tst = d__2 * d__2;
                d__1 = d[m];
                d__2 = d[m - 1];
                if (tst <= eps2 * abs(d__1) * abs(d__2) + safmin)
                {
                    goto L110;
                }
            }
        }

        m = lend;

L110:
        if (m > lend)
        {
            e[m - 1] = 0.0;
        }
        p = d[l];
        if (m == l)
        {
            goto L130;
        }

        /*        if remaining matrix is 2-by-2, use dlae2 or dlaev2 */
        /*        to compute its eigensystem. */

        if (m == l - 1)
        {
            if (icompz > 0)
            {
                dlaev2_(&d[l - 1], &e[l - 1], &d[l], &rt1, &rt2, &c, &s);
                /* $$$               work( m ) = c */
                /* $$$               work( n-1+m ) = s */
                /* $$$               call dlasr( 'r', 'v', 'f', n, 2, work( m ), */
                /* $$$     $                     work( n-1+m ), z( 1, l-1 ), ldz ) */

                /*               *** New starting with version 2.5 *** */

                tst = z[l];
                z[l] = c * tst - s * z[l - 1];
                z[l - 1] = s * tst + c * z[l - 1];
                /*               ************************************* */
            }
            else
            {
                dlae2_(&d[l - 1], &e[l - 1], &d[l], &rt1, &rt2);
            }
            d[l - 1] = rt1;
            d[l] = rt2;
            e[l - 1] = 0.0;
            l += -2;
            if (l >= lend)
            {
                goto L90;
            }
            goto L140;
        }

        if (jtot == nmaxit)
        {
            goto L140;
        }
        ++jtot;

        /*        form shift. */

        g = (d[l - 1] - p) / (e[l - 1] * 2.);
        r = dlapy2_(&g, &d_one);
        g = d[m] - p + e[l - 1] / (g + d_sign(&r, &g));

        s = 1.0;
        c = 1.0;
        p = 0.0;

        /*        inner loop */

        lm1 = l - 1;
        i__1 = lm1;
        for (i = m; i <= i__1; ++i)
        {
            f = s * e[i];
            b = c * e[i];
            dlartg_(&g, &f, &c, &s, &r);
            if (i != m)
            {
                e[i - 1] = r;
            }
            g = d[i] - p;
            r = (d[i + 1] - g) * s + c * 2. * b;
            p = s * r;
            d[i] = g + p;
            g = c * r - b;

            /*           if eigenvectors are desired, then save rotations. */

            if (icompz > 0)
            {
                work[i] = c;
                work[*n - 1 + i] = s;
            }

        }

        /*        if eigenvectors are desired, then apply saved rotations. */

        if (icompz > 0)
        {
            mm = l - m + 1;
            /* $$$            call dlasr( 'r', 'v', 'f', n, mm, work( m ), work( n-1+m ), */
            /* $$$     $                  z( 1, m ), ldz ) */

            /*           *** New starting with version 2.5 *** */

            dlasr_("r", "v", "f", &c__1, &mm, &work[m], &work[*n - 1 + m], &z[m], &c__1);
            /*           ************************************* */
        }

        d[l] -= p;
        e[lm1] = g;
        goto L90;

        /*        eigenvalue found. */

L130:
        d[l] = p;

        --l;
        if (l >= lend)
        {
            goto L90;
        }
        goto L140;
    }

    /*     undo scaling if necessary */

L140:
    if (iscale == 1)
    {
        i__1 = lendsv - lsv + 1;
        dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d[lsv], n, info);
        i__1 = lendsv - lsv;
        dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &e[lsv], n, info);
    }
    else if (iscale == 2)
    {
        i__1 = lendsv - lsv + 1;
        dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d[lsv], n, info);
        i__1 = lendsv - lsv;
        dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &e[lsv], n, info);
    }

    /*     check for no convergence to an eigenvalue after a total */
    /*     of n*maxit iterations. */

    if (jtot < nmaxit)
    {
        goto L10;
    }
    i__1 = *n - 1;
    for (i = 1; i <= i__1; ++i)
    {
        if (e[i] != 0.0)
        {
            ++(*info);
        }
    }
    goto L190;

    /*     order eigenvalues and eigenvectors. */

L160:
    if (icompz == 0)
    {
        /*        use quick sort */

        dlasrt_("i", n, &d[1], info);
    }
    else
    {
        /*        use selection sort to minimize swaps of eigenvectors */

        i__1 = *n;
        for (ii = 2; ii <= i__1; ++ii)
        {
            i = ii - 1;
            k = i;
            p = d[i];
            i__2 = *n;
            for (j = ii; j <= i__2; ++j)
            {
                if (d[j] < p)
                {
                    k = j;
                    p = d[j];
                }
            }
            if (k != i)
            {
                d[k] = d[i];
                d[i] = p;
                /* $$$               call dswap( n, z( 1, i ), 1, z( 1, k ), 1 ) */
                /*           *** New starting with version 2.5 *** */

                p = z[k];
                z[k] = z[i];
                z[i] = p;
                /*           ************************************* */
            }

        }
    }

L190:
    return 0;

    /* ------------- */
    /* End of dstqrb */
    /* ------------- */

} /* dstqrb_ */

