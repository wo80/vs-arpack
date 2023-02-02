/* arpack-ng\SRC\sstqrb.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \BeginDoc
 *
 * \Name: sstqrb
 *
 * \Description:
 *  Computes all eigenvalues and the last component of the eigenvectors
 *  of a symmetric tridiagonal matrix using the implicit QL or QR method.
 *
 *  This is mostly a modification of the LAPACK routine ssteqr.
 *  See Remarks.
 *
 * \Usage:
 *  call sstqrb
 *     ( N, D, E, Z, WORK, INFO )
 *
 * \Arguments
 *  N       Integer.  (INPUT)
 *          The number of rows and columns in the matrix.  N >= 0.
 *
 *  D       Real array, dimension (N).  (INPUT/OUTPUT)
 *          On entry, D contains the diagonal elements of the
 *          tridiagonal matrix.
 *          On exit, D contains the eigenvalues, in ascending order.
 *          If an error exit is made, the eigenvalues are correct
 *          for indices 1,2,...,INFO-1, but they are unordered and
 *          may not be the smallest eigenvalues of the matrix.
 *
 *  E       Real array, dimension (N-1).  (INPUT/OUTPUT)
 *          On entry, E contains the subdiagonal elements of the
 *          tridiagonal matrix in positions 1 through N-1.
 *          On exit, E has been destroyed.
 *
 *  Z       Real array, dimension (N).  (OUTPUT)
 *          On exit, Z contains the last row of the orthonormal
 *          eigenvector matrix of the symmetric tridiagonal matrix.
 *          If an error exit is made, Z contains the last row of the
 *          eigenvector matrix associated with the stored eigenvalues.
 *
 *  WORK    Real array, dimension (max(1,2*N-2)).  (WORKSPACE)
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
 *     saxpy   Level 1 BLAS that computes a vector triad.
 *     scopy   Level 1 BLAS that copies one vector to another.
 *     sswap   Level 1 BLAS that swaps the contents of two vectors.
 *     lsame   LAPACK character comparison routine.
 *     slae2   LAPACK routine that computes the eigenvalues of a 2-by-2
 *             symmetric matrix.
 *     slaev2  LAPACK routine that eigendecomposition of a 2-by-2 symmetric
 *             matrix.
 *     slamch  LAPACK routine that determines machine constants.
 *     slanst  LAPACK routine that computes the norm of a matrix.
 *     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     slartg  LAPACK Givens rotation construction routine.
 *     slascl  LAPACK routine for careful scaling of a matrix.
 *     slaset  LAPACK matrix initialization routine.
 *     slasr   LAPACK routine that applies an orthogonal transformation to
 *             a matrix.
 *     slasrt  LAPACK sorting routine.
 *     ssteqr  LAPACK routine that computes eigenvalues and eigenvectors
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
 *        only commented out and new lines inserted.
 *        All lines commented out have "c$$$" at the beginning.
 *        Note that the LAPACK version 1.0 subroutine SSTEQR contained
 *        bugs.
 *
 * \EndLib
 */
int sstqrb_(int *n, float *d, float *e, float *z, float *
            work, int *info)
{
    /* System generated locals */
    int i__1, i__2;
    float r__1, r__2;

    /* Builtin functions */
    double sqrt(double), r_sign(float *, float *);

    /* Local variables */
    float b, c, f, g;
    int i, j, k, l, m;
    float p, r, s;
    int l1, ii, mm, lm1, mm1, nm1;
    float rt1, rt2, eps;
    int lsv;
    float tst, eps2;
    int lend, jtot;
    float anorm;
    int lendm1, lendp1;
    int iscale;
    float safmin, safmax;
    int lendsv;
    float ssfmin;
    int nmaxit, icompz;
    float ssfmax;

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
            z[1] = 1.0f;
        }
        return 0;
    }

    /*     determine the unit roundoff and over/underflow thresholds. */

    eps = slamch_("E");
    /* Computing 2nd power */
    r__1 = eps;
    eps2 = r__1 * r__1;
    safmin = slamch_("S");
    safmax = 1.0f / safmin;
    ssfmax = sqrt(safmax) / 3.0f;
    ssfmin = sqrt(safmin) / eps2;

    /*     compute the eigenvalues and eigenvectors of the tridiagonal */
    /*     matrix. */

    /* $$      if( icompz.eq.2 ) */
    /* $$$     $   call slaset( 'full', n, n, zero, one, z, ldz ) */

    /*     *** New starting with version 2.5 *** */

    if (icompz == 2)
    {
        i__1 = *n - 1;
        for (j = 1; j <= i__1; ++j)
        {
            z[j] = 0.0f;

        }
        z[*n] = 1.0f;
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
        e[l1 - 1] = 0.0f;
    }
    if (l1 <= nm1)
    {
        for (m = l1; m <= nm1; ++m)
        {
            r__1 = e[m];
            tst = dabs(r__1);
            if (tst == 0.0f)
            {
                goto L30;
            }
            r__1 = d[m];
            r__2 = d[m + 1];
            if (tst <= sqrt(dabs(r__1)) * sqrt(dabs(r__2)) * eps)
            {
                e[m] = 0.0f;
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
    anorm = slanst_("i", &i__1, &d[l], &e[l]);
    iscale = 0;
    if (anorm == 0.0f)
    {
        goto L10;
    }
    if (anorm > ssfmax)
    {
        iscale = 1;
        i__1 = lend - l + 1;
        slascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d[l], n, info);
        i__1 = lend - l;
        slascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, info);
    }
    else if (anorm < ssfmin)
    {
        iscale = 2;
        i__1 = lend - l + 1;
        slascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d[l], n, info);
        i__1 = lend - l;
        slascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, info);
    }

    /*     choose between ql and qr iteration */

    if ((r__1 = d[lend], dabs(r__1)) < (r__2 = d[l], dabs(r__2)))
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
            for (m = l; m <= lendm1; ++m)
            {
                /* Computing 2nd power */
                r__1 = e[m];
                r__2 = dabs(r__1);
                tst = r__2 * r__2;
                r__1 = d[m];
                r__2 = d[m + 1];
                if (tst <= eps2 * dabs(r__1) * dabs(r__2) + safmin)
                {
                    goto L60;
                }
            }
        }

        m = lend;

L60:
        if (m < lend)
        {
            e[m] = 0.0f;
        }
        p = d[l];
        if (m == l)
        {
            goto L80;
        }

        /*        if remaining matrix is 2-by-2, use slae2 or slaev2 */
        /*        to compute its eigensystem. */

        if (m == l + 1)
        {
            if (icompz > 0)
            {
                slaev2_(&d[l], &e[l], &d[l + 1], &rt1, &rt2, &c, &s);
                work[l] = c;
                work[*n - 1 + l] = s;
                /* $$$               call slasr( 'r', 'v', 'b', n, 2, work( l ), */
                /* $$$     $                     work( n-1+l ), z( 1, l ), ldz ) */

                /*              *** New starting with version 2.5 *** */

                tst = z[l + 1];
                z[l + 1] = c * tst - s * z[l];
                z[l] = s * tst + c * z[l];
                /*              ************************************* */
            }
            else
            {
                slae2_(&d[l], &e[l], &d[l + 1], &rt1, &rt2);
            }
            d[l] = rt1;
            d[l + 1] = rt2;
            e[l] = 0.0f;
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

        g = (d[l + 1] - p) / (e[l] * 2.0f);
        r = slapy2_(&g, &s_one);
        g = d[m] - p + e[l] / (g + r_sign(&r, &g));

        s = 1.0f;
        c = 1.0f;
        p = 0.0f;

        /*        inner loop */

        mm1 = m - 1;
        for (i = mm1; i >= l; --i)
        {
            f = s * e[i];
            b = c * e[i];
            slartg_(&g, &f, &c, &s, &r);
            if (i != m - 1)
            {
                e[i + 1] = r;
            }
            g = d[i + 1] - p;
            r = (d[i] - g) * s + c * 2.0f * b;
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
            /* $$$            call slasr( 'r', 'v', 'b', n, mm, work( l ), work( n-1+l ), */
            /* $$$     $                  z( 1, l ), ldz ) */

            /*             *** New starting with version 2.5 *** */

            slasr_("r", "v", "b", &c__1, &mm, &work[l], &work[*n - 1 + l], &z[l], &c__1);
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
            for (m = l; m >= lendp1; --m)
            {
                /* Computing 2nd power */
                r__1 = e[m - 1];
                r__2 = dabs(r__1);
                tst = r__2 * r__2;
                r__1 = d[m];
                r__2 = d[m - 1];
                if (tst <= eps2 * dabs(r__1) * dabs(r__2) + safmin)
                {
                    goto L110;
                }
            }
        }

        m = lend;

L110:
        if (m > lend)
        {
            e[m - 1] = 0.0f;
        }
        p = d[l];
        if (m == l)
        {
            goto L130;
        }

        /*        if remaining matrix is 2-by-2, use slae2 or slaev2 */
        /*        to compute its eigensystem. */

        if (m == l - 1)
        {
            if (icompz > 0)
            {
                slaev2_(&d[l - 1], &e[l - 1], &d[l], &rt1, &rt2, &c, &s);
                /* $$$               work( m ) = c */
                /* $$$               work( n-1+m ) = s */
                /* $$$               call slasr( 'r', 'v', 'f', n, 2, work( m ), */
                /* $$$     $                     work( n-1+m ), z( 1, l-1 ), ldz ) */

                /*               *** New starting with version 2.5 *** */

                tst = z[l];
                z[l] = c * tst - s * z[l - 1];
                z[l - 1] = s * tst + c * z[l - 1];
                /*               ************************************* */
            }
            else
            {
                slae2_(&d[l - 1], &e[l - 1], &d[l], &rt1, &rt2);
            }
            d[l - 1] = rt1;
            d[l] = rt2;
            e[l - 1] = 0.0f;
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

        g = (d[l - 1] - p) / (e[l - 1] * 2.0f);
        r = slapy2_(&g, &s_one);
        g = d[m] - p + e[l - 1] / (g + r_sign(&r, &g));

        s = 1.0f;
        c = 1.0f;
        p = 0.0f;

        /*        inner loop */

        lm1 = l - 1;
        for (i = m; i <= lm1; ++i)
        {
            f = s * e[i];
            b = c * e[i];
            slartg_(&g, &f, &c, &s, &r);
            if (i != m)
            {
                e[i - 1] = r;
            }
            g = d[i] - p;
            r = (d[i + 1] - g) * s + c * 2.0f * b;
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
            /* $$$            call slasr( 'r', 'v', 'f', n, mm, work( m ), work( n-1+m ), */
            /* $$$     $                  z( 1, m ), ldz ) */

            /*           *** New starting with version 2.5 *** */

            slasr_("r", "v", "f", &c__1, &mm, &work[m], &work[*n - 1 + m], &z[m], &c__1);
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
        slascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d[lsv], n, info);
        i__1 = lendsv - lsv;
        slascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &e[lsv], n, info);
    }
    else if (iscale == 2)
    {
        i__1 = lendsv - lsv + 1;
        slascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d[lsv], n, info);
        i__1 = lendsv - lsv;
        slascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &e[lsv], n, info);
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
        if (e[i] != 0.0f)
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

        slasrt_("i", n, &d[1], info);
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
                /* $$$               call sswap( n, z( 1, i ), 1, z( 1, k ), 1 ) */
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
    /* End of sstqrb */
    /* ------------- */

} /* sstqrb_ */

