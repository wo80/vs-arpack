/* arpack-ng\SRC\snapps.f -- translated by f2c (version 20100827). */

#include "arpack_internal.h"

/**
 * \BeginDoc
 *
 * \Name: snapps
 *
 * \Description:
 *  Given the Arnoldi factorization
 *
 *     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
 *
 *  apply NP implicit shifts resulting in
 *
 *     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
 *
 *  where Q is an orthogonal matrix which is the product of rotations
 *  and reflections resulting from the NP bulge change sweeps.
 *  The updated Arnoldi factorization becomes:
 *
 *     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
 *
 * \Usage:
 *  call snapps
 *     ( N, KEV, NP, SHIFTR, SHIFTI, V, LDV, H, LDH, RESID, Q, LDQ,
 *       WORKL, WORKD )
 *
 * \Arguments
 *  N       Integer.  (INPUT)
 *          Problem size, i.e. size of matrix A.
 *
 *  KEV     Integer.  (INPUT/OUTPUT)
 *          KEV+NP is the size of the input matrix H.
 *          KEV is the size of the updated matrix HNEW.  KEV is only
 *          updated on output when fewer than NP shifts are applied in
 *          order to keep the conjugate pair together.
 *
 *  NP      Integer.  (INPUT)
 *          Number of implicit shifts to be applied.
 *
 *  SHIFTR, Real array of length NP.  (INPUT)
 *  SHIFTI  Real and imaginary part of the shifts to be applied.
 *          Upon, entry to snapps, the shifts must be sorted so that the
 *          conjugate pairs are in consecutive locations.
 *
 *  V       Real N by (KEV+NP) array.  (INPUT/OUTPUT)
 *          On INPUT, V contains the current KEV+NP Arnoldi vectors.
 *          On OUTPUT, V contains the updated KEV Arnoldi vectors
 *          in the first KEV columns of V.
 *
 *  LDV     Integer.  (INPUT)
 *          Leading dimension of V exactly as declared in the calling
 *          program.
 *
 *  H       Real (KEV+NP) by (KEV+NP) array.  (INPUT/OUTPUT)
 *          On INPUT, H contains the current KEV+NP by KEV+NP upper
 *          Hessenber matrix of the Arnoldi factorization.
 *          On OUTPUT, H contains the updated KEV by KEV upper Hessenberg
 *          matrix in the KEV leading submatrix.
 *
 *  LDH     Integer.  (INPUT)
 *          Leading dimension of H exactly as declared in the calling
 *          program.
 *
 *  RESID   Real array of length N.  (INPUT/OUTPUT)
 *          On INPUT, RESID contains the residual vector r_{k+p}.
 *          On OUTPUT, RESID is the update residual vector rnew_{k}
 *          in the first KEV locations.
 *
 *  Q       Real KEV+NP by KEV+NP work array.  (WORKSPACE)
 *          Work array used to accumulate the rotations and reflections
 *          during the bulge chase sweep.
 *
 *  LDQ     Integer.  (INPUT)
 *          Leading dimension of Q exactly as declared in the calling
 *          program.
 *
 *  WORKL   Real work array of length (KEV+NP).  (WORKSPACE)
 *          Private (replicated) array on each PE or array allocated on
 *          the front end.
 *
 *  WORKD   Real work array of length 2*N.  (WORKSPACE)
 *          Distributed array used in the application of the accumulated
 *          orthogonal matrix Q.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Local variables:
 *     xxxxxx  real
 *
 * \References:
 *  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
 *     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
 *     pp 357-385.
 *
 * \Routines called:
 *     ivout   ARPACK utility routine that prints integers.
 *     arscnd  ARPACK utility routine for timing.
 *     smout   ARPACK utility routine that prints matrices.
 *     svout   ARPACK utility routine that prints vectors.
 *     slabad  LAPACK routine that computes machine constants.
 *     slacpy  LAPACK matrix copy routine.
 *     slamch  LAPACK routine that determines machine constants.
 *     slanhs  LAPACK routine that computes various norms of a matrix.
 *     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     slarf   LAPACK routine that applies Householder reflection to
 *             a matrix.
 *     slarfg  LAPACK Householder reflection construction routine.
 *     slartg  LAPACK Givens rotation construction routine.
 *     slaset  LAPACK matrix initialization routine.
 *     sgemv   Level 2 BLAS routine for matrix vector multiplication.
 *     saxpy   Level 1 BLAS that computes a vector triad.
 *     scopy   Level 1 BLAS that copies one vector to another .
 *     sscal   Level 1 BLAS that scales a vector.
 *
 * \Author
 *     Danny Sorensen               Phuong Vu
 *     Richard Lehoucq              CRPC / Rice University
 *     Dept. of Computational &     Houston, Texas
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * \Revision history:
 *     xx/xx/92: Version ' 2.4'
 *
 * \SCCS Information: @(#)
 * FILE: napps.F   SID: 2.4   DATE OF SID: 3/28/97   RELEASE: 2
 *
 * \Remarks
 *  1. In this version, each shift is applied to all the sublocks of
 *     the Hessenberg matrix H and not just to the submatrix that it
 *     comes from. Deflation as in LAPACK routine slahqr (QR algorithm
 *     for upper Hessenberg matrices ) is used.
 *     The subdiagonals of H are enforced to be non-negative.
 *
 * \EndLib
 */
int snapps_(int *n, int *kev, int *np, float *
            shiftr, float *shifti, float *v, int *ldv, float *h, int *ldh,
            float *resid, float *q, int *ldq, float *workl, float *workd)
{
    /* Initialized data */

    static bool first = true;

    /* System generated locals */
    int h_dim, h_offset, v_dim, v_offset, q_dim, q_offset, i__1, i__2,
            i__3, i__4;
    float r__1, r__2;

    /* Local variables */
    float c, f, g;
    int i, j;
    float r, s, t, u[3];
    static float t0, t1;
    float h11, h12, h21, h22, h32;
    int jj, ir, nr;
    float tau;
    static float ulp;
    float tst1;
    int iend;
    static float unfl, ovfl;
    bool cconj;
    float sigmai;
    float sigmar;
    int istart, kplusp, msglvl;
    static float smlnum;

    /* Parameter adjustments */
    v_dim = *ldv;
    v_offset = 1 + v_dim;
    v -= v_offset;
    h_dim = *ldh;
    h_offset = 1 + h_dim;
    h -= h_offset;
    q_dim = *ldq;
    q_offset = 1 + q_dim;
    q -= q_offset;

    /* Function Body */

    if (first)
    {
        /* --------------------------------------------- */
        /* Set machine-dependent constants for the       */
        /* stopping criterion. If norm(H) <= sqrt(OVFL), */
        /* overflow should not occur.                    */
        /* REFERENCE: LAPACK subroutine slahqr           */
        /* --------------------------------------------- */

        unfl = slamch_("S");
        ovfl = 1.0f / unfl;
        slabad_(&unfl, &ovfl);
        ulp = slamch_("P");
        smlnum = unfl * (*n / ulp);
        first = false;
    }

    /* ----------------------------- */
    /* Initialize timing statistics  */
    /* & message level for debugging */
    /* ----------------------------- */

#ifndef NO_TIMER
    arscnd_(&t0);
#endif

    msglvl = debug_1.mnapps;
    kplusp = *kev + *np;

    /* ------------------------------------------ */
    /* Initialize Q to the identity to accumulate */
    /* the rotations and reflections              */
    /* ------------------------------------------ */

    slaset_("A", &kplusp, &kplusp, &s_zero, &s_one, &q[q_offset], ldq);

    /* -------------------------------------------- */
    /* Quick return if there are no shifts to apply */
    /* -------------------------------------------- */

    if (*np == 0)
    {
        goto L9000;
    }

    /* -------------------------------------------- */
    /* Chase the bulge with the application of each */
    /* implicit shift. Each shift is applied to the */
    /* whole matrix including each block.           */
    /* -------------------------------------------- */

    cconj = false;
    i__1 = *np;
    for (jj = 1; jj <= i__1; ++jj)
    {
        /* TODO: fix jj - 1 */
        sigmar = shiftr[jj - 1];
        sigmai = shifti[jj - 1];

#ifndef NO_TRACE
        if (msglvl > 2)
        {
            ivout_(&c__1, &jj, &debug_1.ndigit, "_napps: shift number.");
            svout_(&c__1, &sigmar, &debug_1.ndigit, "_napps: The float part of the shift ");
            svout_(&c__1, &sigmai, &debug_1.ndigit, "_napps: The imaginary part of the shift ");
        }
#endif

        /* ----------------------------------------------- */
        /* The following set of conditionals is necessary  */
        /* in order that complex conjugate pairs of shifts */
        /* are applied together or not at all.             */
        /* ----------------------------------------------- */

        if (cconj)
        {
            /* --------------------------------------- */
            /* cconj = .true. means the previous shift */
            /* had non-zero imaginary part.            */
            /* --------------------------------------- */

            cconj = false;
            goto L110;
        }
        else if (jj < *np && dabs(sigmai) > 0.0f)
        {
            /* ---------------------------------- */
            /* Start of a complex conjugate pair. */
            /* ---------------------------------- */

            cconj = true;
        }
        else if (jj == *np && dabs(sigmai) > 0.0f)
        {
            /* -------------------------------------------- */
            /* The last shift has a nonzero imaginary part. */
            /* Don't apply it; thus the order of the        */
            /* compressed H is order KEV+1 since only np-1  */
            /* were applied.                                */
            /* -------------------------------------------- */

            ++(*kev);
            goto L110;
        }
        istart = 1;
L20:

        /* ------------------------------------------------ */
        /* if sigmai = 0 then                               */
        /*    Apply the jj-th shift ...                     */
        /* else                                             */
        /*    Apply the jj-th and (jj+1)-th together ...    */
        /*    (Note that jj < np at this point in the code) */
        /* end                                              */
        /* to the current block of H. The next do loop      */
        /* determines the current block ;                   */
        /* ------------------------------------------------ */

        i__2 = kplusp - 1;
        for (i = istart; i <= i__2; ++i)
        {
            /* -------------------------------------- */
            /* Check for splitting and deflation. Use */
            /* a standard test as in the QR algorithm */
            /* REFERENCE: LAPACK subroutine slahqr    */
            /* -------------------------------------- */

            r__1 = h[i + i * h_dim];
            r__2 = h[i + 1 + (i + 1) * h_dim];
            tst1 = dabs(r__1) + dabs(r__2);
            if (tst1 == 0.0f)
            {
                i__3 = kplusp - jj + 1;
                tst1 = slanhs_("1", &i__3, &h[h_offset], ldh, workl);
            }
            /* Computing MAX */
            r__1 = h[i + 1 + i * h_dim];
            r__2 = ulp * tst1;
            if (dabs(r__1) <= dmax(r__2, smlnum))
            {
#ifndef NO_TRACE
                if (msglvl > 0)
                {
                    ivout_(&c__1, &i, &debug_1.ndigit, "_napps: matrix splitting at row/column no.");
                    ivout_(&c__1, &jj, &debug_1.ndigit, "_napps: matrix splitting with shift number.");
                    svout_(&c__1, &h[i + 1 + i * h_dim], &debug_1.ndigit, "_napps: off diagonal element.");
                }
#endif

                iend = i;
                h[i + 1 + i * h_dim] = 0.0f;
                goto L40;
            }
        }
        iend = kplusp;
L40:

#ifndef NO_TRACE
        if (msglvl > 2)
        {
            ivout_(&c__1, &istart, &debug_1.ndigit, "_napps: Start of current block ");
            ivout_(&c__1, &iend, &debug_1.ndigit, "_napps: End of current block ");
        }
#endif

        /* ---------------------------------------------- */
        /* No reason to apply a shift to block of order 1 */
        /* ---------------------------------------------- */

        if (istart == iend)
        {
            goto L100;
        }

        /* ---------------------------------------------------- */
        /* If istart + 1 = iend then no reason to apply a       */
        /* complex conjugate pair of shifts on a 2 by 2 matrix. */
        /* ---------------------------------------------------- */

        if (istart + 1 == iend && dabs(sigmai) > 0.0f)
        {
            goto L100;
        }

        h11 = h[istart + istart * h_dim];
        h21 = h[istart + 1 + istart * h_dim];
        if (dabs(sigmai) <= 0.0f)
        {
            /* ------------------------------------------- */
            /* Real-valued shift ==> apply single shift QR */
            /* ------------------------------------------- */

            f = h11 - sigmar;
            g = h21;

            i__2 = iend - 1;
            for (i = istart; i <= i__2; ++i)
            {
                /* --------------------------------------------------- */
                /* Construct the plane rotation G to zero out the bulge */
                /* --------------------------------------------------- */

                slartg_(&f, &g, &c, &s, &r);
                if (i > istart)
                {
                    /* ----------------------------------------- */
                    /* The following ensures that h(1:iend-1,1), */
                    /* the first iend-2 off diagonal of elements */
                    /* H, remain non negative.                   */
                    /* ----------------------------------------- */

                    if (r < 0.0f)
                    {
                        r = -r;
                        c = -c;
                        s = -s;
                    }
                    h[i + (i - 1) * h_dim] = r;
                    h[i + 1 + (i - 1) * h_dim] = 0.0f;
                }

                /* ------------------------------------------- */
                /* Apply rotation to the left of H;  H <- G'*H */
                /* ------------------------------------------- */

                i__3 = kplusp;
                for (j = i; j <= i__3; ++j)
                {
                    t = c * h[i + j * h_dim] + s * h[i + 1 + j * h_dim];
                    h[i + 1 + j * h_dim] = -s * h[i + j * h_dim] + c * h[i + 1 + j * h_dim];
                    h[i + j * h_dim] = t;
                }

                /* ------------------------------------------- */
                /* Apply rotation to the right of H;  H <- H*G */
                /* ------------------------------------------- */

                /* Computing MIN */
                i__4 = i + 2;
                i__3 = min(i__4,iend);
                for (j = 1; j <= i__3; ++j)
                {
                    t = c * h[j + i * h_dim] + s * h[j + (i + 1) * h_dim];
                    h[j + (i + 1) * h_dim] = -s * h[j + i * h_dim] + c * h[j + (i + 1) * h_dim];
                    h[j + i * h_dim] = t;
                }

                /* -------------------------------------------------- */
                /* Accumulate the rotation in the matrix Q;  Q <- Q*G */
                /* -------------------------------------------------- */

                /* Computing MIN */
                i__4 = i + jj;
                i__3 = min(i__4,kplusp);
                for (j = 1; j <= i__3; ++j)
                {
                    t = c * q[j + i * q_dim] + s * q[j + (i + 1) * q_dim];
                    q[j + (i + 1) * q_dim] = -s * q[j + i * q_dim] + c * q[j + (i + 1) * q_dim];
                    q[j + i * q_dim] = t;
                }

                /* ------------------------- */
                /* Prepare for next rotation */
                /* ------------------------- */

                if (i < iend - 1)
                {
                    f = h[i + 1 + i * h_dim];
                    g = h[i + 2 + i * h_dim];
                }
            }

            /* --------------------------------- */
            /* Finished applying the real shift. */
            /* --------------------------------- */

        }
        else
        {
            /* -------------------------------------------------- */
            /* Complex conjugate shifts ==> apply double shift QR */
            /* -------------------------------------------------- */

            h12 = h[istart + (istart + 1) * h_dim];
            h22 = h[istart + 1 + (istart + 1) * h_dim];
            h32 = h[istart + 2 + (istart + 1) * h_dim];

            /* ------------------------------------------------------- */
            /* Compute 1st column of (H - shift*I)*(H - conj(shift)*I) */
            /* ------------------------------------------------------- */

            s = sigmar * 2.0f;
            t = slapy2_(&sigmar, &sigmai);
            u[0] = (h11 * (h11 - s) + t * t) / h21 + h12;
            u[1] = h11 + h22 - s;
            u[2] = h32;

            i__2 = iend - 1;
            for (i = istart; i <= i__2; ++i)
            {
                /* Computing MIN */
                i__3 = 3, i__4 = iend - i + 1;
                nr = min(i__3,i__4);

                /* --------------------------------------------------- */
                /* Construct Householder reflector G to zero out u(1). */
                /* G is of the form I - tau*( 1 u )' * ( 1 u' ).       */
                /* --------------------------------------------------- */

                slarfg_(&nr, u, &u[1], &c__1, &tau);

                if (i > istart)
                {
                    h[i + (i - 1) * h_dim] = u[0];
                    h[i + 1 + (i - 1) * h_dim] = 0.0f;
                    if (i < iend - 1)
                    {
                        h[i + 2 + (i - 1) * h_dim] = 0.0f;
                    }
                }
                u[0] = 1.0f;

                /* ------------------------------------ */
                /* Apply the reflector to the left of H */
                /* ------------------------------------ */

                i__3 = kplusp - i + 1;
                slarf_("L", &nr, &i__3, u, &c__1, &tau, &h[i + i * h_dim], ldh, workl);

                /* ------------------------------------- */
                /* Apply the reflector to the right of H */
                /* ------------------------------------- */

                /* Computing MIN */
                i__3 = i + 3;
                ir = min(i__3,iend);
                slarf_("R", &ir, &nr, u, &c__1, &tau, &h[i * h_dim + 1], ldh, workl);

                /* --------------------------------------------------- */
                /* Accumulate the reflector in the matrix Q;  Q <- Q*G */
                /* --------------------------------------------------- */

                slarf_("R", &kplusp, &nr, u, &c__1, &tau, &q[i * q_dim + 1], ldq, workl);

                /* -------------------------- */
                /* Prepare for next reflector */
                /* -------------------------- */

                if (i < iend - 1)
                {
                    u[0] = h[i + 1 + i * h_dim];
                    u[1] = h[i + 2 + i * h_dim];
                    if (i < iend - 2)
                    {
                        u[2] = h[i + 3 + i * h_dim];
                    }
                }
            }

            /* ------------------------------------------ */
            /* Finished applying a complex pair of shifts */
            /* to the current block                       */
            /* ------------------------------------------ */

        }

L100:

        /* ------------------------------------------------------- */
        /* Apply the same shift to the next block if there is any. */
        /* ------------------------------------------------------- */

        istart = iend + 1;
        if (iend < kplusp)
        {
            goto L20;
        }

        /* ------------------------------------------- */
        /* Loop back to the top to get the next shift. */
        /* ------------------------------------------- */

L110:
        ;
    }

    /* ------------------------------------------------ */
    /* Perform a similarity transformation that makes   */
    /* sure that H will have non negative sub diagonals */
    /* ------------------------------------------------ */

    i__1 = *kev;
    for (j = 1; j <= i__1; ++j)
    {
        if (h[j + 1 + j * h_dim] < 0.0f)
        {
            i__2 = kplusp - j + 1;
            sscal_(&i__2, &s_m1, &h[j + 1 + j * h_dim], ldh);
            /* Computing MIN */
            i__3 = j + 2;
            i__2 = min(i__3,kplusp);
            sscal_(&i__2, &s_m1, &h[(j + 1) * h_dim + 1], &c__1);
            /* Computing MIN */
            i__3 = j + *np + 1;
            i__2 = min(i__3,kplusp);
            sscal_(&i__2, &s_m1, &q[(j + 1) * q_dim + 1], &c__1);
        }
    }

    i__1 = *kev;
    for (i = 1; i <= i__1; ++i)
    {
        /* ------------------------------------------ */
        /* Final check for splitting and deflation.   */
        /* Use a standard test as in the QR algorithm */
        /* REFERENCE: LAPACK subroutine slahqr        */
        /* ------------------------------------------ */

        r__1 = h[i + i * h_dim];
        r__2 = h[i + 1 + (i + 1) * h_dim];
        tst1 = dabs(r__1) + dabs(r__2);
        if (tst1 == 0.0f)
        {
            tst1 = slanhs_("1", kev, &h[h_offset], ldh, workl);
        }
        /* Computing MAX */
        r__1 = ulp * tst1;
        if (h[i + 1 + i * h_dim] <= dmax(r__1,smlnum))
        {
            h[i + 1 + i * h_dim] = 0.0f;
        }
    }

    /* ----------------------------------------------- */
    /* Compute the (kev+1)-st column of (V*Q) and      */
    /* temporarily store the result in WORKD(N+1:2*N). */
    /* This is needed in the residual update since we  */
    /* cannot GUARANTEE that the corresponding entry   */
    /* of H would be zero as in exact arithmetic.      */
    /* ----------------------------------------------- */

    if (h[*kev + 1 + *kev * h_dim] > 0.0f)
    {
        sgemv_("N", n, &kplusp, &s_one, &v[v_offset], ldv, &q[(*kev + 1) * q_dim + 1], &c__1, &s_zero, &workd[*n], &c__1);
    }

    /* -------------------------------------------------------- */
    /* Compute column 1 to kev of (V*Q) in backward order       */
    /* taking advantage of the upper Hessenberg structure of Q. */
    /* -------------------------------------------------------- */

    i__1 = *kev;
    for (i = 1; i <= i__1; ++i)
    {
        i__2 = kplusp - i + 1;
        sgemv_("N", n, &i__2, &s_one, &v[v_offset], ldv, &q[(*kev - i + 1) * q_dim + 1], &c__1, &s_zero, workd, &c__1);
        scopy_(n, workd, &c__1, &v[(kplusp - i + 1) * v_dim + 1], &c__1);
    }

    /* ----------------------------------------------- */
    /*  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). */
    /* ----------------------------------------------- */

    slacpy_("A", n, kev, &v[(kplusp - *kev + 1) * v_dim + 1], ldv, &v[v_offset], ldv);

    /* ------------------------------------------------------------ */
    /* Copy the (kev+1)-st column of (V*Q) in the appropriate place */
    /* ------------------------------------------------------------ */

    if (h[*kev + 1 + *kev * h_dim] > 0.0f)
    {
        scopy_(n, &workd[*n], &c__1, &v[(*kev + 1) * v_dim + 1], &c__1);
    }

    /* ----------------------------------- */
    /* Update the residual vector:         */
    /*    r <- sigmak*r + betak*v(:,kev+1) */
    /* where                               */
    /*    sigmak = (e_{kplusp}'*Q)*e_{kev} */
    /*    betak = e_{kev+1}'*H*e_{kev}     */
    /* ----------------------------------- */

    sscal_(n, &q[kplusp + *kev * q_dim], resid, &c__1);
    if (h[*kev + 1 + *kev * h_dim] > 0.0f)
    {
        saxpy_(n, &h[*kev + 1 + *kev * h_dim], &v[(*kev + 1) * v_dim + 1],&c__1, resid, &c__1);
    }

#ifndef NO_TRACE
    if (msglvl > 1)
    {
        svout_(&c__1, &q[kplusp + *kev * q_dim], &debug_1.ndigit, "_napps: sigmak = (e_{kev+p}^T*Q)*e_{kev}");
        svout_(&c__1, &h[*kev + 1 + *kev * h_dim], &debug_1.ndigit, "_napps: betak = e_{kev+1}^T*H*e_{kev}");
        ivout_(&c__1, kev, &debug_1.ndigit, "_napps: Order of the final Hessenberg matrix ");
        if (msglvl > 2)
        {
            smout_(kev, kev, &h[h_offset], ldh, &debug_1.ndigit, "_napps: updated Hessenberg matrix H for next iteration");
        }
    }
#endif

L9000:
#ifndef NO_TIMER
    arscnd_(&t1);
    timing_1.tnapps += t1 - t0;
#endif

    return 0;

    /* ------------- */
    /* End of snapps */
    /* ------------- */

} /* snapps_ */

