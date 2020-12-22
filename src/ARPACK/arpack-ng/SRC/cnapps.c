/* arpack-ng\SRC\cnapps.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \BeginDoc
 *
 * \Name: cnapps
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
 *  call cnapps
 *     ( N, KEV, NP, SHIFT, V, LDV, H, LDH, RESID, Q, LDQ,
 *       WORKL, WORKD )
 *
 * \Arguments
 *  N       Integer.  (INPUT)
 *          Problem size, i.e. size of matrix A.
 *
 *  KEV     Integer.  (INPUT/OUTPUT)
 *          KEV+NP is the size of the input matrix H.
 *          KEV is the size of the updated matrix HNEW.
 *
 *  NP      Integer.  (INPUT)
 *          Number of implicit shifts to be applied.
 *
 *  SHIFT   Complex array of length NP.  (INPUT)
 *          The shifts to be applied.
 *
 *  V       Complex N by (KEV+NP) array.  (INPUT/OUTPUT)
 *          On INPUT, V contains the current KEV+NP Arnoldi vectors.
 *          On OUTPUT, V contains the updated KEV Arnoldi vectors
 *          in the first KEV columns of V.
 *
 *  LDV     Integer.  (INPUT)
 *          Leading dimension of V exactly as declared in the calling
 *          program.
 *
 *  H       Complex (KEV+NP) by (KEV+NP) array.  (INPUT/OUTPUT)
 *          On INPUT, H contains the current KEV+NP by KEV+NP upper
 *          Hessenberg matrix of the Arnoldi factorization.
 *          On OUTPUT, H contains the updated KEV by KEV upper Hessenberg
 *          matrix in the KEV leading submatrix.
 *
 *  LDH     Integer.  (INPUT)
 *          Leading dimension of H exactly as declared in the calling
 *          program.
 *
 *  RESID   Complex array of length N.  (INPUT/OUTPUT)
 *          On INPUT, RESID contains the the residual vector r_{k+p}.
 *          On OUTPUT, RESID is the update residual vector rnew_{k}
 *          in the first KEV locations.
 *
 *  Q       Complex KEV+NP by KEV+NP work array.  (WORKSPACE)
 *          Work array used to accumulate the rotations and reflections
 *          during the bulge chase sweep.
 *
 *  LDQ     Integer.  (INPUT)
 *          Leading dimension of Q exactly as declared in the calling
 *          program.
 *
 *  WORKL   Complex work array of length (KEV+NP).  (WORKSPACE)
 *          Private (replicated) array on each PE or array allocated on
 *          the front end.
 *
 *  WORKD   Complex work array of length 2*N.  (WORKSPACE)
 *          Distributed array used in the application of the accumulated
 *          orthogonal matrix Q.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Local variables:
 *     xxxxxx  Complex
 *
 * \References:
 *  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
 *     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
 *     pp 357-385.
 *
 * \Routines called:
 *     ivout   ARPACK utility routine that prints integers.
 *     arscnd  ARPACK utility routine for timing.
 *     cmout   ARPACK utility routine that prints matrices
 *     cvout   ARPACK utility routine that prints vectors.
 *     clacpy  LAPACK matrix copy routine.
 *     clanhs  LAPACK routine that computes various norms of a matrix.
 *     clartg  LAPACK Givens rotation construction routine.
 *     claset  LAPACK matrix initialization routine.
 *     slabad  LAPACK routine for defining the underflow and overflow
 *             limits.
 *     slamch  LAPACK routine that determines machine constants.
 *     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     cgemv   Level 2 BLAS routine for matrix vector multiplication.
 *     caxpy   Level 1 BLAS that computes a vector triad.
 *     ccopy   Level 1 BLAS that copies one vector to another.
 *     cscal   Level 1 BLAS that scales a vector.
 *
 * \Author
 *     Danny Sorensen               Phuong Vu
 *     Richard Lehoucq              CRPC / Rice University
 *     Dept. of Computational &     Houston, Texas
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * \SCCS Information: @(#)
 * FILE: napps.F   SID: 2.3   DATE OF SID: 3/28/97   RELEASE: 2
 *
 * \Remarks
 *  1. In this version, each shift is applied to all the sublocks of
 *     the Hessenberg matrix H and not just to the submatrix that it
 *     comes from. Deflation as in LAPACK routine clahqr (QR algorithm
 *     for upper Hessenberg matrices ) is used.
 *     Upon output, the subdiagonals of H are enforced to be non-negative
 *     real numbers.
 *
 * \EndLib
 */
int cnapps_(int *n, int *kev, int *np, complex *
            shift, complex *v, int *ldv, complex *h, int *ldh, complex *
            resid, complex *q, int *ldq, complex *workl, complex *workd)
{
    /* Initialized data */

    static bool first = true;

    /* System generated locals */
    int h_dim, h_offset, v_dim, v_offset, q_dim, q_offset, i__1, i__2,
            i__3, i__4, i__5, i__6;
    float r__1, r__2, r__3, r__4;
    complex q__1, q__2, q__3, q__4, q__5;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);

    /* Local variables */
    float c;
    complex f, g;
    int i, j;
    complex r, s, t;
    static float t0, t1;
    complex h11, h21;
    int jj;
    static float ulp;
    float tst1;
    int iend;
    static float unfl, ovfl;
    complex sigma;
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
        /* REFERENCE: LAPACK subroutine clahqr           */
        /* --------------------------------------------- */

        unfl = slamch_("S");
        q__1.r = 1.0f / unfl, q__1.i = 0.0f / unfl;
        ovfl = q__1.r;
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

    msglvl = debug_1.mcapps;

    kplusp = *kev + *np;

    /* ------------------------------------------ */
    /* Initialize Q to the identity to accumulate */
    /* the rotations and reflections              */
    /* ------------------------------------------ */

    claset_("A", &kplusp, &kplusp, &c_zero, &c_one, &q[q_offset], ldq);

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

    i__1 = *np;
    for (jj = 1; jj <= i__1; ++jj)
    {
        /* TODO: fix jj - 1 */
        sigma.r = shift[jj - 1].r, sigma.i = shift[jj - 1].i;

#ifndef NO_TRACE
        if (msglvl > 2)
        {
            ivout_(&c__1, &jj, &debug_1.ndigit, "_napps: shift number.");
            cvout_(&c__1, &sigma, &debug_1.ndigit, "_napps: Value of the shift ");
        }
#endif

        istart = 1;
L20:

        i__2 = kplusp - 1;
        for (i = istart; i <= i__2; ++i)
        {
            /* -------------------------------------- */
            /* Check for splitting and deflation. Use */
            /* a standard test as in the QR algorithm */
            /* REFERENCE: LAPACK subroutine clahqr    */
            /* -------------------------------------- */

            i__3 = i + i * h_dim;
            i__4 = i + 1 + (i + 1) * h_dim;
            r__1 = h[i__3].r;
            r__2 = h[i__3].i;
            r__3 = h[i__4].r;
            r__4 = h[i__4].i;
            tst1 = dabs(r__1) + dabs(r__2) + dabs(r__3) + dabs(r__4);
            if (tst1 == 0.0f)
            {
                i__3 = kplusp - jj + 1;
                tst1 = clanhs_("1", &i__3, &h[h_offset], ldh, workl);
            }
            i__3 = i + 1 + i * h_dim;
            /* Computing MAX */
            r__2 = ulp * tst1;
            if ((r__1 = h[i__3].r, dabs(r__1)) <= dmax(r__2,smlnum))
            {
#ifndef NO_TRACE
                if (msglvl > 0)
                {
                    ivout_(&c__1, &i, &debug_1.ndigit, "_napps: matrix splitting at row/column no.");
                    ivout_(&c__1, &jj, &debug_1.ndigit, "_napps: matrix splitting with shift number.");
                    cvout_(&c__1, &h[i__3], &debug_1.ndigit, "_napps: off diagonal element.");
                }
#endif

                iend = i;
                h[i__3].r = 0.0f, h[i__3].i = 0.0f;
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
        /* or if the current block starts after the point */
        /* of compression since we'll discard this stuff  */
        /* ---------------------------------------------- */

        if (istart == iend || istart > *kev)
        {
            goto L100;
        }

        i__2 = istart + istart * h_dim;
        h11.r = h[i__2].r, h11.i = h[i__2].i;
        i__2 = istart + 1 + istart * h_dim;
        h21.r = h[i__2].r, h21.i = h[i__2].i;
        q__1.r = h11.r - sigma.r, q__1.i = h11.i - sigma.i;
        f.r = q__1.r, f.i = q__1.i;
        g.r = h21.r, g.i = h21.i;

        i__2 = iend - 1;
        for (i = istart; i <= i__2; ++i)
        {
            /* ---------------------------------------------------- */
            /* Construct the plane rotation G to zero out the bulge */
            /* ---------------------------------------------------- */

            clartg_(&f, &g, &c, &s, &r);
            if (i > istart)
            {
                i__3 = i + (i - 1) * h_dim;
                h[i__3].r = r.r, h[i__3].i = r.i;
                i__3 = i + 1 + (i - 1) * h_dim;
                h[i__3].r = 0.0f, h[i__3].i = 0.0f;
            }

            /* ------------------------------------------- */
            /* Apply rotation to the left of H;  H <- G'*H */
            /* ------------------------------------------- */

            for (j = i; j <= kplusp; ++j)
            {
                /* t        =  c*h(i,j) + s*h(i+1,j)        */
                /* h(i+1,j) = -conjg(s)*h(i,j) + c*h(i+1,j) */
                /* h(i,j)   = t                             */

                i__4 = i + j * h_dim;
                i__5 = i + 1 + j * h_dim;

                t.r = c * h[i__4].r + (s.r * h[i__5].r - s.i * h[i__5].i);
                t.i = c * h[i__4].i + (s.r * h[i__5].i + s.i * h[i__5].r);

                h[i__5].r = -(s.r * h[i__4].r + s.i * h[i__4].i) + c * h[i__5].r;
                h[i__5].i = -(s.r * h[i__4].i - s.i * h[i__4].r) + c * h[i__5].i;

                h[i__4].r = t.r;
                h[i__4].i = t.i;
            }

            /* ------------------------------------------- */
            /* Apply rotation to the right of H;  H <- H*G */
            /* ------------------------------------------- */

            /* Computing MIN */
            i__4 = i + 2;
            i__3 = min(i__4,iend);
            for (j = 1; j <= i__3; ++j)
            {
                /* t        =  c*h(j,i) + conjg(s)*h(j,i+1) */
                /* h(j,i+1) = -s*h(j,i) + c*h(j,i+1)        */
                /* h(j,i)   = t                             */

                i__4 = j + i * h_dim;
                i__5 = j + (i + 1) * h_dim;

                t.r = c * h[i__4].r + (s.r * h[i__5].r + s.i * h[i__5].i);
                t.i = c * h[i__4].i + (s.r * h[i__5].i - s.i * h[i__5].r);

                h[i__5].r = -(s.r * h[i__4].r - s.i * h[i__4].i) + c * h[i__5].r;
                h[i__5].i = -(s.r * h[i__4].i + s.i * h[i__4].r) + c * h[i__5].i;

                h[i__4].r = t.r;
                h[i__4].i = t.i;
            }

            /* --------------------------------------------------- */
            /* Accumulate the rotation in the matrix Q;  Q <- Q*G' */
            /* --------------------------------------------------- */

            /* Computing MIN */
            i__4 = i + jj;
            i__3 = min(i__4,kplusp);
            for (j = 1; j <= i__3; ++j)
            {
                /* t        =   c*q(j,i) + conjg(s)*q(j,i+1) */
                /* q(j,i+1) = - s*q(j,i) + c*q(j,i+1)        */
                /* q(j,i)   = t                              */

                i__4 = j + i * q_dim;
                i__5 = j + (i + 1) * q_dim;

                t.r = c * q[i__4].r + (s.r * q[i__5].r + s.i * q[i__5].i);
                t.i = c * q[i__4].i + (s.r * q[i__5].i - s.i * q[i__5].r);

                q[i__5].r = -(s.r * q[i__4].r - s.i * q[i__4].i) + c * q[i__5].r;
                q[i__5].i = -(s.r * q[i__4].i + s.i * q[i__4].r) + c * q[i__5].i;

                q[i__4].r = t.r;
                q[i__4].i = t.i;
            }

            /* ------------------------- */
            /* Prepare for next rotation */
            /* ------------------------- */

            if (i < iend - 1)
            {
                i__3 = i + 1 + i * h_dim;
                f.r = h[i__3].r, f.i = h[i__3].i;
                i__3 = i + 2 + i * h_dim;
                g.r = h[i__3].r, g.i = h[i__3].i;
            }
        }

        /* ----------------------------- */
        /* Finished applying the shift.  */
        /* ----------------------------- */

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

    }

    /* ------------------------------------------------- */
    /* Perform a similarity transformation that makes    */
    /* sure that the compressed H will have non-negative */
    /* real subdiagonal elements.                        */
    /* ------------------------------------------------- */

    i__1 = *kev;
    for (j = 1; j <= i__1; ++j)
    {
        i__2 = j + 1 + j * h_dim;
        if (h[i__2].r < 0.0f || h[i__2].i != 0.0f)
        {
            /* t = h(j+1,j) / slapy2(real(h(j+1,j)),aimag(h(j+1,j))) */
            i__2 = j + 1 + j * h_dim;
            r__1 = slapy2_(&h[i__2].r, &h[i__2].i);
            t.r = h[i__2].r / r__1;
            t.i = h[i__2].i / r__1;

            i__2 = kplusp - j + 1;
            i__3 = j + 1 + j * h_dim;
            r_cnjg(&q__1, &t);
            cscal_(&i__2, &q__1, &h[i__3], ldh);
            /* Computing MIN */
            i__3 = j + 2;
            i__2 = min(i__3,kplusp);
            cscal_(&i__2, &t, &h[(j + 1) * h_dim + 1], &c__1);
            /* Computing MIN */
            i__3 = j + *np + 1;
            i__2 = min(i__3,kplusp);
            cscal_(&i__2, &t, &q[(j + 1) * q_dim + 1], &c__1);

            /* h(j+1,j) = cmplx( real( h(j+1,j) ), rzero ) */
            i__2 = j + 1 + j * h_dim;
            h[i__2].i = 0.0f;
        }
    }

    i__1 = *kev;
    for (i = 1; i <= i__1; ++i)
    {
        /* ------------------------------------------ */
        /* Final check for splitting and deflation.   */
        /* Use a standard test as in the QR algorithm */
        /* REFERENCE: LAPACK subroutine clahqr.       */
        /* Note: Since the subdiagonals of the        */
        /* compressed H are nonnegative real numbers, */
        /* we take advantage of this.                 */
        /* ------------------------------------------ */

        i__2 = i + i * h_dim;
        i__3 = i + 1 + (i + 1) * h_dim;
        r__1 = h[i__2].r;
        r__2 = h[i__2].i;
        r__3 = h[i__3].r;
        r__4 = h[i__3].i;
        tst1 = dabs(r__1) + dabs(r__2) + dabs(r__3) + dabs(r__4);
        if (tst1 == 0.0f)
        {
            tst1 = clanhs_("1", kev, &h[h_offset], ldh, workl);
        }
        i__2 = i + 1 + i * h_dim;
        /* Computing MAX */
        r__1 = ulp * tst1;
        if (h[i__2].r <= dmax(r__1,smlnum))
        {
            h[i__2].r = 0.0f, h[i__2].i = 0.0f;
        }
    }

    /* ----------------------------------------------- */
    /* Compute the (kev+1)-st column of (V*Q) and      */
    /* temporarily store the result in WORKD(N+1:2*N). */
    /* This is needed in the residual update since we  */
    /* cannot GUARANTEE that the corresponding entry   */
    /* of H would be zero as in exact arithmetic.      */
    /* ----------------------------------------------- */

    i__1 = *kev + 1 + *kev * h_dim;
    if (h[i__1].r > 0.0f)
    {
        cgemv_("N", n, &kplusp, &c_one, &v[v_offset], ldv, &q[(*kev + 1) * q_dim + 1], &c__1, &c_zero, &workd[*n], &c__1);
    }

    /* -------------------------------------------------------- */
    /* Compute column 1 to kev of (V*Q) in backward order       */
    /* taking advantage of the upper Hessenberg structure of Q. */
    /* -------------------------------------------------------- */

    i__1 = *kev;
    for (i = 1; i <= i__1; ++i)
    {
        i__2 = kplusp - i + 1;
        cgemv_("N", n, &i__2, &c_one, &v[v_offset], ldv, &q[(*kev - i + 1) * q_dim + 1], &c__1, &c_zero, workd, &c__1);
        ccopy_(n, workd, &c__1, &v[(kplusp - i + 1) * v_dim + 1], &c__1);
    }

    /* ----------------------------------------------- */
    /*  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). */
    /* ----------------------------------------------- */

    clacpy_("A", n, kev, &v[(kplusp - *kev + 1) * v_dim + 1], ldv, &v[v_offset], ldv);

    /* ------------------------------------------------------------ */
    /* Copy the (kev+1)-st column of (V*Q) in the appropriate place */
    /* ------------------------------------------------------------ */

    i__1 = *kev + 1 + *kev * h_dim;
    if (h[i__1].r > 0.0f)
    {
        ccopy_(n, &workd[*n], &c__1, &v[(*kev + 1) * v_dim + 1], &c__1);
    }

    /* ----------------------------------- */
    /* Update the residual vector:         */
    /*    r <- sigmak*r + betak*v(:,kev+1) */
    /* where                               */
    /*    sigmak = (e_{kev+p}'*Q)*e_{kev}  */
    /*    betak = e_{kev+1}'*H*e_{kev}     */
    /* ----------------------------------- */

    cscal_(n, &q[kplusp + *kev * q_dim], resid, &c__1);
    i__1 = *kev + 1 + *kev * h_dim;
    if (h[i__1].r > 0.0f)
    {
        caxpy_(n, &h[i__1], &v[(*kev + 1) * v_dim + 1],&c__1, resid, &c__1);
    }

#ifndef NO_TRACE
    if (msglvl > 1)
    {
        cvout_(&c__1, &q[kplusp + *kev * q_dim], &debug_1.ndigit, "_napps: sigmak = (e_{kev+p}^T*Q)*e_{kev}");
        cvout_(&c__1, &h[i__1], &debug_1.ndigit, "_napps: betak = e_{kev+1}^T*H*e_{kev}");
        ivout_(&c__1, kev, &debug_1.ndigit, "_napps: Order of the final Hessenberg matrix ");
        if (msglvl > 2)
        {
            cmout_(kev, kev, &h[h_offset], ldh, &debug_1.ndigit, "_napps: updated Hessenberg matrix H for next iteration");
        }
    }
#endif

L9000:
#ifndef NO_TIMER
    arscnd_(&t1);
    timing_1.tcapps += t1 - t0;
#endif

    return 0;

    /* ------------- */
    /* End of cnapps */
    /* ------------- */

} /* cnapps_ */

