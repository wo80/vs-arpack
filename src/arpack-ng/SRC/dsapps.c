/* arpack-ng\SRC\dsapps.f -- translated by f2c (version 20100827). */

#include "arpack_internal.h"

/* Constants */
static int i_one = 1;
static double d_zero = 0.0;
static double d_one  = 1.0;
static double d_m1  = -1.0;

/**
 * \BeginDoc
 *
 * \Name: dsapps
 *
 * \Description:
 *  Given the Arnoldi factorization
 *
 *     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
 *
 *  apply NP shifts implicitly resulting in
 *
 *     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
 *
 *  where Q is an orthogonal matrix of order KEV+NP. Q is the product of
 *  rotations resulting from the NP bulge chasing sweeps.  The updated Arnoldi
 *  factorization becomes:
 *
 *     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
 *
 * \Usage:
 *  call dsapps
 *     ( N, KEV, NP, SHIFT, V, LDV, H, LDH, RESID, Q, LDQ, WORKD )
 *
 * \Arguments
 *  N       Integer.  (INPUT)
 *          Problem size, i.e. dimension of matrix A.
 *
 *  KEV     Integer.  (INPUT)
 *          INPUT: KEV+NP is the size of the input matrix H.
 *          OUTPUT: KEV is the size of the updated matrix HNEW.
 *
 *  NP      Integer.  (INPUT)
 *          Number of implicit shifts to be applied.
 *
 *  SHIFT   Double precision array of length NP.  (INPUT)
 *          The shifts to be applied.
 *
 *  V       Double precision N by (KEV+NP) array.  (INPUT/OUTPUT)
 *          INPUT: V contains the current KEV+NP Arnoldi vectors.
 *          OUTPUT: VNEW = V(1:n,1:KEV); the updated Arnoldi vectors
 *          are in the first KEV columns of V.
 *
 *  LDV     Integer.  (INPUT)
 *          Leading dimension of V exactly as declared in the calling
 *          program.
 *
 *  H       Double precision (KEV+NP) by 2 array.  (INPUT/OUTPUT)
 *          INPUT: H contains the symmetric tridiagonal matrix of the
 *          Arnoldi factorization with the subdiagonal in the 1st column
 *          starting at H(2,1) and the main diagonal in the 2nd column.
 *          OUTPUT: H contains the updated tridiagonal matrix in the
 *          KEV leading submatrix.
 *
 *  LDH     Integer.  (INPUT)
 *          Leading dimension of H exactly as declared in the calling
 *          program.
 *
 *  RESID   Double precision array of length (N).  (INPUT/OUTPUT)
 *          INPUT: RESID contains the residual vector r_{k+p}.
 *          OUTPUT: RESID is the updated residual vector rnew_{k}.
 *
 *  Q       Double precision KEV+NP by KEV+NP work array.  (WORKSPACE)
 *          Work array used to accumulate the rotations during the bulge
 *          chase sweep.
 *
 *  LDQ     Integer.  (INPUT)
 *          Leading dimension of Q exactly as declared in the calling
 *          program.
 *
 *  WORKD   Double precision work array of length 2*N.  (WORKSPACE)
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
 *  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
 *     Restarted Arnoldi Iteration", Rice University Technical Report
 *     TR95-13, Department of Computational and Applied Mathematics.
 *
 * \Routines called:
 *     ivout   ARPACK utility routine that prints integers.
 *     arscnd  ARPACK utility routine for timing.
 *     dvout   ARPACK utility routine that prints vectors.
 *     dlamch  LAPACK routine that determines machine constants.
 *     dlartg  LAPACK Givens rotation construction routine.
 *     dlacpy  LAPACK matrix copy routine.
 *     dlaset  LAPACK matrix initialization routine.
 *     dgemv   Level 2 BLAS routine for matrix vector multiplication.
 *     daxpy   Level 1 BLAS that computes a vector triad.
 *     dcopy   Level 1 BLAS that copies one vector to another.
 *     dscal   Level 1 BLAS that scales a vector.
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
 *     12/16/93: Version ' 2.4'
 *
 * \SCCS Information: @(#)
 * FILE: sapps.F   SID: 2.6   DATE OF SID: 3/28/97   RELEASE: 2
 *
 * \Remarks
 *  1. In this version, each shift is applied to all the subblocks of
 *     the tridiagonal matrix H and not just to the submatrix that it
 *     comes from. This routine assumes that the subdiagonal elements
 *     of H that are stored in h(1:kev+np,1) are nonegative upon input
 *     and enforce this condition upon output. This version incorporates
 *     deflation. See code for documentation.
 *
 * \EndLib
 */
int dsapps_(int *n, int *kev, int *np,
            double *shift, double *v, int *ldv, double *h,
            int *ldh, double *resid, double *q, int *ldq,
            double *workd)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    int h_dim, h_offset, q_dim, q_offset, v_dim, v_offset, i__1, i__2,
            i__3, i__4;
    double d__1, d__2;

    /* Local variables */
    double c, f, g;
    int i, j;
    double r, s, a1, a2, a3, a4;
    static float t0, t1;
    int jj;
    double big;
    int iend, itop;
    static double epsmch;
    int istart, kplusp, msglvl;

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
        epsmch = dlamch_("E");
        first = FALSE_;
    }
    itop = 1;

    /* ----------------------------- */
    /* Initialize timing statistics  */
    /* & message level for debugging */
    /* ----------------------------- */

#ifndef NO_TIMER
    arscnd_(&t0);
#endif

    msglvl = debug_1.msapps;

    kplusp = *kev + *np;

    /* -------------------------------------------- */
    /* Initialize Q to the identity matrix of order */
    /* kplusp used to accumulate the rotations.     */
    /* -------------------------------------------- */

    dlaset_("A", &kplusp, &kplusp, &d_zero, &d_one, &q[q_offset], ldq);

    /* -------------------------------------------- */
    /* Quick return if there are no shifts to apply */
    /* -------------------------------------------- */

    if (*np == 0)
    {
        goto L9000;
    }

    /* -------------------------------------------------------- */
    /* Apply the np shifts implicitly. Apply each shift to the  */
    /* whole matrix and not just to the submatrix from which it */
    /* comes.                                                   */
    /* -------------------------------------------------------- */

    i__1 = *np;
    for (jj = 1; jj <= i__1; ++jj)
    {
        istart = itop;

        /* -------------------------------------------------------- */
        /* Check for splitting and deflation. Currently we consider */
        /* an off-diagonal element h(i+1,1) negligible if           */
        /*         h(i+1,1) .le. epsmch*( |h(i,2)| + |h(i+1,2)| )   */
        /* for i=1:KEV+NP-1.                                        */
        /* If above condition tests true then we set h(i+1,1) = 0.  */
        /* Note that h(1:KEV+NP,1) are assumed to be non negative.  */
        /* -------------------------------------------------------- */

L20:

        /* ---------------------------------------------- */
        /* The following loop exits early if we encounter */
        /* a negligible off diagonal element.             */
        /* ---------------------------------------------- */

        i__2 = kplusp - 1;
        for (i = istart; i <= i__2; ++i)
        {
            d__1 = h[i + (h_dim << 1)];
            d__2 = h[i + 1 + (h_dim << 1)];
            big = abs(d__1) + abs(d__2);
            if (h[i + 1 + h_dim] <= epsmch * big)
            {
#ifndef NO_TRACE
                if (msglvl > 0)
                {
                    ivout_(1, &i, debug_1.ndigit, "_sapps: deflation at row/column no.");
                    ivout_(1, &jj, debug_1.ndigit, "_sapps: occurred before shift number.");
                    dvout_(1, &h[i + 1 + h_dim], debug_1.ndigit, "_sapps: the corresponding off diagonal element");
                }
#endif

                h[i + 1 + h_dim] = 0.0;
                iend = i;
                goto L40;
            }

        }
        iend = kplusp;
L40:

        if (istart < iend)
        {
            /* ------------------------------------------------------ */
            /* Construct the plane rotation G'(istart,istart+1,theta) */
            /* that attempts to drive h(istart+1,1) to zero.          */
            /* ------------------------------------------------------ */

            /* TODO: fix jj - 1 */
            f = h[istart + (h_dim << 1)] - shift[jj - 1];
            g = h[istart + 1 + h_dim];
            dlartg_(&f, &g, &c, &s, &r);

            /* ----------------------------------------------------- */
            /* Apply rotation to the left and right of H;            */
            /* H <- G' * H * G,  where G = G(istart,istart+1,theta). */
            /* This will create a "bulge".                           */
            /* ----------------------------------------------------- */

            a1 = c * h[istart + (h_dim << 1)] + s * h[istart + 1 + h_dim];
            a2 = c * h[istart + 1 + h_dim] + s * h[istart + 1 + (h_dim << 1)];
            a4 = c * h[istart + 1 + (h_dim << 1)] - s * h[istart + 1 + h_dim];
            a3 = c * h[istart + 1 + h_dim] - s * h[istart + (h_dim <<
                                                    1)];
            h[istart + (h_dim << 1)] = c * a1 + s * a2;
            h[istart + 1 + (h_dim << 1)] = c * a4 - s * a3;
            h[istart + 1 + h_dim] = c * a3 + s * a4;

            /* -------------------------------------------------- */
            /* Accumulate the rotation in the matrix Q;  Q <- Q*G */
            /* -------------------------------------------------- */

            /* Computing MIN */
            i__3 = istart + jj;
            i__2 = min(i__3,kplusp);
            for (j = 1; j <= i__2; ++j)
            {
                a1 = c * q[j + istart * q_dim] + s * q[j + (istart + 1) * q_dim];
                q[j + (istart + 1) * q_dim] = -s * q[j + istart * q_dim] + c * q[j + (istart + 1) * q_dim];
                q[j + istart * q_dim] = a1;
            }

            /* -------------------------------------------- */
            /* The following loop chases the bulge created. */
            /* Note that the previous rotation may also be  */
            /* done within the following loop. But it is    */
            /* kept separate to make the distinction among  */
            /* the bulge chasing sweeps and the first plane */
            /* rotation designed to drive h(istart+1,1) to  */
            /* zero.                                        */
            /* -------------------------------------------- */

            i__2 = iend - 1;
            for (i = istart + 1; i <= i__2; ++i)
            {
                /* -------------------------------------------- */
                /* Construct the plane rotation G'(i,i+1,theta) */
                /* that zeros the i-th bulge that was created   */
                /* by G(i-1,i,theta). g represents the bulge.   */
                /* -------------------------------------------- */

                f = h[i + h_dim];
                g = s * h[i + 1 + h_dim];

                /* -------------------------------- */
                /* Final update with G(i-1,i,theta) */
                /* -------------------------------- */

                h[i + 1 + h_dim] = c * h[i + 1 + h_dim];
                dlartg_(&f, &g, &c, &s, &r);

                /* ----------------------------------------- */
                /* The following ensures that h(1:iend-1,1), */
                /* the first iend-2 off diagonal of elements */
                /* H, remain non negative.                   */
                /* ----------------------------------------- */

                if (r < 0.0)
                {
                    r = -r;
                    c = -c;
                    s = -s;
                }

                /* ------------------------------------------ */
                /* Apply rotation to the left and right of H; */
                /* H <- G * H * G',  where G = G(i,i+1,theta) */
                /* ------------------------------------------ */

                h[i + h_dim] = r;

                a1 = c * h[i + (h_dim << 1)] + s * h[i + 1 + h_dim];
                a2 = c * h[i + 1 + h_dim] + s * h[i + 1 + (h_dim << 1)];
                a3 = c * h[i + 1 + h_dim] - s * h[i + (h_dim << 1)];
                a4 = c * h[i + 1 + (h_dim << 1)] - s * h[i + 1 + h_dim];

                h[i + (h_dim << 1)] = c * a1 + s * a2;
                h[i + 1 + (h_dim << 1)] = c * a4 - s * a3;
                h[i + 1 + h_dim] = c * a3 + s * a4;

                /* -------------------------------------------------- */
                /* Accumulate the rotation in the matrix Q;  Q <- Q*G */
                /* -------------------------------------------------- */

                /* Computing MIN */
                i__4 = i + jj;
                i__3 = min(i__4,kplusp);
                for (j = 1; j <= i__3; ++j)
                {
                    a1 = c * q[j + i * q_dim] + s * q[j + (i + 1) *
                                                       q_dim];
                    q[j + (i + 1) * q_dim] = -s * q[j + i * q_dim] +
                                              c * q[j + (i + 1) * q_dim];
                    q[j + i * q_dim] = a1;
                }
            }
        }

        /* ------------------------ */
        /* Update the block pointer */
        /* ------------------------ */

        istart = iend + 1;

        /* ---------------------------------------- */
        /* Make sure that h(iend,1) is non-negative */
        /* If not then set h(iend,1) <-- -h(iend,1) */
        /* and negate the last column of Q.         */
        /* We have effectively carried out a        */
        /* similarity on transformation H           */
        /* ---------------------------------------- */

        if (h[iend + h_dim] < 0.0)
        {
            h[iend + h_dim] = -h[iend + h_dim];
            dscal_(&kplusp, &d_m1, &q[iend * q_dim + 1], &i_one);
        }

        /* ------------------------------------------------------ */
        /* Apply the same shift to the next block if there is any */
        /* ------------------------------------------------------ */

        if (iend < kplusp)
        {
            goto L20;
        }

        /* --------------------------------------------------- */
        /* Check if we can increase the the start of the block */
        /* --------------------------------------------------- */

        i__2 = kplusp - 1;
        for (i = itop; i <= i__2; ++i)
        {
            if (h[i + 1 + h_dim] > 0.0)
            {
                goto L90;
            }
            ++itop;

        }

        /* --------------------------------- */
        /* Finished applying the jj-th shift */
        /* --------------------------------- */

L90:
        ;
    }

    /* ---------------------------------------- */
    /* All shifts have been applied. Check for  */
    /* more possible deflation that might occur */
    /* after the last shift is applied.         */
    /* ---------------------------------------- */

    i__1 = kplusp - 1;
    for (i = itop; i <= i__1; ++i)
    {
        d__1 = h[i + (h_dim << 1)];
        d__2 = h[i + 1 + (h_dim << 1)];
        big = abs(d__1) + abs(d__2);
        if (h[i + 1 + h_dim] <= epsmch * big)
        {
#ifndef NO_TRACE
            if (msglvl > 0)
            {
                ivout_(1, &i, debug_1.ndigit, "_sapps: deflation at row/column no.");
                dvout_(1, &h[i + 1 + h_dim], debug_1.ndigit, "_sapps: the corresponding off diagonal element");
            }
#endif

            h[i + 1 + h_dim] = 0.0;
        }
    }

    /* ----------------------------------------------- */
    /* Compute the (kev+1)-st column of (V*Q) and      */
    /* temporarily store the result in WORKD(N+1:2*N). */
    /* This is not necessary if h(kev+1,1) = 0.        */
    /* ----------------------------------------------- */

    if (h[*kev + 1 + h_dim] > 0.0)
    {
        dgemv_("N", n, &kplusp, &d_one, &v[v_offset], ldv, &q[(*kev + 1) * q_dim + 1], &i_one, &d_zero, &workd[*n], &i_one);
    }

    /* ----------------------------------------------------- */
    /* Compute column 1 to kev of (V*Q) in backward order    */
    /* taking advantage that Q is an upper triangular matrix */
    /* with lower bandwidth np.                              */
    /* Place results in v(:,kplusp-kev:kplusp) temporarily.  */
    /* ----------------------------------------------------- */

    i__1 = *kev;
    for (i = 1; i <= i__1; ++i)
    {
        i__2 = kplusp - i + 1;
        dgemv_("N", n, &i__2, &d_one, &v[v_offset], ldv, &q[(*kev - i + 1) * q_dim + 1], &i_one, &d_zero, workd, &i_one);
        dcopy_(n, workd, &i_one, &v[(kplusp - i + 1) * v_dim + 1], &i_one);
    }

    /* ----------------------------------------------- */
    /*  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). */
    /* ----------------------------------------------- */

    i__1 = *kev;
    for (i = 1; i <= i__1; ++i)
    {
        dcopy_(n, &v[(*np + i) * v_dim + 1], &i_one, &v[i * v_dim + 1], &i_one);
    }

    /* ------------------------------------------ */
    /* Copy the (kev+1)-st column of (V*Q) in the */
    /* appropriate place if h(kev+1,1) .ne. zero. */
    /* ------------------------------------------ */

    if (h[*kev + 1 + h_dim] > 0.0)
    {
        dcopy_(n, &workd[*n], &i_one, &v[(*kev + 1) * v_dim + 1], &i_one);
    }

    /* ----------------------------------- */
    /* Update the residual vector:         */
    /*    r <- sigmak*r + betak*v(:,kev+1) */
    /* where                               */
    /*    sigmak = (e_{kev+p}'*Q)*e_{kev}  */
    /*    betak = e_{kev+1}'*H*e_{kev}     */
    /* ----------------------------------- */

    dscal_(n, &q[kplusp + *kev * q_dim], resid, &i_one);
    if (h[*kev + 1 + h_dim] > 0.0)
    {
        daxpy_(n, &h[*kev + 1 + h_dim], &v[(*kev + 1) * v_dim + 1], &i_one,resid, &i_one);
    }

#ifndef NO_TRACE
    if (msglvl > 1)
    {
        dvout_(1, &q[kplusp + *kev * q_dim], debug_1.ndigit, "_sapps: sigmak of the updated residual vector");
        dvout_(1, &h[*kev + 1 + h_dim], debug_1.ndigit, "_sapps: betak of the updated residual vector");
        dvout_(*kev, &h[(h_dim << 1) + 1], debug_1.ndigit, "_sapps: updated main diagonal of H for next iteration");
        if (*kev > 1)
        {
            i__1 = *kev - 1;
            dvout_(1, &h[h_dim + 2], debug_1.ndigit, "_sapps: updated sub diagonal of H for next iteration");
        }
    }
#endif

#ifndef NO_TIMER
    arscnd_(&t1);
    timing_1.tsapps += t1 - t0;
#endif

L9000:
    return 0;

    /* ------------- */
    /* End of dsapps */
    /* ------------- */

} /* dsapps_ */

