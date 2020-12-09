/* arpack-ng\SRC\snaitr.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \BeginDoc
 *
 * \Name: snaitr
 *
 * \Description:
 *  Reverse communication interface for applying NP additional steps to
 *  a K step nonsymmetric Arnoldi factorization.
 *
 *  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T
 *
 *          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.
 *
 *  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T
 *
 *          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.
 *
 *  where OP and B are as in snaupd.  The B-norm of r_{k+p} is also
 *  computed and returned.
 *
 * \Usage:
 *  call snaitr
 *     ( IDO, BMAT, N, K, NP, NB, RESID, RNORM, V, LDV, H, LDH,
 *       IPNTR, WORKD, INFO )
 *
 * \Arguments
 *  IDO     Integer.  (INPUT/OUTPUT)
 *          Reverse communication flag.
 *          -------------------------------------------------------------
 *          IDO =  0: first call to the reverse communication interface
 *          IDO = -1: compute  Y = OP * X  where
 *                    IPNTR(1) is the pointer into WORK for X,
 *                    IPNTR(2) is the pointer into WORK for Y.
 *                    This is for the restart phase to force the new
 *                    starting vector into the range of OP.
 *          IDO =  1: compute  Y = OP * X  where
 *                    IPNTR(1) is the pointer into WORK for X,
 *                    IPNTR(2) is the pointer into WORK for Y,
 *                    IPNTR(3) is the pointer into WORK for B * X.
 *          IDO =  2: compute  Y = B * X  where
 *                    IPNTR(1) is the pointer into WORK for X,
 *                    IPNTR(2) is the pointer into WORK for Y.
 *          IDO = 99: done
 *          -------------------------------------------------------------
 *          When the routine is used in the "shift-and-invert" mode, the
 *          vector B * Q is already available and do not need to be
 *          recompute in forming OP * Q.
 *
 *  BMAT    Character*1.  (INPUT)
 *          BMAT specifies the type of the matrix B that defines the
 *          semi-inner product for the operator OP.  See snaupd.
 *          B = 'I' -> standard eigenvalue problem A*x = lambda*x
 *          B = 'G' -> generalized eigenvalue problem A*x = lambda*M**x
 *
 *  N       Integer.  (INPUT)
 *          Dimension of the eigenproblem.
 *
 *  K       Integer.  (INPUT)
 *          Current size of V and H.
 *
 *  NP      Integer.  (INPUT)
 *          Number of additional Arnoldi steps to take.
 *
 *  NB      Integer.  (INPUT)
 *          Blocksize to be used in the recurrence.
 *          Only work for NB = 1 right now.  The goal is to have a
 *          program that implement both the block and non-block method.
 *
 *  RESID   Real array of length N.  (INPUT/OUTPUT)
 *          On INPUT:  RESID contains the residual vector r_{k}.
 *          On OUTPUT: RESID contains the residual vector r_{k+p}.
 *
 *  RNORM   Real scalar.  (INPUT/OUTPUT)
 *          B-norm of the starting residual on input.
 *          B-norm of the updated residual r_{k+p} on output.
 *
 *  V       Real N by K+NP array.  (INPUT/OUTPUT)
 *          On INPUT:  V contains the Arnoldi vectors in the first K
 *          columns.
 *          On OUTPUT: V contains the new NP Arnoldi vectors in the next
 *          NP columns.  The first K columns are unchanged.
 *
 *  LDV     Integer.  (INPUT)
 *          Leading dimension of V exactly as declared in the calling
 *          program.
 *
 *  H       Real (K+NP) by (K+NP) array.  (INPUT/OUTPUT)
 *          H is used to store the generated upper Hessenberg matrix.
 *
 *  LDH     Integer.  (INPUT)
 *          Leading dimension of H exactly as declared in the calling
 *          program.
 *
 *  IPNTR   Integer array of length 3.  (OUTPUT)
 *          Pointer to mark the starting locations in the WORK for
 *          vectors used by the Arnoldi iteration.
 *          -------------------------------------------------------------
 *          IPNTR(1): pointer to the current operand vector X.
 *          IPNTR(2): pointer to the current result vector Y.
 *          IPNTR(3): pointer to the vector B * X when used in the
 *                    shift-and-invert mode.  X is the current operand.
 *          -------------------------------------------------------------
 *
 *  WORKD   Real work array of length 3*N.  (REVERSE COMMUNICATION)
 *          Distributed array to be used in the basic Arnoldi iteration
 *          for reverse communication.  The calling program should not
 *          use WORKD as temporary workspace during the iteration !!!!!!
 *          On input, WORKD(1:N) = B*RESID and is used to save some
 *          computation at the first step.
 *
 *  INFO    Integer.  (OUTPUT)
 *          = 0: Normal exit.
 *          > 0: Size of the spanning invariant subspace of OP found.
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
 *     sgetv0  ARPACK routine to generate the initial vector.
 *     ivout   ARPACK utility routine that prints integers.
 *     arscnd  ARPACK utility routine for timing.
 *     smout   ARPACK utility routine that prints matrices
 *     svout   ARPACK utility routine that prints vectors.
 *     slabad  LAPACK routine that computes machine constants.
 *     slamch  LAPACK routine that determines machine constants.
 *     slascl  LAPACK routine for careful scaling of a matrix.
 *     slanhs  LAPACK routine that computes various norms of a matrix.
 *     sgemv   Level 2 BLAS routine for matrix vector multiplication.
 *     saxpy   Level 1 BLAS that computes a vector triad.
 *     sscal   Level 1 BLAS that scales a vector.
 *     scopy   Level 1 BLAS that copies one vector to another .
 *     sdot    Level 1 BLAS that computes the scalar product of two vectors.
 *     snrm2   Level 1 BLAS that computes the norm of a vector.
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
 * FILE: naitr.F   SID: 2.4   DATE OF SID: 8/27/96   RELEASE: 2
 *
 * \Remarks
 *  The algorithm implemented is:
 *
 *  restart = .false.
 *  Given V_{k} = [v_{1}, ..., v_{k}], r_{k};
 *  r_{k} contains the initial residual vector even for k = 0;
 *  Also assume that rnorm = || B*r_{k} || and B*r_{k} are already
 *  computed by the calling program.
 *
 *  betaj = rnorm ; p_{k+1} = B*r_{k} ;
 *  For  j = k+1, ..., k+np  Do
 *     1) if ( betaj < tol ) stop or restart depending on j.
 *        ( At present tol is zero )
 *        if ( restart ) generate a new starting vector.
 *     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];
 *        p_{j} = p_{j}/betaj
 *     3) r_{j} = OP*v_{j} where OP is defined as in snaupd
 *        For shift-invert mode p_{j} = B*v_{j} is already available.
 *        wnorm = || OP*v_{j} ||
 *     4) Compute the j-th step residual vector.
 *        w_{j} =  V_{j}^T * B * OP * v_{j}
 *        r_{j} =  OP*v_{j} - V_{j} * w_{j}
 *        H(:,j) = w_{j};
 *        H(j,j-1) = rnorm
 *        rnorm = || r_(j) ||
 *        If (rnorm > 0.717*wnorm) accept step and go back to 1)
 *     5) Re-orthogonalization step:
 *        s = V_{j}'*B*r_{j}
 *        r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} ||
 *        alphaj = alphaj + s_{j};
 *     6) Iterative refinement step:
 *        If (rnorm1 > 0.717*rnorm) then
 *           rnorm = rnorm1
 *           accept step and go back to 1)
 *        Else
 *           rnorm = rnorm1
 *           If this is the first time in step 6), go to 5)
 *           Else r_{j} lies in the span of V_{j} numerically.
 *              Set r_{j} = 0 and rnorm = 0; go to 1)
 *        EndIf
 *  End Do
 *
 * \EndLib
 */
int snaitr_(int *ido, char *bmat, int *n, int *k,int *np, int *nb,
            float *resid, float *rnorm, float *v, int *ldv, float *h, int *ldh,
            int *ipntr, float *workd, int *info)
{
    /* Initialized data */

    static bool first = true;

    /* System generated locals */
    int h_dim1, h_offset, v_dim1, v_offset, i__1, i__2;
    float r__1, r__2;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    int i;
    static int j;
    static float t0, t1, t2, t3, t4, t5;
    int jj;
    static int ipj, irj, ivj;
    static float ulp;
    float tst1;
    static int ierr, iter;
    static float unfl, ovfl;
    static int itry;
    float temp1;
    static bool orth1, orth2, step3, step4;
    static float betaj;
    int infol;
    float xtemp[2];
    static float wnorm;
    static float rnorm1;
    static bool rstart;
    static int msglvl;
    static float smlnum;


    /* Parameter adjustments */
    --workd;
    --resid;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h -= h_offset;
    --ipntr;

    /* Function Body */

    if (first)
    {
        /* --------------------------------------- */
        /* Set machine-dependent constants for the */
        /* the splitting and deflation criterion.  */
        /* If norm(H) <= sqrt(OVFL),               */
        /* overflow should not occur.              */
        /* REFERENCE: LAPACK subroutine slahqr     */
        /* --------------------------------------- */

        unfl = slamch_("S");
        ovfl = 1.f / unfl;
        slabad_(&unfl, &ovfl);
        ulp = slamch_("P");
        smlnum = unfl * (*n / ulp);
        first = false;
    }

    if (*ido == 0)
    {
        /* ----------------------------- */
        /* Initialize timing statistics  */
        /* & message level for debugging */
        /* ----------------------------- */

#ifndef NO_TIMER
        arscnd_(&t0);
#endif

        msglvl = debug_1.mnaitr;

        /* ---------------------------- */
        /* Initial call to this routine */
        /* ---------------------------- */

        *info = 0;
        step3 = false;
        step4 = false;
        rstart = false;
        orth1 = false;
        orth2 = false;
        j = *k + 1;
        ipj = 1;
        irj = ipj + *n;
        ivj = irj + *n;
    }

    /* ----------------------------------------------- */
    /* When in reverse communication mode one of:      */
    /* STEP3, STEP4, ORTH1, ORTH2, RSTART              */
    /* will be .true. when ....                        */
    /* STEP3: return from computing OP*v_{j}.          */
    /* STEP4: return from computing B-norm of OP*v_{j} */
    /* ORTH1: return from computing B-norm of r_{j+1}  */
    /* ORTH2: return from computing B-norm of          */
    /*        correction to the residual vector.       */
    /* RSTART: return from OP computations needed by   */
    /*         sgetv0.                                 */
    /* ----------------------------------------------- */

    if (step3)
    {
        goto L50;
    }
    if (step4)
    {
        goto L60;
    }
    if (orth1)
    {
        goto L70;
    }
    if (orth2)
    {
        goto L90;
    }
    if (rstart)
    {
        goto L30;
    }

    /* --------------------------- */
    /* Else this is the first step */
    /* --------------------------- */

    /* ------------------------------------------------------------ */
    /*                                                              */
    /*        A R N O L D I     I T E R A T I O N     L O O P       */
    /*                                                              */
    /* Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) */
    /* ------------------------------------------------------------ */
L1000:

#ifndef NO_TRACE
    if (msglvl > 1)
    {
        ivout_(&c__1, &j, &debug_1.ndigit, "_naitr: generating Arnoldi vector number");
        svout_(&c__1, rnorm, &debug_1.ndigit, "_naitr: B-norm of the current residual is");
    }
#endif

    /* ------------------------------------------------- */
    /* STEP 1: Check if the B norm of j-th residual      */
    /* vector is zero. Equivalent to determining whether */
    /* an exact j-step Arnoldi factorization is present. */
    /* ------------------------------------------------- */

    betaj = *rnorm;
    if (*rnorm > 0.f)
    {
        goto L40;
    }

    /* ------------------------------------------------- */
    /* Invariant subspace found, generate a new starting */
    /* vector which is orthogonal to the current Arnoldi */
    /* basis and continue the iteration.                 */
    /* ------------------------------------------------- */

#ifndef NO_TRACE
    if (msglvl > 0)
    {
        ivout_(&c__1, &j, &debug_1.ndigit, "_naitr: ****** RESTART AT STEP ******");
    }
#endif

    /* ------------------------------------------- */
    /* ITRY is the loop variable that controls the */
    /* maximum amount of times that a restart is   */
    /* attempted. NRSTRT is used by stat.h         */
    /* ------------------------------------------- */

    betaj = 0.f;
    ++timing_1.nrstrt;
    itry = 1;
L20:
    rstart = true;
    *ido = 0;
L30:

    /* ------------------------------------ */
    /* If in reverse communication mode and */
    /* RSTART = .true. flow returns here.   */
    /* ------------------------------------ */

    sgetv0_(ido, bmat, &itry, &c_false, n, &j, &v[v_offset], ldv, &resid[1], rnorm, &ipntr[1], &workd[1], &ierr);
    if (*ido != 99)
    {
        goto L9000;
    }
    if (ierr < 0)
    {
        ++itry;
        if (itry <= 3)
        {
            goto L20;
        }

        /* ---------------------------------------------- */
        /* Give up after several restart attempts.        */
        /* Set INFO to the size of the invariant subspace */
        /* which spans OP and exit.                       */
        /* ---------------------------------------------- */

        *info = j - 1;
#ifndef NO_TIMER
        arscnd_(&t1);
        timing_1.tnaitr += t1 - t0;
#endif

        *ido = 99;
        goto L9000;
    }

L40:

    /* ------------------------------------------------------- */
    /* STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  */
    /* Note that p_{j} = B*r_{j-1}. In order to avoid overflow */
    /* when reciprocating a small RNORM, test against lower    */
    /* machine bound.                                          */
    /* ------------------------------------------------------- */

    scopy_(n, &resid[1], &c__1, &v[j * v_dim1 + 1], &c__1);
    if (*rnorm >= unfl)
    {
        temp1 = 1.f / *rnorm;
        sscal_(n, &temp1, &v[j * v_dim1 + 1], &c__1);
        sscal_(n, &temp1, &workd[ipj], &c__1);
    }
    else
    {
        /* --------------------------------------- */
        /* To scale both v_{j} and p_{j} carefully */
        /* use LAPACK routine SLASCL               */
        /* --------------------------------------- */

        slascl_("G", &i, &i, rnorm, &s_one, n, &c__1, &v[j * v_dim1 + 1], n, &infol);
        slascl_("G", &i, &i, rnorm, &s_one, n, &c__1, &workd[ipj], n, &infol);
    }

    /* ---------------------------------------------------- */
    /* STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} */
    /* Note that this is not quite yet r_{j}. See STEP 4    */
    /* ---------------------------------------------------- */

    step3 = true;
    ++timing_1.nopx;
#ifndef NO_TIMER
    arscnd_(&t2);
#endif

    scopy_(n, &v[j * v_dim1 + 1], &c__1, &workd[ivj], &c__1);
    ipntr[1] = ivj;
    ipntr[2] = irj;
    ipntr[3] = ipj;
    *ido = 1;

    /* --------------------------------- */
    /* Exit in order to compute OP*v_{j} */
    /* --------------------------------- */

    goto L9000;
L50:

    /* -------------------------------- */
    /* Back from reverse communication; */
    /* WORKD(IRJ:IRJ+N-1) := OP*v_{j}   */
    /* if step3 = .true.                */
    /* -------------------------------- */

#ifndef NO_TIMER
    arscnd_(&t3);
    timing_1.tmvopx += t3 - t2;
#endif

    step3 = false;

    /* ---------------------------------------- */
    /* Put another copy of OP*v_{j} into RESID. */
    /* ---------------------------------------- */

    scopy_(n, &workd[irj], &c__1, &resid[1], &c__1);

    /* ------------------------------------- */
    /* STEP 4:  Finish extending the Arnoldi */
    /*          factorization to length j.   */
    /* ------------------------------------- */

#ifndef NO_TIMER
    arscnd_(&t2);
#endif

    if (*bmat == 'G')
    {
        ++timing_1.nbx;
        step4 = true;
        ipntr[1] = irj;
        ipntr[2] = ipj;
        *ido = 2;

        /* ----------------------------------- */
        /* Exit in order to compute B*OP*v_{j} */
        /* ----------------------------------- */

        goto L9000;
    }
    else if (*bmat == 'I')
    {
        scopy_(n, &resid[1], &c__1, &workd[ipj], &c__1);
    }
L60:

    /* -------------------------------- */
    /* Back from reverse communication; */
    /* WORKD(IPJ:IPJ+N-1) := B*OP*v_{j} */
    /* if step4 = .true.                */
    /* -------------------------------- */

#ifndef NO_TIMER
    if (*bmat == 'G')
    {
        arscnd_(&t3);
        timing_1.tmvbx += t3 - t2;
    }
#endif

    step4 = false;

    /* ----------------------------------- */
    /* The following is needed for STEP 5. */
    /* Compute the B-norm of OP*v_{j}.     */
    /* ----------------------------------- */

    if (*bmat == 'G')
    {
        wnorm = sdot_(n, &resid[1], &c__1, &workd[ipj], &c__1);
        wnorm = sqrt((dabs(wnorm)));
    }
    else if (*bmat == 'I')
    {
        wnorm = snrm2_(n, &resid[1], &c__1);
    }

    /* --------------------------------------- */
    /* Compute the j-th residual corresponding */
    /* to the j step factorization.            */
    /* Use Classical Gram Schmidt and compute: */
    /* w_{j} <-  V_{j}^T * B * OP * v_{j}      */
    /* r_{j} <-  OP*v_{j} - V_{j} * w_{j}      */
    /* --------------------------------------- */

    /* ---------------------------------------- */
    /* Compute the j Fourier coefficients w_{j} */
    /* WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  */
    /* ---------------------------------------- */

    sgemv_("T", n, &j, &s_one, &v[v_offset], ldv, &workd[ipj], &c__1, &s_zero, &h[j * h_dim1 + 1], &c__1);

    /* ------------------------------------ */
    /* Orthogonalize r_{j} against V_{j}.   */
    /* RESID contains OP*v_{j}. See STEP 3. */
    /* ------------------------------------ */

    sgemv_("N", n, &j, &s_m1, &v[v_offset], ldv, &h[j * h_dim1 + 1], &c__1,&s_one, &resid[1], &c__1);

    if (j > 1)
    {
        h[j + (j - 1) * h_dim1] = betaj;
    }

#ifndef NO_TIMER
    arscnd_(&t4);
    arscnd_(&t2);
#endif

    orth1 = true;

    if (*bmat == 'G')
    {
        ++timing_1.nbx;
        scopy_(n, &resid[1], &c__1, &workd[irj], &c__1);
        ipntr[1] = irj;
        ipntr[2] = ipj;
        *ido = 2;

        /* -------------------------------- */
        /* Exit in order to compute B*r_{j} */
        /* -------------------------------- */

        goto L9000;
    }
    else if (*bmat == 'I')
    {
        scopy_(n, &resid[1], &c__1, &workd[ipj], &c__1);
    }
L70:

    /* ------------------------------------------------- */
    /* Back from reverse communication if ORTH1 = .true. */
    /* WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    */
    /* ------------------------------------------------- */

#ifndef NO_TIMER
    if (*bmat == 'G')
    {
        arscnd_(&t3);
        timing_1.tmvbx += t3 - t2;
    }
#endif

    orth1 = false;

    /* ---------------------------- */
    /* Compute the B-norm of r_{j}. */
    /* ---------------------------- */

    if (*bmat == 'G')
    {
        *rnorm = sdot_(n, &resid[1], &c__1, &workd[ipj], &c__1);
        *rnorm = sqrt((dabs(*rnorm)));
    }
    else if (*bmat == 'I')
    {
        *rnorm = snrm2_(n, &resid[1], &c__1);
    }

    /* --------------------------------------------------------- */
    /* STEP 5: Re-orthogonalization / Iterative refinement phase */
    /* Maximum NITER_ITREF tries.                                */
    /*                                                           */
    /*          s      = V_{j}^T * B * r_{j}                     */
    /*          r_{j}  = r_{j} - V_{j}*s                         */
    /*          alphaj = alphaj + s_{j}                          */
    /*                                                           */
    /* The stopping criteria used for iterative refinement is    */
    /* discussed in Parlett's book SEP, page 107 and in Gragg &  */
    /* Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         */
    /* Determine if we need to correct the residual. The goal is */
    /* to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  */
    /* The following test determines whether the sine of the     */
    /* angle between  OP*x and the computed residual is less     */
    /* than or equal to 0.717.                                   */
    /* --------------------------------------------------------- */

    if (*rnorm > wnorm * .717f)
    {
        goto L100;
    }
    iter = 0;
    ++timing_1.nrorth;

    /* ------------------------------------------------- */
    /* Enter the Iterative refinement phase. If further  */
    /* refinement is necessary, loop back here. The loop */
    /* variable is ITER. Perform a step of Classical     */
    /* Gram-Schmidt using all the Arnoldi vectors V_{j}  */
    /* ------------------------------------------------- */

L80:

#ifndef NO_TRACE
    if (msglvl > 2)
    {
        xtemp[0] = wnorm;
        xtemp[1] = *rnorm;
        svout_(&c__2, xtemp, &debug_1.ndigit, "_naitr: re-orthonalization; wnorm and rnorm are");
        svout_(&j, &h[j * h_dim1 + 1], &debug_1.ndigit, "_naitr: j-th column of H");
    }
#endif

    /* -------------------------------------------------- */
    /* Compute V_{j}^T * B * r_{j}.                       */
    /* WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). */
    /* -------------------------------------------------- */

    sgemv_("T", n, &j, &s_one, &v[v_offset], ldv, &workd[ipj], &c__1, &s_zero, &workd[irj], &c__1);

    /* ------------------------------------------- */
    /* Compute the correction to the residual:     */
    /* r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1). */
    /* The correction to H is v(:,1:J)*H(1:J,1:J)  */
    /* + v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j.         */
    /* ------------------------------------------- */

    sgemv_("N", n, &j, &s_m1, &v[v_offset], ldv, &workd[irj], &c__1, &s_one, &resid[1], &c__1);
    saxpy_(&j, &s_one, &workd[irj], &c__1, &h[j * h_dim1 + 1], &c__1);

    orth2 = true;
#ifndef NO_TIMER
    arscnd_(&t2);
#endif

    if (*bmat == 'G')
    {
        ++timing_1.nbx;
        scopy_(n, &resid[1], &c__1, &workd[irj], &c__1);
        ipntr[1] = irj;
        ipntr[2] = ipj;
        *ido = 2;

        /* --------------------------------- */
        /* Exit in order to compute B*r_{j}. */
        /* r_{j} is the corrected residual.  */
        /* --------------------------------- */

        goto L9000;
    }
    else if (*bmat == 'I')
    {
        scopy_(n, &resid[1], &c__1, &workd[ipj], &c__1);
    }
L90:

    /* ------------------------------------------------- */
    /* Back from reverse communication if ORTH2 = .true. */
    /* ------------------------------------------------- */

#ifndef NO_TIMER
    if (*bmat == 'G')
    {
        arscnd_(&t3);
        timing_1.tmvbx += t3 - t2;
    }
#endif

    /* --------------------------------------------------- */
    /* Compute the B-norm of the corrected residual r_{j}. */
    /* --------------------------------------------------- */

    if (*bmat == 'G')
    {
        rnorm1 = sdot_(n, &resid[1], &c__1, &workd[ipj], &c__1);
        rnorm1 = sqrt((dabs(rnorm1)));
    }
    else if (*bmat == 'I')
    {
        rnorm1 = snrm2_(n, &resid[1], &c__1);
    }

#ifndef NO_TRACE
    if (msglvl > 0 && iter > 0)
    {
        ivout_(&c__1, &j, &debug_1.ndigit, "_naitr: Iterative refinement for Arnoldi residual");
        if (msglvl > 2)
        {
            xtemp[0] = *rnorm;
            xtemp[1] = rnorm1;
            svout_(&c__2, xtemp, &debug_1.ndigit, "_naitr: iterative refinement ; rnorm and rnorm1 are");
        }
    }
#endif

    /* --------------------------------------- */
    /* Determine if we need to perform another */
    /* step of re-orthogonalization.           */
    /* --------------------------------------- */

    if (rnorm1 > *rnorm * .717f)
    {
        /* ------------------------------------- */
        /* No need for further refinement.       */
        /* The cosine of the angle between the   */
        /* corrected residual vector and the old */
        /* residual vector is greater than 0.717 */
        /* In other words the corrected residual */
        /* and the old residual vector share an  */
        /* angle of less than arcCOS(0.717)      */
        /* ------------------------------------- */

        *rnorm = rnorm1;
    }
    else
    {
        /* ----------------------------------------- */
        /* Another step of iterative refinement step */
        /* is required. NITREF is used by stat.h     */
        /* ----------------------------------------- */

        ++timing_1.nitref;
        *rnorm = rnorm1;
        ++iter;
        if (iter <= 1)
        {
            goto L80;
        }

        /* ----------------------------------------------- */
        /* Otherwise RESID is numerically in the span of V */
        /* ----------------------------------------------- */

        i__1 = *n;
        for (jj = 1; jj <= i__1; ++jj)
        {
            resid[jj] = 0.f;

        }
        *rnorm = 0.f;
    }

    /* -------------------------------------------- */
    /* Branch here directly if iterative refinement */
    /* wasn't necessary or after at most NITER_REF  */
    /* steps of iterative refinement.               */
    /* -------------------------------------------- */

L100:

    rstart = false;
    orth2 = false;

#ifndef NO_TIMER
    arscnd_(&t5);
    timing_1.titref += t5 - t4;
#endif

    /* ---------------------------------- */
    /* STEP 6: Update  j = j+1;  Continue */
    /* ---------------------------------- */

    ++j;
    if (j > *k + *np)
    {
#ifndef NO_TIMER
        arscnd_(&t1);
        timing_1.tnaitr += t1 - t0;
#endif

        *ido = 99;
        i__1 = *k + *np - 1;
        for (i = max(1,*k); i <= i__1; ++i)
        {
            /* ------------------------------------------ */
            /* Check for splitting and deflation.         */
            /* Use a standard test as in the QR algorithm */
            /* REFERENCE: LAPACK subroutine slahqr        */
            /* ------------------------------------------ */

            tst1 = (r__1 = h[i + i * h_dim1], dabs(r__1)) + (r__2 = h[
                        i + 1 + (i + 1) * h_dim1], dabs(r__2));
            if (tst1 == 0.f)
            {
                i__2 = *k + *np;
                tst1 = slanhs_("1", &i__2, &h[h_offset], ldh, &workd[*n + 1]);
            }
            /* Computing MAX */
            r__2 = ulp * tst1;
            if ((r__1 = h[i + 1 + i * h_dim1], dabs(r__1)) <= dmax(r__2,
                    smlnum))
            {
                h[i + 1 + i * h_dim1] = 0.f;
            }

        }

#ifndef NO_TRACE
        if (msglvl > 2)
        {
            i__1 = *k + *np;
            i__2 = *k + *np;
            smout_(&i__1, &i__2, &h[h_offset], ldh, &debug_1.ndigit, "_naitr: Final upper Hessenberg matrix H of order K+NP");
        }
#endif

        goto L9000;
    }

    /* ------------------------------------------------------ */
    /* Loop back to extend the factorization by another step. */
    /* ------------------------------------------------------ */

    goto L1000;

    /* ------------------------------------------------------------- */
    /*                                                               */
    /*  E N D     O F     M A I N     I T E R A T I O N     L O O P  */
    /*                                                               */
    /* ------------------------------------------------------------- */

L9000:
    return 0;

    /* ------------- */
    /* End of snaitr */
    /* ------------- */

} /* snaitr_ */

