/* arpack-ng\SRC\dnaupd.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \BeginDoc
 *
 * \Name: dnaupd
 *
 * \Description:
 *  Reverse communication interface for the Implicitly Restarted Arnoldi
 *  iteration. This subroutine computes approximations to a few eigenpairs
 *  of a linear operator "OP" with respect to a semi-inner product defined by
 *  a symmetric positive semi-definite real matrix B. B may be the identity
 *  matrix. NOTE: If the linear operator "OP" is real and symmetric
 *  with respect to the real positive semi-definite symmetric matrix B,
 *  i.e. B*OP = (OP`)*B, then subroutine dsaupd  should be used instead.
 *
 *  The computed approximate eigenvalues are called Ritz values and
 *  the corresponding approximate eigenvectors are called Ritz vectors.
 *
 *  dnaupd  is usually called iteratively to solve one of the
 *  following problems:
 *
 *  Mode 1:  A*x = lambda*x.
 *           ===> OP = A  and  B = I.
 *
 *  Mode 2:  A*x = lambda*M*x, M symmetric positive definite
 *           ===> OP = inv[M]*A  and  B = M.
 *           ===> (If M can be factored see remark 3 below)
 *
 *  Mode 3:  A*x = lambda*M*x, M symmetric semi-definite
 *           ===> OP = Real_Part{ inv[A - sigma*M]*M }  and  B = M.
 *           ===> shift-and-invert mode (in real arithmetic)
 *           If OP*x = amu*x, then
 *           amu = 1/2 * [ 1/(lambda-sigma) + 1/(lambda-conjg(sigma)) ].
 *           Note: If sigma is real, i.e. imaginary part of sigma is zero;
 *                 Real_Part{ inv[A - sigma*M]*M } == inv[A - sigma*M]*M
 *                 amu == 1/(lambda-sigma).
 *
 *  Mode 4:  A*x = lambda*M*x, M symmetric semi-definite
 *           ===> OP = Imaginary_Part{ inv[A - sigma*M]*M }  and  B = M.
 *           ===> shift-and-invert mode (in real arithmetic)
 *           If OP*x = amu*x, then
 *           amu = 1/2i * [ 1/(lambda-sigma) - 1/(lambda-conjg(sigma)) ].
 *
 *  Both mode 3 and 4 give the same enhancement to eigenvalues close to
 *  the (complex) shift sigma.  However, as lambda goes to infinity,
 *  the operator OP in mode 4 dampens the eigenvalues more strongly than
 *  does OP defined in mode 3.
 *
 *  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
 *        should be accomplished either by a direct method
 *        using a sparse matrix factorization and solving
 *
 *           [A - sigma*M]*w = v  or M*w = v,
 *
 *        or through an iterative method for solving these
 *        systems.  If an iterative method is used, the
 *        convergence test must be more stringent than
 *        the accuracy requirements for the eigenvalue
 *        approximations.
 *
 * \Usage:
 *  call dnaupd
 *     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
 *       IPNTR, WORKD, WORKL, LWORKL, INFO )
 *
 * \Arguments
 *  IDO     Integer.  (INPUT/OUTPUT)
 *          Reverse communication flag.  IDO must be zero on the first
 *          call to dnaupd .  IDO will be set internally to
 *          indicate the type of operation to be performed.  Control is
 *          then given back to the calling routine which has the
 *          responsibility to carry out the requested operation and call
 *          dnaupd  with the result.  The operand is given in
 *          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
 *          -------------------------------------------------------------
 *          IDO =  0: first call to the reverse communication interface
 *          IDO = -1: compute  Y = OP * X  where
 *                    IPNTR(1) is the pointer into WORKD for X,
 *                    IPNTR(2) is the pointer into WORKD for Y.
 *                    This is for the initialization phase to force the
 *                    starting vector into the range of OP.
 *          IDO =  1: compute  Y = OP * X  where
 *                    IPNTR(1) is the pointer into WORKD for X,
 *                    IPNTR(2) is the pointer into WORKD for Y.
 *                    In mode 3 and 4, the vector B * X is already
 *                    available in WORKD(ipntr(3)).  It does not
 *                    need to be recomputed in forming OP * X.
 *          IDO =  2: compute  Y = B * X  where
 *                    IPNTR(1) is the pointer into WORKD for X,
 *                    IPNTR(2) is the pointer into WORKD for Y.
 *          IDO =  3: compute the IPARAM(8) real and imaginary parts
 *                    of the shifts where INPTR(14) is the pointer
 *                    into WORKL for placing the shifts. See Remark
 *                    5 below.
 *          IDO = 99: done
 *          -------------------------------------------------------------
 *
 *  BMAT    Character*1.  (INPUT)
 *          BMAT specifies the type of the matrix B that defines the
 *          semi-inner product for the operator OP.
 *          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
 *          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
 *
 *  N       Integer.  (INPUT)
 *          Dimension of the eigenproblem.
 *
 *  WHICH   Character*2.  (INPUT)
 *          'LM' -> want the NEV eigenvalues of largest magnitude.
 *          'SM' -> want the NEV eigenvalues of smallest magnitude.
 *          'LR' -> want the NEV eigenvalues of largest real part.
 *          'SR' -> want the NEV eigenvalues of smallest real part.
 *          'LI' -> want the NEV eigenvalues of largest imaginary part.
 *          'SI' -> want the NEV eigenvalues of smallest imaginary part.
 *
 *  NEV     Integer.  (INPUT)
 *          Number of eigenvalues of OP to be computed. 0 < NEV < N-1.
 *
 *  TOL     Double precision  scalar.  (INPUT/OUTPUT)
 *          Stopping criterion: the relative accuracy of the Ritz value
 *          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I))
 *          where ABS(RITZ(I)) is the magnitude when RITZ(I) is complex.
 *          DEFAULT = DLAMCH ('EPS')  (machine precision as computed
 *                    by the LAPACK auxiliary subroutine DLAMCH ).
 *
 *  RESID   Double precision  array of length N.  (INPUT/OUTPUT)
 *          On INPUT:
 *          If INFO .EQ. 0, a random initial residual vector is used.
 *          If INFO .NE. 0, RESID contains the initial residual vector,
 *                          possibly from a previous run.
 *          On OUTPUT:
 *          RESID contains the final residual vector.
 *
 *  NCV     Integer.  (INPUT)
 *          Number of columns of the matrix V. NCV must satisfy the two
 *          inequalities 2 <= NCV-NEV and NCV <= N.
 *          This will indicate how many Arnoldi vectors are generated
 *          at each iteration.  After the startup phase in which NEV
 *          Arnoldi vectors are generated, the algorithm generates
 *          approximately NCV-NEV Arnoldi vectors at each subsequent update
 *          iteration. Most of the cost in generating each Arnoldi vector is
 *          in the matrix-vector operation OP*x.
 *          NOTE: 2 <= NCV-NEV in order that complex conjugate pairs of Ritz
 *          values are kept together. (See remark 4 below)
 *
 *  V       Double precision  array N by NCV.  (OUTPUT)
 *          Contains the final set of Arnoldi basis vectors.
 *
 *  LDV     Integer.  (INPUT)
 *          Leading dimension of V exactly as declared in the calling program.
 *
 *  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
 *          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
 *          The shifts selected at each iteration are used to restart
 *          the Arnoldi iteration in an implicit fashion.
 *          -------------------------------------------------------------
 *          ISHIFT = 0: the shifts are provided by the user via
 *                      reverse communication.  The real and imaginary
 *                      parts of the NCV eigenvalues of the Hessenberg
 *                      matrix H are returned in the part of the WORKL
 *                      array corresponding to RITZR and RITZI. See remark
 *                      5 below.
 *          ISHIFT = 1: exact shifts with respect to the current
 *                      Hessenberg matrix H.  This is equivalent to
 *                      restarting the iteration with a starting vector
 *                      that is a linear combination of approximate Schur
 *                      vectors associated with the "wanted" Ritz values.
 *          -------------------------------------------------------------
 *
 *          IPARAM(2) = No longer referenced.
 *
 *          IPARAM(3) = MXITER
 *          On INPUT:  maximum number of Arnoldi update iterations allowed.
 *          On OUTPUT: actual number of Arnoldi update iterations taken.
 *
 *          IPARAM(4) = NB: blocksize to be used in the recurrence.
 *          The code currently works only for NB = 1.
 *
 *          IPARAM(5) = NCONV: number of "converged" Ritz values.
 *          This represents the number of Ritz values that satisfy
 *          the convergence criterion.
 *
 *          IPARAM(6) = IUPD
 *          No longer referenced. Implicit restarting is ALWAYS used.
 *
 *          IPARAM(7) = MODE
 *          On INPUT determines what type of eigenproblem is being solved.
 *          Must be 1,2,3,4; See under \Description of dnaupd  for the
 *          four modes available.
 *
 *          IPARAM(8) = NP
 *          When ido = 3 and the user provides shifts through reverse
 *          communication (IPARAM(1)=0), dnaupd  returns NP, the number
 *          of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark
 *          5 below.
 *
 *          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
 *          OUTPUT: NUMOP  = total number of OP*x operations,
 *                  NUMOPB = total number of B*x operations if BMAT='G',
 *                  NUMREO = total number of steps of re-orthogonalization.
 *
 *  IPNTR   Integer array of length 14.  (OUTPUT)
 *          Pointer to mark the starting locations in the WORKD and WORKL
 *          arrays for matrices/vectors used by the Arnoldi iteration.
 *          -------------------------------------------------------------
 *          IPNTR(1): pointer to the current operand vector X in WORKD.
 *          IPNTR(2): pointer to the current result vector Y in WORKD.
 *          IPNTR(3): pointer to the vector B * X in WORKD when used in
 *                    the shift-and-invert mode.
 *          IPNTR(4): pointer to the next available location in WORKL
 *                    that is untouched by the program.
 *          IPNTR(5): pointer to the NCV by NCV upper Hessenberg matrix
 *                    H in WORKL.
 *          IPNTR(6): pointer to the real part of the ritz value array
 *                    RITZR in WORKL.
 *          IPNTR(7): pointer to the imaginary part of the ritz value array
 *                    RITZI in WORKL.
 *          IPNTR(8): pointer to the Ritz estimates in array WORKL associated
 *                    with the Ritz values located in RITZR and RITZI in WORKL.
 *
 *          IPNTR(14): pointer to the NP shifts in WORKL. See Remark 5 below.
 *
 *          Note: IPNTR(9:13) is only referenced by dneupd . See Remark 2 below.
 *
 *          IPNTR(9):  pointer to the real part of the NCV RITZ values of the
 *                     original system.
 *          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of
 *                     the original system.
 *          IPNTR(11): pointer to the NCV corresponding error bounds.
 *          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
 *                     Schur matrix for H.
 *          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
 *                     of the upper Hessenberg matrix H. Only referenced by
 *                     dneupd  if RVEC = .TRUE. See Remark 2 below.
 *          -------------------------------------------------------------
 *
 *  WORKD   Double precision  work array of length 3*N.  (REVERSE COMMUNICATION)
 *          Distributed array to be used in the basic Arnoldi iteration
 *          for reverse communication.  The user should not use WORKD
 *          as temporary workspace during the iteration. Upon termination
 *          WORKD(1:N) contains B*RESID(1:N). If an invariant subspace
 *          associated with the converged Ritz values is desired, see remark
 *          2 below, subroutine dneupd  uses this output.
 *          See Data Distribution Note below.
 *
 *  WORKL   Double precision  work array of length LWORKL.  (OUTPUT/WORKSPACE)
 *          Private (replicated) array on each PE or array allocated on
 *          the front end.  See Data Distribution Note below.
 *
 *  LWORKL  Integer.  (INPUT)
 *          LWORKL must be at least 3*NCV**2 + 6*NCV.
 *
 *  INFO    Integer.  (INPUT/OUTPUT)
 *          If INFO .EQ. 0, a randomly initial residual vector is used.
 *          If INFO .NE. 0, RESID contains the initial residual vector,
 *                          possibly from a previous run.
 *          Error flag on output.
 *          =  0: Normal exit.
 *          =  1: Maximum number of iterations taken.
 *                All possible eigenvalues of OP has been found. IPARAM(5)
 *                returns the number of wanted converged Ritz values.
 *          =  2: No longer an informational error. Deprecated starting
 *                with release 2 of ARPACK.
 *          =  3: No shifts could be applied during a cycle of the
 *                Implicitly restarted Arnoldi iteration. One possibility
 *                is to increase the size of NCV relative to NEV.
 *                See remark 4 below.
 *          = -1: N must be positive.
 *          = -2: NEV must be positive.
 *          = -3: NCV-NEV >= 2 and less than or equal to N.
 *          = -4: The maximum number of Arnoldi update iteration
 *                must be greater than zero.
 *          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
 *          = -6: BMAT must be one of 'I' or 'G'.
 *          = -7: Length of private work array is not sufficient.
 *          = -8: Error return from LAPACK eigenvalue calculation;
 *          = -9: Starting vector is zero.
 *          = -10: IPARAM(7) must be 1,2,3,4.
 *          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
 *          = -12: IPARAM(1) must be equal to 0 or 1.
 *          = -9999: Could not build an Arnoldi factorization.
 *                   IPARAM(5) returns the size of the current Arnoldi
 *                   factorization.
 *
 * \Remarks
 *  1. The computed Ritz values are approximate eigenvalues of OP. The
 *     selection of WHICH should be made with this in mind when
 *     Mode = 3 and 4.  After convergence, approximate eigenvalues of the
 *     original problem may be obtained with the ARPACK subroutine dneupd .
 *
 *  2. If a basis for the invariant subspace corresponding to the converged Ritz
 *     values is needed, the user must call dneupd  immediately following
 *     completion of dnaupd . This is new starting with release 2 of ARPACK.
 *
 *  3. If M can be factored into a Cholesky factorization M = LL`
 *     then Mode = 2 should not be selected.  Instead one should use
 *     Mode = 1 with  OP = inv(L)*A*inv(L`).  Appropriate triangular
 *     linear systems should be solved with L and L` rather
 *     than computing inverses.  After convergence, an approximate
 *     eigenvector z of the original problem is recovered by solving
 *     L`z = x  where x is a Ritz vector of OP.
 *
 *  4. At present there is no a-priori analysis to guide the selection
 *     of NCV relative to NEV.  The only formal requrement is that NCV > NEV + 2.
 *     However, it is recommended that NCV .ge. 2*NEV+1.  If many problems of
 *     the same type are to be solved, one should experiment with increasing
 *     NCV while keeping NEV fixed for a given test problem.  This will
 *     usually decrease the required number of OP*x operations but it
 *     also increases the work and storage required to maintain the orthogonal
 *     basis vectors.  The optimal "cross-over" with respect to CPU time
 *     is problem dependent and must be determined empirically.
 *     See Chapter 8 of Reference 2 for further information.
 *
 *  5. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the
 *     NP = IPARAM(8) real and imaginary parts of the shifts in locations
 *         real part                  imaginary part
 *         -----------------------    --------------
 *     1   WORKL(IPNTR(14))           WORKL(IPNTR(14)+NP)
 *     2   WORKL(IPNTR(14)+1)         WORKL(IPNTR(14)+NP+1)
 *                        .                          .
 *                        .                          .
 *                        .                          .
 *     NP  WORKL(IPNTR(14)+NP-1)      WORKL(IPNTR(14)+2*NP-1).
 *
 *     Only complex conjugate pairs of shifts may be applied and the pairs
 *     must be placed in consecutive locations. The real part of the
 *     eigenvalues of the current upper Hessenberg matrix are located in
 *     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1) and the imaginary part
 *     in WORKL(IPNTR(7)) through WORKL(IPNTR(7)+NCV-1). They are ordered
 *     according to the order defined by WHICH. The complex conjugate
 *     pairs are kept together and the associated Ritz estimates are located in
 *     WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
 *
 * -----------------------------------------------------------------------
 *
 * \Data Distribution Note:
 *
 *  Fortran-D syntax:
 *  ================
 *  Double precision  resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
 *  decompose  d1(n), d2(n,ncv)
 *  align      resid(i) with d1(i)
 *  align      v(i,j)   with d2(i,j)
 *  align      workd(i) with d1(i)     range (1:n)
 *  align      workd(i) with d1(i-n)   range (n+1:2*n)
 *  align      workd(i) with d1(i-2*n) range (2*n+1:3*n)
 *  distribute d1(block), d2(block,:)
 *  replicated workl(lworkl)
 *
 *  Cray MPP syntax:
 *  ===============
 *  Double precision   resid(n), v(ldv,ncv), workd(n,3), workl(lworkl)
 *  shared     resid(block), v(block,:), workd(block,:)
 *  replicated workl(lworkl)
 *
 *  CM2/CM5 syntax:
 *  ==============
 *
 * -----------------------------------------------------------------------
 *
 *     include   'ex-nonsym.doc'
 *
 * -----------------------------------------------------------------------
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
 *  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
 *     Real Matrices", Linear Algebra and its Applications, vol 88/89,
 *     pp 575-595, (1987).
 *
 * \Routines called:
 *     dnaup2   ARPACK routine that implements the Implicitly Restarted
 *             Arnoldi Iteration.
 *     ivout   ARPACK utility routine that prints integers.
 *     arscnd  ARPACK utility routine for timing.
 *     dvout    ARPACK utility routine that prints vectors.
 *     dlamch   LAPACK routine that determines machine constants.
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
 *     12/16/93: Version '1.1'
 *
 * \SCCS Information: @(#)
 * FILE: naupd.F   SID: 2.8   DATE OF SID: 04/10/01   RELEASE: 2
 *
 * \Remarks
 *
 * \EndLib
 */
int dnaupd_(int *ido, char *bmat, int *n, char *which, int *nev, double *tol,
            double *resid, int *ncv,double *v, int *ldv, int *iparam, int *ipntr,
            double *workd, double *workl, int *lworkl, int *info)
{

    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    int j, ierr;
    static float t0, t1;
    static int mode, np, nev0;
    static int ishift, mxiter;
    int iupd = 1; /* not used */

    int ldh = *ncv;
    int ldq = *ncv;
    int ritzr =  ldh * *ncv;
    int ritzi = ritzr + *ncv;
    int bounds = ritzi + *ncv;
    int iq = bounds + *ncv;
    int iw = iq + ldq * *ncv;

    /* Function Body */
    if (*ido == 0)
    {
        /* ----------------------------- */
        /* Initialize timing statistics  */
        /* & message level for debugging */
        /* ----------------------------- */

        dstatn_();
#ifndef NO_TIMER
        arscnd_(&t0);
#endif

        /* -------------- */
        /* Error checking */
        /* -------------- */

        ierr = 0;
        ishift = iparam[0];
        /* levec  = iparam(2) */
        mxiter = iparam[2];
        /* nb     = iparam(4) */

        /* ------------------------------------------ */
        /* Revision 2 performs only implicit restart. */
        /* ------------------------------------------ */

        mode = iparam[6];

        if (*n <= 0)
        {
            ierr = -1;
        }
        else if (*nev <= 0)
        {
            ierr = -2;
        }
        else if (*ncv <= *nev + 1 || *ncv > *n)
        {
            ierr = -3;
        }
        else if (mxiter <= 0)
        {
            ierr = -4;
        }
        else if (strcmp(which, "LM") != 0 && strcmp(which, "SM") != 0 &&
                 strcmp(which, "LR") != 0 && strcmp(which, "SR") != 0 &&
                 strcmp(which, "LI") != 0 && strcmp(which, "SI") != 0)
        {
            ierr = -5;
        }
        else if (*bmat != 'I' && *bmat != 'G')
        {
            ierr = -6;
        }
        else /* if(complicated condition) */
        {
            /* Computing 2nd power */
            i__1 = *ncv;
            if (*lworkl < i__1 * i__1 * 3 + *ncv * 6)
            {
                ierr = -7;
            }
            else if (mode < 1 || mode > 4)
            {
                ierr = -10;
            }
            else if (mode == 1 && *bmat == 'G')
            {
                ierr = -11;
            }
            else if (ishift < 0 || ishift > 1)
            {
                ierr = -12;
            }
        }

        /* ---------- */
        /* Error Exit */
        /* ---------- */

        if (ierr != 0)
        {
            *info = ierr;
            *ido = 99;
            goto L9000;
        }

        /* ---------------------- */
        /* Set default parameters */
        /* ---------------------- */

        if (*tol <= 0.0)
        {
            *tol = dlamch_("E");
        }

        /* -------------------------------------------- */
        /* NP is the number of additional steps to      */
        /* extend the length NEV Lanczos factorization. */
        /* NEV0 is the local variable designating the   */
        /* size of the invariant subspace desired.      */
        /* -------------------------------------------- */

        np = *ncv - *nev;
        nev0 = *nev;

        /* --------------------------- */
        /* Zero out internal workspace */
        /* --------------------------- */

        /* Computing 2nd power */
        i__2 = *ncv;
        i__1 = i__2 * i__2 * 3 + *ncv * 6;
        for (j = 0; j < i__1; ++j)
        {
            workl[j] = 0.0;
        }

        /* ----------------------------------------------------------- */
        /* Pointer into WORKL for address of H, RITZ, BOUNDS, Q        */
        /* etc... and the remaining workspace.                         */
        /* Also update pointer to be used on output.                   */
        /* Memory is laid out as follows:                              */
        /* workl(1:ncv*ncv) := generated Hessenberg matrix             */
        /* workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary        */
        /*                                   parts of ritz values      */
        /* workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds        */
        /* workl(ncv*ncv+3*ncv+1:2*ncv*ncv+3*ncv) := rotation matrix Q */
        /* workl(2*ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) := workspace       */
        /* The final workspace is needed by subroutine dneigh  called  */
        /* by dnaup2 . Subroutine dneigh  calls LAPACK routines for    */
        /* calculating eigenvalues and the last row of the eigenvector */
        /* matrix.                                                     */
        /* ----------------------------------------------------------- */


        /* Computing 2nd power */
        i__1 = *ncv;
        /* TODO: subtract 1 !!! */
        ipntr[3] = 1 + iw + i__1 * i__1 + i__1 * 3;
        ipntr[4] = 1;
        ipntr[5] = 1 + ritzr;
        ipntr[6] = 1 + ritzi;
        ipntr[7] = 1 + bounds;
        ipntr[13] = 1 + iw;
    }

    /* ----------------------------------------------------- */
    /* Carry out the Implicitly restarted Arnoldi Iteration. */
    /* ----------------------------------------------------- */

    dnaup2_(ido, bmat, n, which, &nev0, &np, tol, resid, &mode, &iupd, &ishift, &mxiter, v, ldv, workl, &ldh, &workl[ritzr], &workl[ritzi], &workl[bounds], &workl[iq], &ldq, &workl[iw], ipntr, workd, info);

    /* ------------------------------------------------ */
    /* ido .ne. 99 implies use of reverse communication */
    /* to compute operations involving OP or shifts.    */
    /* ------------------------------------------------ */

    if (*ido == 3)
    {
        iparam[7] = np;
    }
    if (*ido != 99)
    {
        goto L9000;
    }

    iparam[2] = mxiter;
    iparam[4] = np;
    iparam[8] = timing_1.nopx;
    iparam[9] = timing_1.nbx;
    iparam[10] = timing_1.nrorth;

    /* ---------------------------------- */
    /* Exit if there was an informational */
    /* error within dnaup2 .              */
    /* ---------------------------------- */

    if (*info < 0)
    {
        goto L9000;
    }
    if (*info == 2)
    {
        *info = 3;
    }

#ifndef NO_TIMER
    arscnd_(&t1);
    timing_1.tnaupd = t1 - t0;
#endif

#ifndef NO_TRACE
    int msglvl = debug_1.mnaupd;

    if (msglvl > 1)
    {
        ivout_(&c__1, &mxiter, &debug_1.ndigit, "_naupd: Number of update iterations taken");
        ivout_(&c__1, &np, &debug_1.ndigit, "_naupd: Number of wanted \"converged\" Ritz values");
        dvout_(&np, &workl[ritzr], &debug_1.ndigit, "_naupd: Real part of the final Ritz values");
        dvout_(&np, &workl[ritzi], &debug_1.ndigit, "_naupd: Imaginary part of the final Ritz values");
        dvout_(&np, &workl[bounds], &debug_1.ndigit, "_naupd: Associated Ritz estimates");
    }

    if (msglvl > 0)
    {
        printf("\n ============================================= ");
        printf("\n = Nonsymmetric implicit Arnoldi update code = ");
        printf("\n = Version Number:     %s                 = ", ARPACK_VERSION);
        printf("\n = Version Date:       %s            = ", ARPACK_DATE);
        printf("\n ============================================= ");
        printf("\n = Summary of timing statistics              = ");
        printf("\n ============================================= ");

        printf("\n Total number update iterations             =  %5d", mxiter);
        printf("\n Total number of OP*x operations            =  %5d", timing_1.nopx);
        printf("\n Total number of B*x operations             =  %5d", timing_1.nbx);
        printf("\n Total number of reorthogonalization steps  =  %5d", timing_1.nrorth);
        printf("\n Total number of iterative refinement steps =  %5d", timing_1.nitref);
        printf("\n Total number of restart steps              =  %5d", timing_1.nrstrt);
#ifndef NO_TIMER
        printf("\n Total time in user OP*x operation          =  %12f", timing_1.tmvopx);
        printf("\n Total time in user B*x operation           =  %12f", timing_1.tmvbx);
        printf("\n Total time in Arnoldi update routine       =  %12f", timing_1.tnaupd);
        printf("\n Total time in naup2 routine                =  %12f", timing_1.tnaup2);
        printf("\n Total time in basic Arnoldi iteration loop =  %12f", timing_1.tnaitr);
        printf("\n Total time in reorthogonalization phase    =  %12f", timing_1.titref);
        printf("\n Total time in (re)start vector generation  =  %12f", timing_1.tgetv0);
        printf("\n Total time in Hessenberg eig. subproblem   =  %12f", timing_1.tneigh);
        printf("\n Total time in getting the shifts           =  %12f", timing_1.tngets);
        printf("\n Total time in applying the shifts          =  %12f", timing_1.tnapps);
        printf("\n Total time in convergence testing          =  %12f", timing_1.tnconv);
        printf("\n Total time in computing final Ritz vectors =  %12f", timing_1.trvec);
#endif
    }
#endif

L9000:

    return 0;

    /* ------------- */
    /* End of dnaupd */
    /* ------------- */

} /* dnaupd_ */

