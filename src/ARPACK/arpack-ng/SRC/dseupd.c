/* arpack-ng\SRC\dseupd.f -- translated by f2c (version 20100827). */

#include <math.h>
#include "arpack_internal.h"

/**
 * \BeginDoc
 *
 * \Name: dseupd
 *
 * \Description:
 *
 *  This subroutine returns the converged approximations to eigenvalues
 *  of A*z = lambda*B*z and (optionally):
 *
 *      (1) the corresponding approximate eigenvectors,
 *
 *      (2) an orthonormal (Lanczos) basis for the associated approximate
 *          invariant subspace,
 *
 *      (3) Both.
 *
 *  There is negligible additional cost to obtain eigenvectors.  An orthonormal
 *  (Lanczos) basis is always computed.  There is an additional storage cost
 *  of n*nev if both are requested (in this case a separate array Z must be
 *  supplied).
 *
 *  These quantities are obtained from the Lanczos factorization computed
 *  by DSAUPD  for the linear operator OP prescribed by the MODE selection
 *  (see IPARAM(7) in DSAUPD  documentation.)  DSAUPD  must be called before
 *  this routine is called. These approximate eigenvalues and vectors are
 *  commonly called Ritz values and Ritz vectors respectively.  They are
 *  referred to as such in the comments that follow.   The computed orthonormal
 *  basis for the invariant subspace corresponding to these Ritz values is
 *  referred to as a Lanczos basis.
 *
 *  See documentation in the header of the subroutine DSAUPD  for a definition
 *  of OP as well as other terms and the relation of computed Ritz values
 *  and vectors of OP with respect to the given problem  A*z = lambda*B*z.
 *
 *  The approximate eigenvalues of the original problem are returned in
 *  ascending algebraic order.  The user may elect to call this routine
 *  once for each desired Ritz vector and store it peripherally if desired.
 *  There is also the option of computing a selected set of these vectors
 *  with a single call.
 *
 * \Usage:
 *  call dseupd
 *     ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, BMAT, N, WHICH, NEV, TOL,
 *       RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO )
 *
 *  RVEC    LOGICAL  (INPUT)
 *          Specifies whether Ritz vectors corresponding to the Ritz value
 *          approximations to the eigenproblem A*z = lambda*B*z are computed.
 *
 *             RVEC = .FALSE.     Compute Ritz values only.
 *
 *             RVEC = .TRUE.      Compute Ritz vectors.
 *
 *  HOWMNY  Character*1  (INPUT)
 *          Specifies how many Ritz vectors are wanted and the form of Z
 *          the matrix of Ritz vectors. See remark 1 below.
 *          = 'A': compute NEV Ritz vectors;
 *          = 'S': compute some of the Ritz vectors, specified
 *                 by the logical array SELECT.
 *
 *  SELECT  Logical array of dimension NCV.  (INPUT/WORKSPACE)
 *          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
 *          computed. To select the Ritz vector corresponding to a
 *          Ritz value D(j), SELECT(j) must be set to .TRUE..
 *          If HOWMNY = 'A' , SELECT is used as a workspace for
 *          reordering the Ritz values.
 *
 *  D       Double precision  array of dimension NEV.  (OUTPUT)
 *          On exit, D contains the Ritz value approximations to the
 *          eigenvalues of A*z = lambda*B*z. The values are returned
 *          in ascending order. If IPARAM(7) = 3,4,5 then D represents
 *          the Ritz values of OP computed by dsaupd  transformed to
 *          those of the original eigensystem A*z = lambda*B*z. If
 *          IPARAM(7) = 1,2 then the Ritz values of OP are the same
 *          as the those of A*z = lambda*B*z.
 *
 *  Z       Double precision  N by NEV array if HOWMNY = 'A'.  (OUTPUT)
 *          On exit, Z contains the B-orthonormal Ritz vectors of the
 *          eigensystem A*z = lambda*B*z corresponding to the Ritz
 *          value approximations.
 *          If  RVEC = .FALSE. then Z is not referenced.
 *          NOTE: The array Z may be set equal to first NEV columns of the
 *          Arnoldi/Lanczos basis array V computed by DSAUPD .
 *
 *  LDZ     Integer.  (INPUT)
 *          The leading dimension of the array Z.  If Ritz vectors are
 *          desired, then  LDZ .ge.  max( 1, N ).  In any case,  LDZ .ge. 1.
 *
 *  SIGMA   Double precision   (INPUT)
 *          If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
 *          IPARAM(7) = 1 or 2.
 *
 *  **** The remaining arguments MUST be the same as for the   ****
 *  **** call to DSAUPD  that was just completed.               ****
 *
 *  NOTE: The remaining arguments
 *
 *           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
 *           WORKD, WORKL, LWORKL, INFO
 *
 *         must be passed directly to DSEUPD  following the last call
 *         to DSAUPD .  These arguments MUST NOT BE MODIFIED between
 *         the the last call to DSAUPD  and the call to DSEUPD .
 *
 *  Two of these parameters (WORKL, INFO) are also output parameters:
 *
 *  WORKL   Double precision  work array of length LWORKL.  (OUTPUT/WORKSPACE)
 *          WORKL(1:4*ncv) contains information obtained in
 *          dsaupd .  They are not changed by dseupd .
 *          WORKL(4*ncv+1:ncv*ncv+8*ncv) holds the
 *          untransformed Ritz values, the computed error estimates,
 *          and the associated eigenvector matrix of H.
 *
 *          Note: IPNTR(8:10) contains the pointer into WORKL for addresses
 *          of the above information computed by dseupd .
 *          -------------------------------------------------------------
 *          IPNTR(8): pointer to the NCV RITZ values of the original system.
 *          IPNTR(9): pointer to the NCV corresponding error bounds.
 *          IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
 *                     of the tridiagonal matrix T. Only referenced by
 *                     dseupd  if RVEC = .TRUE. See Remarks.
 *          -------------------------------------------------------------
 *
 *  INFO    Integer.  (OUTPUT)
 *          Error flag on output.
 *          =  0: Normal exit.
 *          = -1: N must be positive.
 *          = -2: NEV must be positive.
 *          = -3: NCV must be greater than NEV and less than or equal to N.
 *          = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
 *          = -6: BMAT must be one of 'I' or 'G'.
 *          = -7: Length of private work WORKL array is not sufficient.
 *          = -8: Error return from trid. eigenvalue calculation;
 *                Information error from LAPACK routine dsteqr .
 *          = -9: Starting vector is zero.
 *          = -10: IPARAM(7) must be 1,2,3,4,5.
 *          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
 *          = -12: NEV and WHICH = 'BE' are incompatible.
 *          = -14: DSAUPD  did not find any eigenvalues to sufficient
 *                 accuracy.
 *          = -15: HOWMNY must be one of 'A' or 'S' if RVEC = .true.
 *          = -16: HOWMNY = 'S' not yet implemented
 *          = -17: DSEUPD  got a different count of the number of converged
 *                 Ritz values than DSAUPD  got.  This indicates the user
 *                 probably made an error in passing data from DSAUPD  to
 *                 DSEUPD  or that the data was modified before entering
 *                 DSEUPD .
 * \EndDoc
 *
 * \BeginLib
 *
 * \References:
 *  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
 *     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
 *     pp 357-385.
 *  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
 *     Restarted Arnoldi Iteration", Rice University Technical Report
 *     TR95-13, Department of Computational and Applied Mathematics.
 *  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
 *     1980.
 *  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
 *     Computer Physics Communications, 53 (1989), pp 169-179.
 *  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
 *     Implement the Spectral Transformation", Math. Comp., 48 (1987),
 *     pp 663-673.
 *  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos
 *     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems",
 *     SIAM J. Matr. Anal. Apps.,  January (1993).
 *  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
 *     for Updating the QR decomposition", ACM TOMS, December 1990,
 *     Volume 16 Number 4, pp 369-377.
 *
 * \Remarks
 *  1. The converged Ritz values are always returned in increasing
 *     (algebraic) order.
 *
 *  2. Currently only HOWMNY = 'A' is implemented. It is included at this
 *     stage for the user who wants to incorporate it.
 *
 * \Routines called:
 *     dsesrt   ARPACK routine that sorts an array X, and applies the
 *             corresponding permutation to a matrix A.
 *     dsortr   dsortr   ARPACK sorting routine.
 *     ivout   ARPACK utility routine that prints integers.
 *     dvout    ARPACK utility routine that prints vectors.
 *     dgeqr2   LAPACK routine that computes the QR factorization of
 *             a matrix.
 *     dlacpy   LAPACK matrix copy routine.
 *     dlamch   LAPACK routine that determines machine constants.
 *     dorm2r   LAPACK routine that applies an orthogonal matrix in
 *             factored form.
 *     dsteqr   LAPACK routine that computes eigenvalues and eigenvectors
 *             of a tridiagonal matrix.
 *     dger     Level 2 BLAS rank one update to a matrix.
 *     dcopy    Level 1 BLAS that copies one vector to another .
 *     dnrm2    Level 1 BLAS that computes the norm of a vector.
 *     dscal    Level 1 BLAS that scales a vector.
 *     dswap    Level 1 BLAS that swaps the contents of two vectors.
 * \Authors
 *     Danny Sorensen               Phuong Vu
 *     Richard Lehoucq              CRPC / Rice University
 *     Chao Yang                    Houston, Texas
 *     Dept. of Computational &
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * \Revision history:
 *     12/15/93: Version ' 2.1'
 *
 * \SCCS Information: @(#)
 * FILE: seupd.F   SID: 2.11   DATE OF SID: 04/10/01   RELEASE: 2
 *
 * \EndLib
 */
int dseupd_(bool *rvec, char *howmny, bool *select, double *d, double *z, int *ldz,
            double *sigma, char *bmat, int *n, char *which, int *nev, double *tol,
            double *resid, int *ncv, double *v, int *ldv, int *iparam, int *ipntr,
            double *workd, double *workl, int *lworkl, int *info)
{
    /* System generated locals */
    int i__1;
    double d__1, d__2, d__3;

    /* Builtin functions */

    /* Local variables */
    int j, k, ih, jj, iq, np, iw, ibd, ihb, ihd, ldh, ldq, irz;
    int mode;
    double eps23;
    int ierr;
    double temp;
    int next;
    char type[7];
    int ritz;
    double temp1;
    bool reord;
    int nconv;
    double rnorm;
    double bnorm2;
    int bounds, msglvl, ishift, numcnv;
    int leftptr, rghtptr;


    /* Parameter adjustments */
    --select;
    --workl;

    /* Function Body */
    msglvl = debug_1.mseupd;
    mode = iparam[6];
    nconv = iparam[4];
    *info = 0;

    /* ------------ */
    /* Quick return */
    /* ------------ */

    if (nconv == 0)
    {
        goto L9000;
    }
    ierr = 0;

    if (nconv <= 0)
    {
        ierr = -14;
    }
    if (*n <= 0)
    {
        ierr = -1;
    }
    if (*nev <= 0)
    {
        ierr = -2;
    }
    if (*ncv <= *nev || *ncv > *n)
    {
        ierr = -3;
    }
    if (strcmp(which, "LM") != 0 && strcmp(which, "SM") != 0 &&
        strcmp(which, "LA") != 0 && strcmp(which, "SA") != 0 &&
        strcmp(which, "BE") != 0)
    {
        ierr = -5;
    }
    if (*bmat != 'I' && *bmat != 'G')
    {
        ierr = -6;
    }
    if (*howmny != 'A' && *howmny != 'P' && *howmny != 'S' && *rvec)
    {
        ierr = -15;
    }
    if (*rvec && *howmny == 'S')
    {
        ierr = -16;
    }

    /* Computing 2nd power */
    i__1 = *ncv;
    if (*rvec && *lworkl < i__1 * i__1 + (*ncv << 3))
    {
        ierr = -7;
    }

    if (mode == 1 || mode == 2)
    {
        strcpy(type, "REGULR");
    }
    else if (mode == 3)
    {
        strcpy(type, "SHIFTI");
    }
    else if (mode == 4)
    {
        strcpy(type, "BUCKLE");
    }
    else if (mode == 5)
    {
        strcpy(type, "CAYLEY");
    }
    else
    {
        ierr = -10;
    }
    if (mode == 1 && *bmat == 'G')
    {
        ierr = -11;
    }
    if (*nev == 1 && strcmp(which, "BE") == 0)
    {
        ierr = -12;
    }

    /* ---------- */
    /* Error Exit */
    /* ---------- */

    if (ierr != 0)
    {
        *info = ierr;
        goto L9000;
    }

    /* ----------------------------------------------------- */
    /* Pointer into WORKL for address of H, RITZ, BOUNDS, Q  */
    /* etc... and the remaining workspace.                   */
    /* Also update pointer to be used on output.             */
    /* Memory is laid out as follows:                        */
    /* workl(1:2*ncv) := generated tridiagonal matrix H      */
    /*       The subdiagonal is stored in workl(2:ncv).      */
    /*       The dead spot is workl(1) but upon exiting      */
    /*       dsaupd  stores the B-norm of the last residual  */
    /*       vector in workl(1). We use this !!!             */
    /* workl(2*ncv+1:2*ncv+ncv) := ritz values               */
    /*       The wanted values are in the first NCONV spots. */
    /* workl(3*ncv+1:3*ncv+ncv) := computed Ritz estimates   */
    /*       The wanted values are in the first NCONV spots. */
    /* NOTE: workl(1:4*ncv) is set by dsaupd  and is not     */
    /*       modified by dseupd .                            */
    /* ----------------------------------------------------- */

    /* ----------------------------------------------------- */
    /* The following is used and set by dseupd .             */
    /* workl(4*ncv+1:4*ncv+ncv) := used as workspace during  */
    /*       computation of the eigenvectors of H. Stores    */
    /*       the diagonal of H. Upon EXIT contains the NCV   */
    /*       Ritz values of the original system. The first   */
    /*       NCONV spots have the wanted values. If MODE =   */
    /*       1 or 2 then will equal workl(2*ncv+1:3*ncv).    */
    /* workl(5*ncv+1:5*ncv+ncv) := used as workspace during  */
    /*       computation of the eigenvectors of H. Stores    */
    /*       the subdiagonal of H. Upon EXIT contains the    */
    /*       NCV corresponding Ritz estimates of the         */
    /*       original system. The first NCONV spots have the */
    /*       wanted values. If MODE = 1,2 then will equal    */
    /*       workl(3*ncv+1:4*ncv).                           */
    /* workl(6*ncv+1:6*ncv+ncv*ncv) := orthogonal Q that is  */
    /*       the eigenvector matrix for H as returned by     */
    /*       dsteqr . Not referenced if RVEC = .False.       */
    /*       Ordering follows that of workl(4*ncv+1:5*ncv)   */
    /* workl(6*ncv+ncv*ncv+1:6*ncv+ncv*ncv+2*ncv) :=         */
    /*       Workspace. Needed by dsteqr  and by dseupd .    */
    /* GRAND total of NCV*(NCV+8) locations.                 */
    /* ----------------------------------------------------- */

    ih = ipntr[4];
    ritz = ipntr[5];
    bounds = ipntr[6];
    ldh = *ncv;
    ldq = *ncv;
    ihd = bounds + ldh;
    ihb = ihd + ldh;
    iq = ihb + ldh;
    iw = iq + ldh * *ncv;
    next = iw + (*ncv << 1);
    ipntr[3] = next;
    ipntr[7] = ihd;
    ipntr[8] = ihb;
    ipntr[9] = iq;

    /* -------------------------------------- */
    /* irz points to the Ritz values computed */
    /*     by _seigt before exiting _saup2.   */
    /* ibd points to the Ritz estimates       */
    /*     computed by _seigt before exiting  */
    /*     _saup2.                            */
    /* -------------------------------------- */

    irz = ipntr[10] + *ncv;
    ibd = irz + *ncv;

    /* ------------------------------- */
    /* Set machine dependent constant. */
    /* ------------------------------- */

    eps23 = dlamch_("E");
    eps23 = pow(eps23, d_23);

    /* ------------------------------------- */
    /* RNORM is B-norm of the RESID(1:N).    */
    /* BNORM2 is the 2 norm of B*RESID(1:N). */
    /* Upon exit of dsaupd  WORKD(1:N) has   */
    /* B*RESID(1:N).                         */
    /* ------------------------------------- */

    rnorm = workl[ih];
    if (*bmat == 'I')
    {
        bnorm2 = rnorm;
    }
    else if (*bmat == 'G')
    {
        bnorm2 = dnrm2_(n, workd, &c__1);
    }

#ifndef NO_TRACE
    if (msglvl > 2)
    {
        dvout_(ncv, &workl[irz], &debug_1.ndigit, "_seupd: Ritz values passed in from _SAUPD.");
        dvout_(ncv, &workl[ibd], &debug_1.ndigit, "_seupd: Ritz estimates passed in from _SAUPD.");
    }
#endif

    if (*rvec)
    {
        reord = false;

        /* ------------------------------------------------- */
        /* Use the temporary bounds array to store indices   */
        /* These will be used to mark the select array later */
        /* ------------------------------------------------- */

        i__1 = *ncv;
        for (j = 1; j <= i__1; ++j)
        {
            workl[bounds + j - 1] = (double) j;
            select[j] = false;
        }

        /* ----------------------------------- */
        /* Select the wanted Ritz values.      */
        /* Sort the Ritz values so that the    */
        /* wanted ones appear at the tailing   */
        /* NEV positions of workl(irr) and     */
        /* workl(iri).  Move the corresponding */
        /* error estimates in workl(bound)     */
        /* accordingly.                        */
        /* ----------------------------------- */

        np = *ncv - *nev;
        ishift = 0;
        dsgets_(&ishift, which, nev, &np, &workl[irz], &workl[bounds], &workl[1]);

#ifndef NO_TRACE
        if (msglvl > 2)
        {
            dvout_(ncv, &workl[irz], &debug_1.ndigit, "_seupd: Ritz values after calling _SGETS.");
            dvout_(ncv, &workl[bounds], &debug_1.ndigit, "_seupd: Ritz value indices after calling _SGETS.");
        }
#endif

        /* --------------------------------------------------- */
        /* Record indices of the converged wanted Ritz values  */
        /* Mark the select array for possible reordering       */
        /* --------------------------------------------------- */

        numcnv = 0;
        i__1 = *ncv;
        for (j = 1; j <= i__1; ++j)
        {
            /* Computing MAX */
            d__2 = eps23, d__3 = (d__1 = workl[irz + *ncv - j], abs(d__1));
            temp1 = max(d__2,d__3);
            jj = (int) workl[bounds + *ncv - j];
            if (numcnv < nconv && workl[ibd + jj - 1] <= *tol * temp1)
            {
                select[jj] = true;
                ++numcnv;
                if (jj > nconv)
                {
                    reord = true;
                }
            }
        }

        /* --------------------------------------------------------- */
        /* Check the count (numcnv) of converged Ritz values with    */
        /* the number (nconv) reported by _saupd.  If these two      */
        /* are different then there has probably been an error       */
        /* caused by incorrect passing of the _saupd data.           */
        /* --------------------------------------------------------- */

#ifndef NO_TRACE
        if (msglvl > 2)
        {
            ivout_(&c__1, &numcnv, &debug_1.ndigit, "_seupd: Number of specified eigenvalues");
            ivout_(&c__1, &nconv, &debug_1.ndigit, "_seupd: Number of \"converged\" eigenvalues");
        }
#endif

        if (numcnv != nconv)
        {
            *info = -17;
            goto L9000;
        }

        /* --------------------------------------------------------- */
        /* Call LAPACK routine _steqr to compute the eigenvalues and */
        /* eigenvectors of the final symmetric tridiagonal matrix H. */
        /* Initialize the eigenvector matrix Q to the identity.      */
        /* --------------------------------------------------------- */

        i__1 = *ncv - 1;
        dcopy_(&i__1, &workl[ih + 1], &c__1, &workl[ihb], &c__1);
        dcopy_(ncv, &workl[ih + ldh], &c__1, &workl[ihd], &c__1);

        dsteqr_("I", ncv, &workl[ihd], &workl[ihb], &workl[iq], &ldq, &workl[iw], &ierr);

        if (ierr != 0)
        {
            *info = -8;
            goto L9000;
        }

#ifndef NO_TRACE
        if (msglvl > 1)
        {
            dcopy_(ncv, &workl[iq + *ncv - 1], &ldq, &workl[iw], &c__1);
            dvout_(ncv, &workl[ihd], &debug_1.ndigit, "_seupd: NCV Ritz values of the final H matrix");
            dvout_(ncv, &workl[iw], &debug_1.ndigit, "_seupd: last row of the eigenvector matrix for H");
        }
#endif

        if (reord)
        {
            /* ------------------------------------------- */
            /* Reordered the eigenvalues and eigenvectors  */
            /* computed by _steqr so that the "converged"  */
            /* eigenvalues appear in the first NCONV       */
            /* positions of workl(ihd), and the associated */
            /* eigenvectors appear in the first NCONV      */
            /* columns.                                    */
            /* ------------------------------------------- */

            leftptr = 1;
            rghtptr = *ncv;

            if (*ncv == 1)
            {
                goto L30;
            }

L20:
            if (select[leftptr])
            {
                /* ----------------------------------------- */
                /* Search, from the left, for the first Ritz */
                /* value that has not converged.             */
                /* ----------------------------------------- */

                ++leftptr;
            }
            else if (! select[rghtptr])
            {
                /* -------------------------------------------- */
                /* Search, from the right, the first Ritz value */
                /* that has converged.                          */
                /* -------------------------------------------- */

                --rghtptr;
            }
            else
            {
                /* -------------------------------------------- */
                /* Swap the Ritz value on the left that has not */
                /* converged with the Ritz value on the right   */
                /* that has converged.  Swap the associated     */
                /* eigenvector of the tridiagonal matrix H as   */
                /* well.                                        */
                /* -------------------------------------------- */

                temp = workl[ihd + leftptr - 1];
                workl[ihd + leftptr - 1] = workl[ihd + rghtptr - 1];
                workl[ihd + rghtptr - 1] = temp;
                dcopy_(ncv, &workl[iq + *ncv * (leftptr - 1)], &c__1, &workl[iw], &c__1);
                dcopy_(ncv, &workl[iq + *ncv * (rghtptr - 1)], &c__1, &workl[iq + *ncv * (leftptr - 1)], &c__1);
                dcopy_(ncv, &workl[iw], &c__1, &workl[iq + *ncv * (rghtptr - 1)], &c__1);
                ++leftptr;
                --rghtptr;
            }

            if (leftptr < rghtptr)
            {
                goto L20;
            }
        }

L30:
#ifndef NO_TRACE
        if (msglvl > 2)
        {
            dvout_(ncv, &workl[ihd], &debug_1.ndigit, "_seupd: The eigenvalues of H--reordered");
        }
#endif

        /* -------------------------------------- */
        /* Load the converged Ritz values into D. */
        /* -------------------------------------- */

        dcopy_(&nconv, &workl[ihd], &c__1, d, &c__1);
    }
    else
    {
        /* --------------------------------------------------- */
        /* Ritz vectors not required. Load Ritz values into D. */
        /* --------------------------------------------------- */

        dcopy_(&nconv, &workl[ritz], &c__1, d, &c__1);
        dcopy_(ncv, &workl[ritz], &c__1, &workl[ihd], &c__1);
    }

    /* ---------------------------------------------------------------- */
    /* Transform the Ritz values and possibly vectors and corresponding */
    /* Ritz estimates of OP to those of A*x=lambda*B*x. The Ritz values */
    /* (and corresponding data) are returned in ascending order.        */
    /* ---------------------------------------------------------------- */

    if (strcmp(type, "REGULR") == 0)
    {
        /* ------------------------------------------------------- */
        /* Ascending sort of wanted Ritz values, vectors and error */
        /* bounds. Not necessary if only Ritz values are desired.  */
        /* ------------------------------------------------------- */

        if (*rvec)
        {
            dsesrt_("LA", rvec, &nconv, d, ncv, &workl[iq], &ldq);
        }
        else
        {
            dcopy_(ncv, &workl[bounds], &c__1, &workl[ihb], &c__1);
        }
    }
    else
    {
        /* ----------------------------------------------------------- */
        /* *  Make a copy of all the Ritz values.                      */
        /* *  Transform the Ritz values back to the original system.   */
        /*    For TYPE = 'SHIFTI' the transformation is                */
        /*             lambda = 1/theta + sigma                        */
        /*    For TYPE = 'BUCKLE' the transformation is                */
        /*             lambda = sigma * theta / ( theta - 1 )          */
        /*    For TYPE = 'CAYLEY' the transformation is                */
        /*             lambda = sigma * (theta + 1) / (theta - 1 )     */
        /*    where the theta are the Ritz values returned by dsaupd . */
        /* NOTES:                                                      */
        /* *The Ritz vectors are not affected by the transformation.   */
        /*  They are only reordered.                                   */
        /* ----------------------------------------------------------- */

        dcopy_(ncv, &workl[ihd], &c__1, &workl[iw], &c__1);
        if (strcmp(type, "SHIFTI") == 0)
        {
            i__1 = *ncv;
            for (k = 0; k < i__1; ++k)
            {
                workl[ihd + k] = 1.0 / workl[ihd + k] + *sigma;
            }
        }
        else if (strcmp(type, "BUCKLE") == 0)
        {
            i__1 = *ncv;
            for (k = 0; k < i__1; ++k)
            {
                workl[ihd + k] = *sigma * workl[ihd + k] / (workl[ihd + k] - 1.0);
            }
        }
        else if (strcmp(type, "CAYLEY") == 0)
        {
            i__1 = *ncv;
            for (k = 0; k < i__1; ++k)
            {
                workl[ihd + k] = *sigma * (workl[ihd + k] + 1.0) / (workl[ihd + k] - 1.0);
            }
        }

        /* ----------------------------------------------------------- */
        /* *  Store the wanted NCONV lambda values into D.             */
        /* *  Sort the NCONV wanted lambda in WORKL(IHD:IHD+NCONV-1)   */
        /*    into ascending order and apply sort to the NCONV theta   */
        /*    values in the transformed system. We will need this to   */
        /*    compute Ritz estimates in the original system.           */
        /* *  Finally sort the lambda`s into ascending order and apply */
        /*    to Ritz vectors if wanted. Else just sort lambda`s into  */
        /*    ascending order.                                         */
        /* NOTES:                                                      */
        /* *workl(iw:iw+ncv-1) contain the theta ordered so that they  */
        /*  match the ordering of the lambda. We`ll use them again for */
        /*  Ritz vector purification.                                  */
        /* ----------------------------------------------------------- */

        dcopy_(&nconv, &workl[ihd], &c__1, d, &c__1);
        dsortr_("LA", &c_true, &nconv, &workl[ihd], &workl[iw]);
        if (*rvec)
        {
            dsesrt_("LA", rvec, &nconv, d, ncv, &workl[iq], &ldq);
        }
        else
        {
            dcopy_(ncv, &workl[bounds], &c__1, &workl[ihb], &c__1);
            d__1 = bnorm2 / rnorm;
            dscal_(ncv, &d__1, &workl[ihb], &c__1);
            dsortr_("LA", &c_true, &nconv, d, &workl[ihb]);
        }
    }

    /* ---------------------------------------------- */
    /* Compute the Ritz vectors. Transform the wanted */
    /* eigenvectors of the symmetric tridiagonal H by */
    /* the Lanczos basis matrix V.                    */
    /* ---------------------------------------------- */

    if (*rvec && *howmny == 'A')
    {
        /* -------------------------------------------------------- */
        /* Compute the QR factorization of the matrix representing  */
        /* the wanted invariant subspace located in the first NCONV */
        /* columns of workl(iq,ldq).                                */
        /* -------------------------------------------------------- */

        dgeqr2_(ncv, &nconv, &workl[iq], &ldq, &workl[iw + *ncv], &workl[ihb],&ierr);

        /* ------------------------------------------------------ */
        /* * Postmultiply V by Q.                                 */
        /* * Copy the first NCONV columns of VQ into Z.           */
        /* The N by NCONV matrix Z is now a matrix representation */
        /* of the approximate invariant subspace associated with  */
        /* the Ritz values in workl(ihd).                         */
        /* ------------------------------------------------------ */

        dorm2r_("R", "N", n, ncv, &nconv, &workl[iq], &ldq, &workl[iw + *ncv], v, ldv, &workd[*n], &ierr);
        dlacpy_("A", n, &nconv, v, ldv, z, ldz);

        /* --------------------------------------------------- */
        /* In order to compute the Ritz estimates for the Ritz */
        /* values in both systems, need the last row of the    */
        /* eigenvector matrix. Remember, it`s in factored form */
        /* --------------------------------------------------- */

        i__1 = *ncv - 1;
        for (j = 1; j <= i__1; ++j)
        {
            workl[ihb + j - 1] = 0.0;
        }
        workl[ihb + *ncv - 1] = 1.0;
        dorm2r_("L", "T", ncv, &c__1, &nconv, &workl[iq], &ldq, &workl[iw + *ncv], &workl[ihb], ncv, &temp, &ierr);

        /* --------------------------------------------------- */
        /* Make a copy of the last row into                    */
        /* workl(iw+ncv:iw+2*ncv), as it is needed again in    */
        /* the Ritz vector purification step below             */
        /* --------------------------------------------------- */

        for (j = 1; j <= nconv; ++j)
        {
            workl[iw + *ncv + j - 1] = workl[ihb + j - 1];
        }
    }
    else if (*rvec && *howmny == 'S')
    {
        /*     Not yet implemented. See remark 2 above. */
    }

    if (strcmp(type, "REGULR") == 0 && *rvec)
    {
        i__1 = *ncv;
        for (j = 1; j <= i__1; ++j)
        {
            d__1 = workl[ihb + j - 1];
            workl[ihb + j - 1] = rnorm * abs(d__1);
        }
    }
    else if (strcmp(type, "REGULR") != 0 && *rvec)
    {
        /* ----------------------------------------------- */
        /* *  Determine Ritz estimates of the theta.       */
        /*    If RVEC = .true. then compute Ritz estimates */
        /*               of the theta.                     */
        /*    If RVEC = .false. then copy Ritz estimates   */
        /*              as computed by dsaupd .            */
        /* *  Determine Ritz estimates of the lambda.      */
        /* ----------------------------------------------- */

        dscal_(ncv, &bnorm2, &workl[ihb], &c__1);
        if (strcmp(type, "SHIFTI") == 0)
        {
            i__1 = *ncv;
            for (k = 0; k < i__1; ++k)
            {
                /* Computing 2nd power */
                d__1 = workl[ihb + k];
                d__2 = workl[iw + k];
                workl[ihb + k] = abs(d__1) / (d__2 * d__2);
            }

        }
        else if (strcmp(type, "BUCKLE") == 0)
        {
            i__1 = *ncv;
            for (k = 0; k < i__1; ++k)
            {
                /* Computing 2nd power */
                d__1 = workl[ihb + k];
                d__2 = workl[iw + k] - 1.0;
                workl[ihb + k] = *sigma * abs(d__1) / (d__2 * d__2);
            }
        }
        else if (strcmp(type, "CAYLEY") == 0)
        {
            i__1 = *ncv;
            for (k = 0; k < i__1; ++k)
            {
                d__1 = workl[ihb + k] / workl[iw + k] * (workl[iw + k] - 1.0);
                workl[ihb + k] = abs(d__1);
            }
        }
    }

#ifndef NO_TRACE
    if (msglvl > 1 && strcmp(type, "REGULR") != 0)
    {
        dvout_(&nconv, d, &debug_1.ndigit, "_seupd: Untransformed converged Ritz values");
        dvout_(&nconv, &workl[ihb], &debug_1.ndigit, "_seupd: Ritz estimates of the untransformed Ritz values");
    }
    else if (msglvl > 1)
    {
        dvout_(&nconv, d, &debug_1.ndigit, "_seupd: Converged Ritz values");
        dvout_(&nconv, &workl[ihb], &debug_1.ndigit, "_seupd: Associated Ritz estimates");
    }
#endif

    /* ----------------------------------------------- */
    /* Ritz vector purification step. Formally perform */
    /* one of inverse subspace iteration. Only used    */
    /* for MODE = 3,4,5. See reference 7               */
    /* ----------------------------------------------- */

    if (*rvec && (strcmp(type, "SHIFTI") == 0 || strcmp(type, "CAYLEY") == 0))
    {
        i__1 = nconv - 1;
        for (k = 0; k <= i__1; ++k)
        {
            workl[iw + k] = workl[iw + *ncv + k] / workl[iw + k];
        }
    }
    else if (*rvec && strcmp(type, "BUCKLE") == 0)
    {
        i__1 = nconv - 1;
        for (k = 0; k <= i__1; ++k)
        {
            workl[iw + k] = workl[iw + *ncv + k] / (workl[iw + k] - 1.0);
        }
    }

    if (*rvec && strcmp(type, "REGULR") != 0)
    {
        dger_(n, &nconv, &d_one, resid, &c__1, &workl[iw], &c__1, z, ldz);
    }

L9000:

    return 0;

    /* ------------- */
    /* End of dseupd */
    /* ------------- */

} /* dseupd_ */

