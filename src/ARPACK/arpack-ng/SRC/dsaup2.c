/* arpack-ng\SRC\dsaup2.f -- translated by f2c (version 20100827). */

#include <math.h>
#include "arpack_internal.h"

/**
 * \BeginDoc
 *
 * \Name: dsaup2
 *
 * \Description:
 *  Intermediate level interface called by dsaupd.
 *
 * \Usage:
 *  call dsaup2
 *     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
 *       ISHIFT, MXITER, V, LDV, H, LDH, RITZ, BOUNDS, Q, LDQ, WORKL,
 *       IPNTR, WORKD, INFO )
 *
 * \Arguments
 *
 *  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in dsaupd.
 *  MODE, ISHIFT, MXITER: see the definition of IPARAM in dsaupd.
 *
 *  NP      Integer.  (INPUT/OUTPUT)
 *          Contains the number of implicit shifts to apply during
 *          each Arnoldi/Lanczos iteration.
 *          If ISHIFT=1, NP is adjusted dynamically at each iteration
 *          to accelerate convergence and prevent stagnation.
 *          This is also roughly equal to the number of matrix-vector
 *          products (involving the operator OP) per Arnoldi iteration.
 *          The logic for adjusting is contained within the current
 *          subroutine.
 *          If ISHIFT=0, NP is the number of shifts the user needs
 *          to provide via reverse communication. 0 < NP < NCV-NEV.
 *          NP may be less than NCV-NEV since a leading block of the current
 *          upper Tridiagonal matrix has split off and contains "unwanted"
 *          Ritz values.
 *          Upon termination of the IRA iteration, NP contains the number
 *          of "converged" wanted Ritz values.
 *
 *  IUPD    Integer.  (INPUT)
 *          IUPD .EQ. 0: use explicit restart instead implicit update.
 *          IUPD .NE. 0: use implicit update.
 *
 *  V       Double precision N by (NEV+NP) array.  (INPUT/OUTPUT)
 *          The Lanczos basis vectors.
 *
 *  LDV     Integer.  (INPUT)
 *          Leading dimension of V exactly as declared in the calling
 *          program.
 *
 *  H       Double precision (NEV+NP) by 2 array.  (OUTPUT)
 *          H is used to store the generated symmetric tridiagonal matrix
 *          The subdiagonal is stored in the first column of H starting
 *          at H(2,1).  The main diagonal is stored in the arscnd column
 *          of H starting at H(1,2). If dsaup2 converges store the
 *          B-norm of the final residual vector in H(1,1).
 *
 *  LDH     Integer.  (INPUT)
 *          Leading dimension of H exactly as declared in the calling
 *          program.
 *
 *  RITZ    Double precision array of length NEV+NP.  (OUTPUT)
 *          RITZ(1:NEV) contains the computed Ritz values of OP.
 *
 *  BOUNDS  Double precision array of length NEV+NP.  (OUTPUT)
 *          BOUNDS(1:NEV) contain the error bounds corresponding to RITZ.
 *
 *  Q       Double precision (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
 *          Private (replicated) work array used to accumulate the
 *          rotation in the shift application step.
 *
 *  LDQ     Integer.  (INPUT)
 *          Leading dimension of Q exactly as declared in the calling
 *          program.
 *
 *  WORKL   Double precision array of length at least 3*(NEV+NP).  (INPUT/WORKSPACE)
 *          Private (replicated) array on each PE or array allocated on
 *          the front end.  It is used in the computation of the
 *          tridiagonal eigenvalue problem, the calculation and
 *          application of the shifts and convergence checking.
 *          If ISHIFT .EQ. O and IDO .EQ. 3, the first NP locations
 *          of WORKL are used in reverse communication to hold the user
 *          supplied shifts.
 *
 *  IPNTR   Integer array of length 3.  (OUTPUT)
 *          Pointer to mark the starting locations in the WORKD for
 *          vectors used by the Lanczos iteration.
 *          -------------------------------------------------------------
 *          IPNTR(1): pointer to the current operand vector X.
 *          IPNTR(2): pointer to the current result vector Y.
 *          IPNTR(3): pointer to the vector B * X when used in one of
 *                    the spectral transformation modes.  X is the current
 *                    operand.
 *          -------------------------------------------------------------
 *
 *  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
 *          Distributed array to be used in the basic Lanczos iteration
 *          for reverse communication.  The user should not use WORKD
 *          as temporary workspace during the iteration !!!!!!!!!!
 *          See Data Distribution Note in dsaupd.
 *
 *  INFO    Integer.  (INPUT/OUTPUT)
 *          If INFO .EQ. 0, a randomly initial residual vector is used.
 *          If INFO .NE. 0, RESID contains the initial residual vector,
 *                          possibly from a previous run.
 *          Error flag on output.
 *          =     0: Normal return.
 *          =     1: All possible eigenvalues of OP has been found.
 *                   NP returns the size of the invariant subspace
 *                   spanning the operator OP.
 *          =     2: No shifts could be applied.
 *          =    -8: Error return from trid. eigenvalue calculation;
 *                   This should never happen.
 *          =    -9: Starting vector is zero.
 *          = -9999: Could not build an Lanczos factorization.
 *                   Size that was built in returned in NP.
 *
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
 * \Routines called:
 *     dgetv0  ARPACK initial vector generation routine.
 *     dsaitr  ARPACK Lanczos factorization routine.
 *     dsapps  ARPACK application of implicit shifts routine.
 *     dsconv  ARPACK convergence of Ritz values routine.
 *     dseigt  ARPACK compute Ritz values and error bounds routine.
 *     dsgets  ARPACK reorder Ritz values and error bounds routine.
 *     dsortr  ARPACK sorting routine.
 *     ivout   ARPACK utility routine that prints integers.
 *     arscnd  ARPACK utility routine for timing.
 *     dvout   ARPACK utility routine that prints vectors.
 *     dlamch  LAPACK routine that determines machine constants.
 *     dcopy   Level 1 BLAS that copies one vector to another.
 *     ddot    Level 1 BLAS that computes the scalar product of two vectors.
 *     dnrm2   Level 1 BLAS that computes the norm of a vector.
 *     dscal   Level 1 BLAS that scales a vector.
 *     dswap   Level 1 BLAS that swaps two vectors.
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
 *     12/15/93: Version ' 2.4'
 *     xx/xx/95: Version ' 2.4'.  (R.B. Lehoucq)
 *
 * \SCCS Information: @(#)
 * FILE: saup2.F   SID: 2.7   DATE OF SID: 5/19/98   RELEASE: 2
 *
 * \EndLib
 */
int dsaup2_(int *ido, char *bmat, int *n, char *which, int *nev, int *np,
            double *tol, double *resid, int *mode, int *iupd, int *ishift, int *mxiter,
            double *v, int *ldv, double *h, int *ldh, double *ritz, double *bounds,
            double *q, int *ldq, double *workl, int *ipntr, double *workd, int *info)
{
    /* System generated locals */
    int i__1, i__2, i__3;
    double d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    int j;
    static float t0, t1, t2, t3;
    int kp[3];
    static int np0, nev0;
    static double eps23;
    int ierr;
    static int iter;
    double temp;
    int nevd2;
    static bool getv0;
    int nevm2;
    static bool cnorm;
    static int nconv;
    static bool initv;
    static double rnorm;
    int nevbef;
    static bool update;
    char wprime[3];
    static bool ushift;
    static int kplusp, msglvl;
    int nptemp;

    /* Function Body */
    if (*ido == 0)
    {
        /* ----------------------------- */
        /* Initialize timing statistics  */
        /* & message level for debugging */
        /* ----------------------------- */

#ifndef NO_TIMER
        arscnd_(&t0);
#endif

        msglvl = debug_1.msaup2;

        /* ------------------------------- */
        /* Set machine dependent constant. */
        /* ------------------------------- */

        eps23 = dlamch_("E");
        eps23 = pow(eps23, d_23);

        /* ----------------------------------- */
        /* nev0 and np0 are integer variables  */
        /* hold the initial values of NEV & NP */
        /* ----------------------------------- */

        nev0 = *nev;
        np0 = *np;

        /* ----------------------------------- */
        /* kplusp is the bound on the largest  */
        /*        Lanczos factorization built. */
        /* nconv is the current number of      */
        /*        "converged" eigenvlues.      */
        /* iter is the counter on the current  */
        /*      iteration step.                */
        /* ----------------------------------- */

        kplusp = nev0 + np0;
        nconv = 0;
        iter = 0;

        /* ------------------------------------------- */
        /* Set flags for computing the first NEV steps */
        /* of the Lanczos factorization.               */
        /* ------------------------------------------- */

        getv0 = true;
        update = false;
        ushift = false;
        cnorm = false;

        if (*info != 0)
        {
            /* ------------------------------------------ */
            /* User provides the initial residual vector. */
            /* ------------------------------------------ */

            initv = true;
            *info = 0;
        }
        else
        {
            initv = false;
        }
    }

    /* ------------------------------------------- */
    /* Get a possibly random starting vector and   */
    /* force it into the range of the operator OP. */
    /* ------------------------------------------- */


    if (getv0)
    {
        dgetv0_(ido, bmat, &c__1, &initv, n, &c__1, v, ldv, resid, &rnorm, ipntr, workd, info);

        if (*ido != 99)
        {
            goto L9000;
        }

        if (rnorm == 0.0)
        {
            /* --------------------------------------- */
            /* The initial vector is zero. Error exit. */
            /* --------------------------------------- */

            *info = -9;
            goto L1200;
        }
        getv0 = false;
        *ido = 0;
    }

    /* ---------------------------------------------------------- */
    /* Back from reverse communication: continue with update step */
    /* ---------------------------------------------------------- */

    if (update)
    {
        goto L20;
    }

    /* ----------------------------------------- */
    /* Back from computing user specified shifts */
    /* ----------------------------------------- */

    if (ushift)
    {
        goto L50;
    }

    /* ----------------------------------- */
    /* Back from computing residual norm   */
    /* at the end of the current iteration */
    /* ----------------------------------- */

    if (cnorm)
    {
        goto L100;
    }

    /* -------------------------------------------------------- */
    /* Compute the first NEV steps of the Lanczos factorization */
    /* -------------------------------------------------------- */

    dsaitr_(ido, bmat, n, &c__0, &nev0, mode, resid, &rnorm, v,ldv, h, ldh, ipntr, workd, info);

    /* ------------------------------------------------- */
    /* ido .ne. 99 implies use of reverse communication  */
    /* to compute operations involving OP and possibly B */
    /* ------------------------------------------------- */

    if (*ido != 99)
    {
        goto L9000;
    }

    if (*info > 0)
    {
        /* --------------------------------------------------- */
        /* dsaitr was unable to build an Lanczos factorization */
        /* of length NEV0. INFO is returned with the size of   */
        /* the factorization built. Exit main loop.            */
        /* --------------------------------------------------- */

        *np = *info;
        *mxiter = iter;
        *info = -9999;
        goto L1200;
    }

    /* ------------------------------------------------------------ */
    /*                                                              */
    /*           M A I N  LANCZOS  I T E R A T I O N  L O O P       */
    /*           Each iteration implicitly restarts the Lanczos     */
    /*           factorization in place.                            */
    /*                                                              */
    /* ------------------------------------------------------------ */

L1000:

    ++iter;

#ifndef NO_TRACE
    if (msglvl > 0)
    {
        ivout_(&c__1, &iter, &debug_1.ndigit, "_saup2: **** Start of major iteration number ****");
    }

    if (msglvl > 1)
    {
        ivout_(&c__1, nev, &debug_1.ndigit, "_saup2: The length of the current Lanczos factorization");
        ivout_(&c__1, np, &debug_1.ndigit, "_saup2: Extend the Lanczos factorization by");
    }
#endif

    /* ---------------------------------------------------------- */
    /* Compute NP additional steps of the Lanczos factorization. */
    /* ---------------------------------------------------------- */

    *ido = 0;
L20:
    update = true;

    dsaitr_(ido, bmat, n, nev, np, mode, resid, &rnorm, v, ldv,h, ldh, ipntr, workd, info);

    /* ------------------------------------------------- */
    /* ido .ne. 99 implies use of reverse communication  */
    /* to compute operations involving OP and possibly B */
    /* ------------------------------------------------- */

    if (*ido != 99)
    {
        goto L9000;
    }

    if (*info > 0)
    {
        /* --------------------------------------------------- */
        /* dsaitr was unable to build an Lanczos factorization */
        /* of length NEV0+NP0. INFO is returned with the size  */
        /* of the factorization built. Exit main loop.         */
        /* --------------------------------------------------- */

        *np = *info;
        *mxiter = iter;
        *info = -9999;
        goto L1200;
    }
    update = false;

#ifndef NO_TRACE
    if (msglvl > 1)
    {
        dvout_(&c__1, &rnorm, &debug_1.ndigit, "_saup2: Current B-norm of residual for factorization");
    }
#endif

    /* ------------------------------------------------------ */
    /* Compute the eigenvalues and corresponding error bounds */
    /* of the current symmetric tridiagonal matrix.           */
    /* ------------------------------------------------------ */

    dseigt_(&rnorm, &kplusp, h, ldh, ritz, bounds, workl, &ierr);

    if (ierr != 0)
    {
        *info = -8;
        goto L1200;
    }

    /* -------------------------------------------------- */
    /* Make a copy of eigenvalues and corresponding error */
    /* bounds obtained from _seigt.                       */
    /* -------------------------------------------------- */

    dcopy_(&kplusp, ritz, &c__1, &workl[kplusp], &c__1);
    dcopy_(&kplusp, bounds, &c__1, &workl[kplusp << 1], &c__1);

    /* ------------------------------------------------- */
    /* Select the wanted Ritz values and their bounds    */
    /* to be used in the convergence test.               */
    /* The selection is based on the requested number of */
    /* eigenvalues instead of the current NEV and NP to  */
    /* prevent possible misconvergence.                  */
    /* * Wanted Ritz values := RITZ(NP+1:NEV+NP)         */
    /* * Shifts := RITZ(1:NP) := WORKL(1:NP)             */
    /* ------------------------------------------------- */

    *nev = nev0;
    *np = np0;
    dsgets_(ishift, which, nev, np, ritz, bounds, workl);

    /* ----------------- */
    /* Convergence test. */
    /* ----------------- */

    dcopy_(nev, &bounds[*np], &c__1, &workl[*np], &c__1);
    dsconv_(nev, &ritz[*np], &workl[*np], tol, &nconv);

#ifndef NO_TRACE
    if (msglvl > 2)
    {
        kp[0] = *nev;
        kp[1] = *np;
        kp[2] = nconv;
        ivout_(&c__3, kp, &debug_1.ndigit, "_saup2: NEV, NP, NCONV are");
        dvout_(&kplusp, ritz, &debug_1.ndigit, "_saup2: The eigenvalues of H");
        dvout_(&kplusp, bounds, &debug_1.ndigit, "_saup2: Ritz estimates of the current NCV Ritz values");
    }
#endif

    /* ------------------------------------------------------- */
    /* Count the number of unwanted Ritz values that have zero */
    /* Ritz estimates. If any Ritz estimates are equal to zero */
    /* then a leading block of H of order equal to at least    */
    /* the number of Ritz values with zero Ritz estimates has  */
    /* split off. None of these Ritz values may be removed by  */
    /* shifting. Decrease NP the number of shifts to apply. If */
    /* no shifts may be applied, then prepare to exit          */
    /* ------------------------------------------------------- */

    nptemp = *np;
    for (j = 0; j < nptemp; ++j)
    {
        if (bounds[j] == 0.0)
        {
            --(*np);
            ++(*nev);
        }
    }

    if (nconv >= nev0 || iter > *mxiter || *np == 0)
    {
        /* ---------------------------------------------- */
        /* Prepare to exit. Put the converged Ritz values */
        /* and corresponding bounds in RITZ(1:NCONV) and  */
        /* BOUNDS(1:NCONV) respectively. Then sort. Be    */
        /* careful when NCONV > NP since we don't want to */
        /* swap overlapping locations.                    */
        /* ---------------------------------------------- */

        if (strcmp(which, "BE") == 0)
        {
            /* --------------------------------------------------- */
            /* Both ends of the spectrum are requested.            */
            /* Sort the eigenvalues into algebraically decreasing  */
            /* order first then swap low end of the spectrum next  */
            /* to high end in appropriate locations.               */
            /* NOTE: when np < floor(nev/2) be careful not to swap */
            /* overlapping locations.                              */
            /* --------------------------------------------------- */

            strcpy(wprime, "SA");
            dsortr_(wprime, &c_true, &kplusp, ritz, bounds);
            nevd2 = nev0 / 2;
            nevm2 = nev0 - nevd2;
            if (*nev > 1)
            {
                *np = kplusp - nev0;
                i__1 = min(nevd2,*np);
                /* Computing MAX */
                i__2 = kplusp - nevd2, i__3 = kplusp - *np;
                dswap_(&i__1, &ritz[nevm2], &c__1, &ritz[max(i__2,i__3)], &c__1);
                dswap_(&i__1, &bounds[nevm2], &c__1, &bounds[max(i__2,i__3)], &c__1);
            }
        }
        else
        {
            /* ------------------------------------------------ */
            /* LM, SM, LA, SA case.                             */
            /* Sort the eigenvalues of H into the an order that */
            /* is opposite to WHICH, and apply the resulting    */
            /* order to BOUNDS.  The eigenvalues are sorted so  */
            /* that the wanted part are always within the first */
            /* NEV locations.                                   */
            /* ------------------------------------------------ */

            if (strcmp(which, "LM") == 0)
            {
                strcpy(wprime, "SM");
            }
            if (strcmp(which, "SM") == 0)
            {
                strcpy(wprime, "LM");
            }
            if (strcmp(which, "LA") == 0)
            {
                strcpy(wprime, "SA");
            }
            if (strcmp(which, "SA") == 0)
            {
                strcpy(wprime, "LA");
            }

            dsortr_(wprime, &c_true, &kplusp, ritz, bounds);
        }

        /* ------------------------------------------------ */
        /* Scale the Ritz estimate of each Ritz value       */
        /* by 1 / max(eps23,magnitude of the Ritz value).   */
        /* ------------------------------------------------ */

        for (j = 0; j < nev0; ++j)
        {
            /* Computing MAX */
            d__2 = eps23, d__3 = (d__1 = ritz[j], abs(d__1));
            temp = max(d__2,d__3);
            bounds[j] /= temp;
        }

        /* -------------------------------------------------- */
        /* Sort the Ritz values according to the scaled Ritz  */
        /* estimates.  This will push all the converged ones  */
        /* towards the front of ritzr, ritzi, bounds          */
        /* (in the case when NCONV < NEV.)                    */
        /* -------------------------------------------------- */

        strcpy(wprime, "LA");
        dsortr_(wprime, &c_true, &nev0, bounds, ritz);

        /* -------------------------------------------- */
        /* Scale the Ritz estimate back to its original */
        /* value.                                       */
        /* -------------------------------------------- */

        for (j = 0; j < nev0; ++j)
        {
            /* Computing MAX */
            d__2 = eps23, d__3 = (d__1 = ritz[j], abs(d__1));
            temp = max(d__2,d__3);
            bounds[j] *= temp;
        }

        /* ------------------------------------------------ */
        /* Sort the "converged" Ritz values again so that   */
        /* the "threshold" values and their associated Ritz */
        /* estimates appear at the appropriate position in  */
        /* ritz and bound.                                  */
        /* ------------------------------------------------ */

        if (strcmp(which, "BE") == 0)
        {
            /* ---------------------------------------------- */
            /* Sort the "converged" Ritz values in increasing */
            /* order.  The "threshold" values are in the      */
            /* middle.                                        */
            /* ---------------------------------------------- */

            strcpy(wprime, "LA");
            dsortr_(wprime, &c_true, &nconv, ritz, bounds);
        }
        else
        {
            /* -------------------------------------------- */
            /* In LM, SM, LA, SA case, sort the "converged" */
            /* Ritz values according to WHICH so that the   */
            /* "threshold" value appears at the front of    */
            /* ritz.                                        */
            /* -------------------------------------------- */
            dsortr_(which, &c_true, &nconv, ritz, bounds);
        }

        /* ---------------------------------------- */
        /*  Use h( 1,1 ) as storage to communicate  */
        /*  rnorm to _seupd if needed               */
        /* ---------------------------------------- */

        h[0] = rnorm;

#ifndef NO_TRACE
        if (msglvl > 1)
        {
            dvout_(&kplusp, ritz, &debug_1.ndigit, "_saup2: Sorted Ritz values.");
            dvout_(&kplusp, bounds, &debug_1.ndigit, "_saup2: Sorted ritz estimates.");
        }
#endif

        /* ---------------------------------- */
        /* Max iterations have been exceeded. */
        /* ---------------------------------- */

        if (iter > *mxiter && nconv < *nev)
        {
            *info = 1;
        }

        /* ------------------- */
        /* No shifts to apply. */
        /* ------------------- */

        if (*np == 0 && nconv < nev0)
        {
            *info = 2;
        }

        *np = nconv;
        goto L1100;
    }
    else if (nconv < *nev && *ishift == 1)
    {
        /* ------------------------------------------------- */
        /* Do not have all the requested eigenvalues yet.    */
        /* To prevent possible stagnation, adjust the number */
        /* of Ritz values and the shifts.                    */
        /* ------------------------------------------------- */

        nevbef = *nev;
        /* Computing MIN */
        i__1 = nconv, i__2 = *np / 2;
        *nev += min(i__1,i__2);
        if (*nev == 1 && kplusp >= 6)
        {
            *nev = kplusp / 2;
        }
        else if (*nev == 1 && kplusp > 2)
        {
            *nev = 2;
        }
        *np = kplusp - *nev;

        /* ------------------------------------- */
        /* If the size of NEV was just increased */
        /* resort the eigenvalues.               */
        /* ------------------------------------- */

        if (nevbef < *nev)
        {
            dsgets_(ishift, which, nev, np, ritz, bounds, workl);
        }
    }

#ifndef NO_TRACE
    if (msglvl > 0)
    {
        ivout_(&c__1, &nconv, &debug_1.ndigit, "_saup2: no. of \"converged\" Ritz values at this iter.");
        if (msglvl > 1)
        {
            kp[0] = *nev;
            kp[1] = *np;
            ivout_(&c__2, kp, &debug_1.ndigit, "_saup2: NEV and NP are");
            dvout_(nev, &ritz[*np], &debug_1.ndigit, "_saup2: \"wanted\" Ritz values.");
            dvout_(nev, &bounds[*np], &debug_1.ndigit, "_saup2: Ritz estimates of the \"wanted\" values ");
        }
    }
#endif

    if (*ishift == 0)
    {
        /* --------------------------------------------------- */
        /* User specified shifts: reverse communication to     */
        /* compute the shifts. They are returned in the first  */
        /* NP locations of WORKL.                              */
        /* --------------------------------------------------- */

        ushift = true;
        *ido = 3;
        goto L9000;
    }

L50:

    /* ---------------------------------- */
    /* Back from reverse communication;   */
    /* User specified shifts are returned */
    /* in WORKL(1:*NP)                    */
    /* ---------------------------------- */

    ushift = false;

    /* ------------------------------------------------------- */
    /* Move the NP shifts to the first NP locations of RITZ to */
    /* free up WORKL.  This is for the non-exact shift case;   */
    /* in the exact shift case, dsgets already handles this.   */
    /* ------------------------------------------------------- */

    if (*ishift == 0)
    {
        dcopy_(np, workl, &c__1, ritz, &c__1);
    }

#ifndef NO_TRACE
    if (msglvl > 2)
    {
        ivout_(&c__1, np, &debug_1.ndigit, "_saup2: The number of shifts to apply ");
        dvout_(np, workl, &debug_1.ndigit, "_saup2: shifts selected");
        if (*ishift == 1)
        {
            dvout_(np, bounds, &debug_1.ndigit, "_saup2: corresponding Ritz estimates");
        }
    }
#endif

    /* ------------------------------------------------------- */
    /* Apply the NP0 implicit shifts by QR bulge chasing.      */
    /* Each shift is applied to the entire tridiagonal matrix. */
    /* The first 2*N locations of WORKD are used as workspace. */
    /* After dsapps is done, we have a Lanczos                 */
    /* factorization of length NEV.                            */
    /* ------------------------------------------------------- */

    dsapps_(n, nev, np, ritz, v, ldv, h, ldh, resid, q, ldq, workd);

    /* ------------------------------------------- */
    /* Compute the B-norm of the updated residual. */
    /* Keep B*RESID in WORKD(1:N) to be used in    */
    /* the first step of the next call to dsaitr.  */
    /* ------------------------------------------- */

    cnorm = true;
#ifndef NO_TIMER
    arscnd_(&t2);
#endif

    if (*bmat == 'G')
    {
        ++timing_1.nbx;
        dcopy_(n, resid, &c__1, &workd[*n], &c__1);
        /* TODO: ipntr index */
        ipntr[0] = *n + 1;
        ipntr[1] = 1;
        *ido = 2;

        /* -------------------------------- */
        /* Exit in order to compute B*RESID */
        /* -------------------------------- */

        goto L9000;
    }
    else if (*bmat == 'I')
    {
        dcopy_(n, resid, &c__1, workd, &c__1);
    }

L100:

    /* -------------------------------- */
    /* Back from reverse communication; */
    /* WORKD(1:N) := B*RESID            */
    /* -------------------------------- */

    if (*bmat == 'G')
    {
#ifndef NO_TIMER
        arscnd_(&t3);
        timing_1.tmvbx += t3 - t2;
#endif
        rnorm = ddot_(n, resid, &c__1, workd, &c__1);
        rnorm = sqrt((abs(rnorm)));
    }
    else if (*bmat == 'I')
    {
        rnorm = dnrm2_(n, resid, &c__1);
    }
    cnorm = false;


#ifndef NO_TRACE
    if (msglvl > 2)
    {
        dvout_(&c__1, &rnorm, &debug_1.ndigit, "_saup2: B-norm of residual for NEV factorization");
        dvout_(nev, &h[*ldh], &debug_1.ndigit,"_saup2: main diagonal of compressed H matrix");
        i__1 = *nev - 1;
        dvout_(&i__1, &h[1], &debug_1.ndigit, "_saup2: subdiagonal of compressed H matrix");
    }
#endif

    goto L1000;

    /* ------------------------------------------------------------- */
    /*                                                               */
    /*  E N D     O F     M A I N     I T E R A T I O N     L O O P  */
    /*                                                               */
    /* ------------------------------------------------------------- */

L1100:

    *mxiter = iter;
    *nev = nconv;

L1200:
    *ido = 99;

    /* ---------- */
    /* Error exit */
    /* ---------- */

#ifndef NO_TIMER
    arscnd_(&t1);
    timing_1.tsaup2 = t1 - t0;
#endif

L9000:
    return 0;

    /* ------------- */
    /* End of dsaup2 */
    /* ------------- */

} /* dsaup2_ */

