/* D:\Projekte\ARPACK\arpack-ng\SRC\cnaup2.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \BeginDoc
 *
 * \Name: cnaup2
 *
 * \Description:
 *  Intermediate level interface called by cnaupd.
 *
 * \Usage:
 *  call cnaup2
 *     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
 *       ISHIFT, MXITER, V, LDV, H, LDH, RITZ, BOUNDS,
 *       Q, LDQ, WORKL, IPNTR, WORKD, RWORK, INFO )
 *
 * \Arguments
 *
 *  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in cnaupd.
 *  MODE, ISHIFT, MXITER: see the definition of IPARAM in cnaupd.
 *
 *  NP      Integer.  (INPUT/OUTPUT)
 *          Contains the number of implicit shifts to apply during
 *          each Arnoldi iteration.
 *          If ISHIFT=1, NP is adjusted dynamically at each iteration
 *          to accelerate convergence and prevent stagnation.
 *          This is also roughly equal to the number of matrix-vector
 *          products (involving the operator OP) per Arnoldi iteration.
 *          The logic for adjusting is contained within the current
 *          subroutine.
 *          If ISHIFT=0, NP is the number of shifts the user needs
 *          to provide via reverse communication. 0 < NP < NCV-NEV.
 *          NP may be less than NCV-NEV since a leading block of the current
 *          upper Hessenberg matrix has split off and contains "unwanted"
 *          Ritz values.
 *          Upon termination of the IRA iteration, NP contains the number
 *          of "converged" wanted Ritz values.
 *
 *  IUPD    Integer.  (INPUT)
 *          IUPD .EQ. 0: use explicit restart instead implicit update.
 *          IUPD .NE. 0: use implicit update.
 *
 *  V       Complex  N by (NEV+NP) array.  (INPUT/OUTPUT)
 *          The Arnoldi basis vectors are returned in the first NEV
 *          columns of V.
 *
 *  LDV     Integer.  (INPUT)
 *          Leading dimension of V exactly as declared in the calling
 *          program.
 *
 *  H       Complex  (NEV+NP) by (NEV+NP) array.  (OUTPUT)
 *          H is used to store the generated upper Hessenberg matrix
 *
 *  LDH     Integer.  (INPUT)
 *          Leading dimension of H exactly as declared in the calling
 *          program.
 *
 *  RITZ    Complex  array of length NEV+NP.  (OUTPUT)
 *          RITZ(1:NEV)  contains the computed Ritz values of OP.
 *
 *  BOUNDS  Complex  array of length NEV+NP.  (OUTPUT)
 *          BOUNDS(1:NEV) contain the error bounds corresponding to
 *          the computed Ritz values.
 *
 *  Q       Complex  (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
 *          Private (replicated) work array used to accumulate the
 *          rotation in the shift application step.
 *
 *  LDQ     Integer.  (INPUT)
 *          Leading dimension of Q exactly as declared in the calling
 *          program.
 *
 *  WORKL   Complex  work array of length at least
 *          (NEV+NP)**2 + 3*(NEV+NP).  (WORKSPACE)
 *          Private (replicated) array on each PE or array allocated on
 *          the front end.  It is used in shifts calculation, shifts
 *          application and convergence checking.
 *
 *
 *  IPNTR   Integer array of length 3.  (OUTPUT)
 *          Pointer to mark the starting locations in the WORKD for
 *          vectors used by the Arnoldi iteration.
 *          -------------------------------------------------------------
 *          IPNTR(1): pointer to the current operand vector X.
 *          IPNTR(2): pointer to the current result vector Y.
 *          IPNTR(3): pointer to the vector B * X when used in the
 *                    shift-and-invert mode.  X is the current operand.
 *          -------------------------------------------------------------
 *
 *  WORKD   Complex  work array of length 3*N.  (WORKSPACE)
 *          Distributed array to be used in the basic Arnoldi iteration
 *          for reverse communication.  The user should not use WORKD
 *          as temporary workspace during the iteration !!!!!!!!!!
 *          See Data Distribution Note in CNAUPD.
 *
 *  RWORK   Real    work array of length  NEV+NP ( WORKSPACE)
 *          Private (replicated) array on each PE or array allocated on
 *          the front end.
 *
 *  INFO    Integer.  (INPUT/OUTPUT)
 *          If INFO .EQ. 0, a randomly initial residual vector is used.
 *          If INFO .NE. 0, RESID contains the initial residual vector,
 *                          possibly from a previous run.
 *          Error flag on output.
 *          =     0: Normal return.
 *          =     1: Maximum number of iterations taken.
 *                   All possible eigenvalues of OP has been found.
 *                   NP returns the number of converged Ritz values.
 *          =     2: No shifts could be applied.
 *          =    -8: Error return from LAPACK eigenvalue calculation;
 *                   This should never happen.
 *          =    -9: Starting vector is zero.
 *          = -9999: Could not build an Arnoldi factorization.
 *                   Size that was built in returned in NP.
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
 *  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
 *     Restarted Arnoldi Iteration", Rice University Technical Report
 *     TR95-13, Department of Computational and Applied Mathematics.
 *
 * \Routines called:
 *     cgetv0  ARPACK initial vector generation routine.
 *     cnaitr  ARPACK Arnoldi factorization routine.
 *     cnapps  ARPACK application of implicit shifts routine.
 *     cneigh  ARPACK compute Ritz values and error bounds routine.
 *     cngets  ARPACK reorder Ritz values and error bounds routine.
 *     csortc  ARPACK sorting routine.
 *     ivout   ARPACK utility routine that prints int32_ts.
 *     arscnd  ARPACK utility routine for timing.
 *     cmout   ARPACK utility routine that prints matrices
 *     cvout   ARPACK utility routine that prints vectors.
 *     svout   ARPACK utility routine that prints vectors.
 *     slamch  LAPACK routine that determines machine constants.
 *     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     ccopy   Level 1 BLAS that copies one vector to another .
 *     cdotc   Level 1 BLAS that computes the scalar product of two vectors.
 *     cswap   Level 1 BLAS that swaps two vectors.
 *     scnrm2  Level 1 BLAS that computes the norm of a vector.
 *
 * \Author
 *     Danny Sorensen               Phuong Vu
 *     Richard Lehoucq              CRPC / Rice Universitya
 *     Chao Yang                    Houston, Texas
 *     Dept. of Computational &
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * \SCCS Information: @(#)
 * FILE: naup2.F   SID: 2.6   DATE OF SID: 06/01/00   RELEASE: 2
 *
 * \Remarks
 *     1. None
 *
 * \EndLib
 */

int cnaup2_(int32_t *ido, char *bmat, int32_t *n, char *which, int32_t *nev, int32_t *np,
     float *tol, complex *resid, int32_t *mode, int32_t *iupd, int32_t *ishift, int32_t *mxiter,
     complex *v, int32_t *ldv, complex *h, int32_t *ldh, complex *ritz, complex *bounds,
     complex *q, int32_t *ldq, complex *workl, int32_t *ipntr, complex *workd, float *rwork,
     int32_t *info)
{
    /* System generated locals */
    int32_t h_dim1, h_offset, q_dim1, q_offset, v_dim1, v_offset, i__1, i__2, 
	    i__3;
    float r__1, r__2, r__3, r__4;
    double d__1;
    complex q__1;

    /* Builtin functions */
    double pow_dd(double *, double *), r_imag(complex *);
    
    double sqrt(double);

    /* Local variables */
    int32_t i, j;
    static float t0, t1, t2, t3;
    int32_t kp[3];
    static int32_t np0, nev0;
    static float eps23;
    int32_t ierr;
    static int32_t iter;
    static bool getv0;
    static bool cnorm;
    static int32_t nconv;
    float rtemp;
    static bool initv;
    static float rnorm;
    static int32_t nevbef;
    static bool update, ushift;
    static int32_t kplusp, msglvl;
    int32_t nptemp;
    char wprime[2];
    complex cmpnorm;

     /* --------------------- */
     /* Executable Statements */
     /* --------------------- */

    /* Parameter adjustments */
    --workd;
    --resid;
    --rwork;
    --workl;
    --bounds;
    --ritz;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h -= h_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --ipntr;

    /* Function Body */
    if (*ido == 0) {

	arscnd_(&t0);

	msglvl = debug_1.mcaup2;

	nev0 = *nev;
	np0 = *np;

        /* ----------------------------------- */
        /* kplusp is the bound on the largest  */
        /*        Lanczos factorization built. */
        /* nconv is the current number of      */
        /*        "converged" eigenvalues.     */
        /* iter is the counter on the current  */
        /*      iteration step.                */
        /* ----------------------------------- */

	kplusp = *nev + *np;
	nconv = 0;
	iter = 0;

        /* ------------------------------- */
        /* Get machine dependent constant. */
        /* ------------------------------- */

	eps23 = slamch_("Epsilon-Machine");
	d__1 = (double) eps23;
	eps23 = pow_dd(&d__1, &d_23);

        /* ------------------------------------- */
        /* Set flags for computing the first NEV */
        /* steps of the Arnoldi factorization.   */
        /* ------------------------------------- */

	getv0 = true;
	update = false;
	ushift = false;
	cnorm = false;

	if (*info != 0) {

           /* ------------------------------------------ */
           /* User provides the initial residual vector. */
           /* ------------------------------------------ */

	    initv = true;
	    *info = 0;
	} else {
	    initv = false;
	}
    }

     /* ------------------------------------------- */
     /* Get a possibly random starting vector and   */
     /* force it into the range of the operator OP. */
     /* ------------------------------------------- */

/* L10: */

    if (getv0) {
	cgetv0_(ido, bmat, &c__1, &initv, n, &c__1, &v[v_offset], ldv, &resid[1], &rnorm, &ipntr[1], &workd[1], info);

	if (*ido != 99) {
	    goto L9000;
	}

	if (rnorm == 0.f) {

           /* --------------------------------------- */
           /* The initial vector is zero. Error exit. */
           /* --------------------------------------- */

	    *info = -9;
	    goto L1100;
	}
	getv0 = false;
	*ido = 0;
    }

     /* --------------------------------- */
     /* Back from reverse communication : */
     /* continue with update step         */
     /* --------------------------------- */

    if (update) {
	goto L20;
    }

     /* ----------------------------------------- */
     /* Back from computing user specified shifts */
     /* ----------------------------------------- */

    if (ushift) {
	goto L50;
    }

     /* ----------------------------------- */
     /* Back from computing residual norm   */
     /* at the end of the current iteration */
     /* ----------------------------------- */

    if (cnorm) {
	goto L100;
    }

     /* -------------------------------------------------------- */
     /* Compute the first NEV steps of the Arnoldi factorization */
     /* -------------------------------------------------------- */

    cnaitr_(ido, bmat, n, &c__0, nev, mode, &resid[1], &rnorm, &v[v_offset], ldv, &h[h_offset], ldh, &ipntr[1], &workd[1], info);

    if (*ido != 99) {
	goto L9000;
    }

    if (*info > 0) {
	*np = *info;
	*mxiter = iter;
	*info = -9999;
	goto L1200;
    }

     /* ------------------------------------------------------------ */
     /*                                                              */
     /*           M A I N  ARNOLDI  I T E R A T I O N  L O O P       */
     /*           Each iteration implicitly restarts the Arnoldi     */
     /*           factorization in place.                            */
     /*                                                              */
     /* ------------------------------------------------------------ */

L1000:

    ++iter;

    if (msglvl > 0) {
	ivout_(&c__1, &iter, &debug_1.ndigit, "_naup2: **** Start of major iteration number ****");
    }

        /* --------------------------------------------------------- */
        /* Compute NP additional steps of the Arnoldi factorization. */
        /* Adjust NP since NEV might have been updated by last call  */
        /* to the shift application routine cnapps.                  */
        /* --------------------------------------------------------- */

    *np = kplusp - *nev;

    if (msglvl > 1) {
	ivout_(&c__1, nev, &debug_1.ndigit, "_naup2: The length of the current Arnoldi factorization");
	ivout_(&c__1, np, &debug_1.ndigit, "_naup2: Extend the Arnoldi factorization by");
    }

        /* --------------------------------------------------------- */
        /* Compute NP additional steps of the Arnoldi factorization. */
        /* --------------------------------------------------------- */

    *ido = 0;
L20:
    update = true;

    cnaitr_(ido, bmat, n, nev, np, mode, &resid[1], &rnorm, &v[v_offset], ldv,&h[h_offset], ldh, &ipntr[1], &workd[1], info);

    if (*ido != 99) {
	goto L9000;
    }

    if (*info > 0) {
	*np = *info;
	*mxiter = iter;
	*info = -9999;
	goto L1200;
    }
    update = false;

    if (msglvl > 1) {
	svout_(&c__1, &rnorm, &debug_1.ndigit, "_naup2: Corresponding B-norm of the residual");
    }

        /* ------------------------------------------------------ */
        /* Compute the eigenvalues and corresponding error bounds */
        /* of the current upper Hessenberg matrix.                */
        /* ------------------------------------------------------ */

    cneigh_(&rnorm, &kplusp, &h[h_offset], ldh, &ritz[1], &bounds[1], &q[q_offset], ldq, &workl[1], &rwork[1], &ierr);

    if (ierr != 0) {
	*info = -8;
	goto L1200;
    }

        /* ------------------------------------------------- */
        /* Select the wanted Ritz values and their bounds    */
        /* to be used in the convergence test.               */
        /* The wanted part of the spectrum and corresponding */
        /* error bounds are in the last NEV loc. of RITZ,    */
        /* and BOUNDS respectively.                          */
        /* ------------------------------------------------- */

    *nev = nev0;
    *np = np0;

        /* ------------------------------------------------ */
        /* Make a copy of Ritz values and the corresponding */
        /* Ritz estimates obtained from cneigh.             */
        /* ------------------------------------------------ */

/* Computing 2nd power */
    i__1 = kplusp;
    ccopy_(&kplusp, &ritz[1], &c__1, &workl[i__1 * i__1 + 1], &c__1);
/* Computing 2nd power */
    i__1 = kplusp;
    ccopy_(&kplusp, &bounds[1], &c__1, &workl[i__1 * i__1 + kplusp + 1], &
	    c__1);

        /* ------------------------------------------------- */
        /* Select the wanted Ritz values and their bounds    */
        /* to be used in the convergence test.               */
        /* The wanted part of the spectrum and corresponding */
        /* bounds are in the last NEV loc. of RITZ           */
        /* BOUNDS respectively.                              */
        /* ------------------------------------------------- */

    cngets_(ishift, which, nev, np, &ritz[1], &bounds[1]);

        /* ---------------------------------------------------------- */
        /* Convergence test: currently we use the following criteria. */
        /* The relative accuracy of a Ritz value is considered        */
        /* acceptable if:                                             */
        /*                                                            */
        /* error_bounds(i) .le. tol*max(eps23, magnitude_of_ritz(i)). */
        /*                                                            */
        /* ---------------------------------------------------------- */

    nconv = 0;

    i__1 = *nev;
    for (i = 1; i <= i__1; ++i) {
/* Computing MAX */
	i__2 = *np + i;
	r__3 = ritz[i__2].r;
	r__4 = r_imag(&ritz[*np + i]);
	r__1 = eps23, r__2 = slapy2_(&r__3, &r__4);
	rtemp = dmax(r__1,r__2);
	i__2 = *np + i;
	r__1 = bounds[i__2].r;
	r__2 = r_imag(&bounds[*np + i]);
	if (slapy2_(&r__1, &r__2) <= *tol * rtemp) {
	    ++nconv;
	}
/* L25: */
    }

    if (msglvl > 2) {
	kp[0] = *nev;
	kp[1] = *np;
	kp[2] = nconv;
	ivout_(&c__3, kp, &debug_1.ndigit, "_naup2: NEV, NP, NCONV are");
	cvout_(&kplusp, &ritz[1], &debug_1.ndigit, "_naup2: The eigenvalues of H");
	cvout_(&kplusp, &bounds[1], &debug_1.ndigit, "_naup2: Ritz estimates of the current NCV Ritz values");
    }

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
    i__1 = nptemp;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	if (bounds[i__2].r == 0.f && bounds[i__2].i == 0.f) {
	    --(*np);
	    ++(*nev);
	}
/* L30: */
    }

    if (nconv >= nev0 || iter > *mxiter || *np == 0) {

	if (msglvl > 4) {
/* Computing 2nd power */
	    i__1 = kplusp;
	    cvout_(&kplusp, &workl[i__1 * i__1 + 1], &debug_1.ndigit, "_naup2: Eigenvalues computed by _neigh:");
/* Computing 2nd power */
	    i__1 = kplusp;
	    cvout_(&kplusp, &workl[i__1 * i__1 + kplusp + 1],&debug_1.ndigit, "_naup2: Ritz estimates computed by _neigh:");
	}

           /* ---------------------------------------------- */
           /* Prepare to exit. Put the converged Ritz values */
           /* and corresponding bounds in RITZ(1:NCONV) and  */
           /* BOUNDS(1:NCONV) respectively. Then sort. Be    */
           /* careful when NCONV > NP                        */
           /* ---------------------------------------------- */

           /* ---------------------------------------- */
           /*  Use h( 3,1 ) as storage to communicate  */
           /*  rnorm to cneupd if needed               */
           /* ---------------------------------------- */
	i__1 = h_dim1 + 3;
	q__1.r = rnorm, q__1.i = 0.f;
	h[i__1].r = q__1.r, h[i__1].i = q__1.i;

           /* -------------------------------------------- */
           /* Sort Ritz values so that converged Ritz      */
           /* values appear within the first NEV locations */
           /* of ritz and bounds, and the most desired one */
           /* appears at the front.                        */
           /* -------------------------------------------- */

	if (strcmp(which, "LM") == 0) {
	    strcpy(wprime, "SM");
	}
	if (strcmp(which, "SM") == 0) {
	    strcpy(wprime, "LM");
	}
	if (strcmp(which, "LR") == 0) {
	    strcpy(wprime, "SR");
	}
	if (strcmp(which, "SR") == 0) {
	    strcpy(wprime, "LR");
	}
	if (strcmp(which, "LI") == 0) {
	    strcpy(wprime, "SI");
	}
	if (strcmp(which, "SI") == 0) {
	    strcpy(wprime, "LI");
	}

	csortc_(wprime, &c_true, &kplusp, &ritz[1], &bounds[1]);

           /* ------------------------------------------------ */
           /* Scale the Ritz estimate of each Ritz value       */
           /* by 1 / max(eps23, magnitude of the Ritz value).  */
           /* ------------------------------------------------ */

	i__1 = nev0;
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	    i__2 = j;
	    r__3 = ritz[i__2].r;
	    r__4 = r_imag(&ritz[j]);
	    r__1 = eps23, r__2 = slapy2_(&r__3, &r__4);
	    rtemp = dmax(r__1,r__2);
	    i__2 = j;
	    i__3 = j;
	    q__1.r = bounds[i__3].r / rtemp, q__1.i = bounds[i__3].i / rtemp;
	    bounds[i__2].r = q__1.r, bounds[i__2].i = q__1.i;
/* L35: */
	}

           /* ------------------------------------------------- */
           /* Sort the Ritz values according to the scaled Ritz */
           /* estimates.  This will push all the converged ones */
           /* towards the front of ritz, bounds (in the case    */
           /* when NCONV < NEV.)                                */
           /* ------------------------------------------------- */

	strcpy(wprime, "LM");
	csortc_(wprime, &c_true, &nev0, &bounds[1], &ritz[1]);

           /* -------------------------------------------- */
           /* Scale the Ritz estimate back to its original */
           /* value.                                       */
           /* -------------------------------------------- */

	i__1 = nev0;
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	    i__2 = j;
	    r__3 = ritz[i__2].r;
	    r__4 = r_imag(&ritz[j]);
	    r__1 = eps23, r__2 = slapy2_(&r__3, &r__4);
	    rtemp = dmax(r__1,r__2);
	    i__2 = j;
	    i__3 = j;
	    q__1.r = rtemp * bounds[i__3].r, q__1.i = rtemp * bounds[i__3].i;
	    bounds[i__2].r = q__1.r, bounds[i__2].i = q__1.i;
/* L40: */
	}

           /* --------------------------------------------- */
           /* Sort the converged Ritz values again so that  */
           /* the "threshold" value appears at the front of */
           /* ritz and bound.                               */
           /* --------------------------------------------- */

	csortc_(which, &c_true, &nconv, &ritz[1], &bounds[1]);

	if (msglvl > 1) {
	    cvout_(&kplusp, &ritz[1], &debug_1.ndigit, "_naup2: Sorted eigenvalues");
	    cvout_(&kplusp, &bounds[1], &debug_1.ndigit, "_naup2: Sorted ritz estimates.");
	}

           /* ---------------------------------- */
           /* Max iterations have been exceeded. */
           /* ---------------------------------- */

	if (iter > *mxiter && nconv < nev0) {
	    *info = 1;
	}

           /* ------------------- */
           /* No shifts to apply. */
           /* ------------------- */

	if (*np == 0 && nconv < nev0) {
	    *info = 2;
	}

	*np = nconv;
	goto L1100;

    } else if (nconv < nev0 && *ishift == 1) {

           /* ----------------------------------------------- */
           /* Do not have all the requested eigenvalues yet.  */
           /* To prevent possible stagnation, adjust the size */
           /* of NEV.                                         */
           /* ----------------------------------------------- */

	nevbef = *nev;
/* Computing MIN */
	i__1 = nconv, i__2 = *np / 2;
	*nev += min(i__1,i__2);
	if (*nev == 1 && kplusp >= 6) {
	    *nev = kplusp / 2;
	} else if (*nev == 1 && kplusp > 3) {
	    *nev = 2;
	}
	*np = kplusp - *nev;

           /* ------------------------------------- */
           /* If the size of NEV was just increased */
           /* resort the eigenvalues.               */
           /* ------------------------------------- */

	if (nevbef < *nev) {
	    cngets_(ishift, which, nev, np, &ritz[1], &bounds[1]);
	}

    }

    if (msglvl > 0) {
	ivout_(&c__1, &nconv, &debug_1.ndigit, "_naup2: no. of \"converged\" Ritz values at this iter.");
	if (msglvl > 1) {
	    kp[0] = *nev;
	    kp[1] = *np;
	    ivout_(&c__2, kp, &debug_1.ndigit, "_naup2: NEV and NP are");
	    cvout_(nev, &ritz[*np + 1], &debug_1.ndigit, "_naup2: \"wanted\" Ritz values ");
	    cvout_(nev, &bounds[*np + 1], &debug_1.ndigit, "_naup2: Ritz estimates of the \"wanted\" values ");
	}
    }

    if (*ishift == 0) {

           /* ----------------------------------------------------- */
           /* User specified shifts: pop back out to get the shifts */
           /* and return them in the first 2*NP locations of WORKL. */
           /* ----------------------------------------------------- */

	ushift = true;
	*ido = 3;
	goto L9000;
    }
L50:
    ushift = false;

    if (*ishift != 1) {

            /* -------------------------------- */
            /* Move the NP shifts from WORKL to */
            /* RITZ, to free up WORKL           */
            /* for non-exact shift case.        */
            /* -------------------------------- */

	ccopy_(np, &workl[1], &c__1, &ritz[1], &c__1);
    }

    if (msglvl > 2) {
	ivout_(&c__1, np, &debug_1.ndigit, "_naup2: The number of shifts to apply ");
	cvout_(np, &ritz[1], &debug_1.ndigit, "_naup2: values of the shifts");
	if (*ishift == 1) {
	    cvout_(np, &bounds[1], &debug_1.ndigit, "_naup2: Ritz estimates of the shifts");
	}
    }

        /* ------------------------------------------------------- */
        /* Apply the NP implicit shifts by QR bulge chasing.       */
        /* Each shift is applied to the whole upper Hessenberg     */
        /* matrix H.                                               */
        /* The first 2*N locations of WORKD are used as workspace. */
        /* ------------------------------------------------------- */

    cnapps_(n, nev, np, &ritz[1], &v[v_offset], ldv, &h[h_offset], ldh, &resid[1], &q[q_offset], ldq, &workl[1], &workd[1]);

        /* ------------------------------------------- */
        /* Compute the B-norm of the updated residual. */
        /* Keep B*RESID in WORKD(1:N) to be used in    */
        /* the first step of the next call to cnaitr.  */
        /* ------------------------------------------- */

    cnorm = true;
    arscnd_(&t2);
    if (*bmat == 'G') {
	++timing_1.nbx;
	ccopy_(n, &resid[1], &c__1, &workd[*n + 1], &c__1);
	ipntr[1] = *n + 1;
	ipntr[2] = 1;
	*ido = 2;

           /* -------------------------------- */
           /* Exit in order to compute B*RESID */
           /* -------------------------------- */

	goto L9000;
    } else if (*bmat == 'I') {
	ccopy_(n, &resid[1], &c__1, &workd[1], &c__1);
    }

L100:

        /* -------------------------------- */
        /* Back from reverse communication; */
        /* WORKD(1:N) := B*RESID            */
        /* -------------------------------- */

    if (*bmat == 'G') {
	arscnd_(&t3);
	timing_1.tmvbx += t3 - t2;
    }

    if (*bmat == 'G') {
	cdotc_(&q__1, n, &resid[1], &c__1, &workd[1], &c__1);
	cmpnorm.r = q__1.r, cmpnorm.i = q__1.i;
	r__1 = cmpnorm.r;
	r__2 = r_imag(&cmpnorm);
	rnorm = sqrt(slapy2_(&r__1, &r__2));
    } else if (*bmat == 'I') {
	rnorm = scnrm2_(n, &resid[1], &c__1);
    }
    cnorm = false;

    if (msglvl > 2) {
	svout_(&c__1, &rnorm, &debug_1.ndigit, "_naup2: B-norm of residual for compressed factorization");
	cmout_(nev, nev, &h[h_offset], ldh, &debug_1.ndigit, "_naup2: Compressed upper Hessenberg matrix H");
    }

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
     /* Error Exit */
     /* ---------- */

    arscnd_(&t1);
    timing_1.tcaup2 = t1 - t0;

L9000:

     /* ------------- */
     /* End of cnaup2 */
     /* ------------- */

    return 0;
} /* cnaup2_ */

