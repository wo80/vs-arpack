/* EXAMPLES\SYM\dsdrv6.f -- translated by f2c (version 20100827). */

#include <stdlib.h>
#include "arpack_internal.h"

/**
 * \BeginDoc
 *
 *     Program to illustrate the idea of reverse communication
 *     in Cayley mode for a generalized symmetric eigenvalue
 *     problem.  The following program uses the two LAPACK subroutines
 *     dgttrf.f and dgttrs.f to factor and solve a tridiagonal system of
 *     equations.
 *
 *     We implement example six of ex-sym.doc in DOCUMENTS directory
 *
 * \Example-6
 *     ... Suppose we want to solve A*x = lambda*M*x in inverse mode,
 *         where A and M are obtained by the finite element of the
 *         1-dimensional discrete Laplacian
 *                             d^2u / dx^2
 *         on the interval [0,1] with zero Dirichlet boundary condition
 *         using piecewise linear elements.
 *
 *     ... OP = (inv[A-sigma*M])*(A+sigma*M)  and  B = M.
 *
 *     ... Use mode 5 of DSAUPD.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Routines called:
 *     dsaupd  ARPACK reverse communication interface routine.
 *     dseupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     dgttrf  LAPACK tridiagonal factorization routine.
 *     dgttrs  LAPACK tridiagonal solve routine.
 *     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
 *     dcopy   Level 1 BLAS that copies one vector to another.
 *     dscal   Level 1 BLAS that scales a vector by a scalar.
 *     dnrm2   Level 1 BLAS that computes the norm of a vector.
 *     av      Matrix vector multiplication routine that computes A*x.
 *     mv      Matrix vector multiplication routine that computes M*x.
 *
 * \EndLib
 */
int dsdrv6()
{
    /* System generated locals */
    double d__1;

    /* Local variables */
    double d[50] /* (2 * MAXNCV) */;
    double r1, r2, h, sigma;

    int j;
    int ierr, nconv;
    int ishfts, maxitr, mode;
    int ipiv[256];
    int ipntr[11];
    int iparam[11];
    bool select[25];
    bool rvec;

    /* Define maximum dimensions for all arrays. */

    const int MAXN   = 256; /* Maximum dimension of the A allowed. */
    const int MAXNEV =  10; /* Maximum NEV allowed */
    const int MAXNCV =  25; /* Maximum NCV allowed */

    /* ------------------------------------------------ */
    /* The number N is the dimension of the matrix. A   */
    /* generalized eigenvalue problem is solved (BMAT = */
    /* 'G'.) NEV is the number of eigenvalues to be     */
    /* approximated.  Since the Cayley mode is used,    */
    /* WHICH is set to 'LM'.  The user can modify NEV,  */
    /* NCV, SIGMA to solve problems of different sizes, */
    /* and to get different parts of the spectrum.      */
    /* However, The following conditions must be        */
    /* satisfied:                                       */
    /*                 N <= MAXN,                       */
    /*               NEV <= MAXNEV,                     */
    /*           NEV + 1 <= NCV <= MAXNCV               */
    /* ------------------------------------------------ */

    int n = 100;
    int nev = 4;
    int ncv = 20;
    if (n > 256)
    {
        printf(" ERROR with _SDRV6: N is greater than MAXN \n");
        return 0;
    }
    else if (nev > 10)
    {
        printf(" ERROR with _SDRV6: NEV is greater than MAXNEV \n");
        return 0;
    }
    else if (ncv > 25)
    {
        printf(" ERROR with _SDRV6: NCV is greater than MAXNCV \n");
        return 0;
    }
    char* bmat = "G";
    char* which = "LM";
    sigma = 150.0;

    /* ------------------------------------------------ */
    /* The work array WORKL is used in DSAUPD as        */
    /* workspace.  Its dimension LWORKL is set as       */
    /* illustrated below.  The parameter TOL determines */
    /* the stopping criterion.  If TOL<=0, machine      */
    /* precision is used.  The variable IDO is used for */
    /* reverse communication and is initially set to 0. */
    /* Setting INFO=0 indicates that a random vector is */
    /* generated in DSAUPD to start the Arnoldi         */
    /* iteration.                                       */
    /* ------------------------------------------------ */

    int lworkl = ncv * (ncv + 8);
    double tol = 0.0;
    int ido = 0;
    int info = 0;

    double* ad = (double*)malloc(n * sizeof(double));
    double* adl = (double*)malloc(n * sizeof(double));
    double* adu = (double*)malloc(n * sizeof(double));
    double* adu2 = (double*)malloc(n * sizeof(double));
    double* temp = (double*)malloc(n * sizeof(double));

    double* ax = (double*)malloc(n * sizeof(double));
    double* mx = (double*)malloc(n * sizeof(double));
    double* resid = (double*)malloc(n * sizeof(double));
    double* v = (double*)malloc(n * ncv * sizeof(double));
    double* workl = (double*)malloc(lworkl * sizeof(double));
    double* workd = (double*)malloc(3 * n * sizeof(double));

    /* ------------------------------------------------- */
    /* This program uses exact shifts with respect to    */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 5 specified in the      */
    /* documentation of DSAUPD is used (IPARAM(7) = 5).  */
    /* All these options may be changed by the user. For */
    /* details, see the documentation in DSAUPD.         */
    /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode = 5;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

    /* ---------------------------------------------------- */
    /* Call LAPACK routine to factor (A-sigma*M).  The      */
    /* stiffness matrix A is the 1-d discrete Laplacian.    */
    /* The mass matrix M is the associated mass matrix      */
    /* arising from using piecewise linear finite elements  */
    /* on the interval [0, 1].                              */
    /* ---------------------------------------------------- */

    h = 1.0 / (double) (n + 1);
    r1 = h * .66666666666666663;
    r2 = h * .16666666666666666;
    for (j = 1; j <= n; ++j)
    {
        ad[j - 1] = 2. / h - sigma * r1;
        adl[j - 1] = -1.0 / h - sigma * r2;
    }
    dcopy_(&n, adl, &c__1, adu, &c__1);
    dgttrf_(&n, adl, ad, adu, adu2, ipiv, &ierr);
    if (ierr != 0)
    {
        printf(" \n");
        printf(" Error with _gttrf in _SDRV6.\n");
        printf(" \n");
        return 0;
    }

    /* ----------------------------------------- */
    /* M A I N   L O O P (Reverse communication) */
    /* ----------------------------------------- */

L10:

    /* ------------------------------------------- */
    /* Repeatedly call the routine DSAUPD and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);

    if (ido == -1)
    {
        /* ----------------------------------------------------- */
        /* Perform  y <--- OP*x = (inv[A-SIGMA*M])*(A+SIGMA*M)*x */
        /* to force starting vector into the range of OP.  The   */
        /* user should provide his/her matrix vector (A*x, M*x)  */
        /* multiplication routines and a linear system solver    */
        /* here.  The matrix vector multiplication routine takes */
        /* workd(ipntr(1)) as the input vector.  The final       */
        /* result is returned to workd(ipntr(2)).                */
        /* ----------------------------------------------------- */

        dsdrv6_av_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
        dsdrv6_mv_(n, &workd[ipntr[0] - 1], temp);
        daxpy_(&n, &sigma, temp, &c__1, &workd[ipntr[1] - 1], &c__1);

        dgttrs_("N", &n, &c__1, adl, ad, adu, adu2, ipiv, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" Error with _gttrs in _SDRV6.\n");
            printf(" \n");
            goto EXIT;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call DSAUPD again. */
        /* --------------------------------------- */

        goto L10;
    }
    else if (ido == 1)
    {
        /* -------------------------------------------------- */
        /* Perform y <-- OP*x = inv[A-SIGMA*M]*(A+SIGMA*M)*x. */
        /* M*x has been saved in workd(ipntr(3)).  The user   */
        /* need the matrix vector multiplication (A*x)        */
        /* routine and a linear system solver here.  The      */
        /* matrix vector multiplication routine takes         */
        /* workd(ipntr(1)) as the input, and the result is    */
        /* combined with workd(ipntr(3)) to form the input    */
        /* for the linear system solver.  The final result is */
        /* returned to workd(ipntr(2)).                       */
        /* -------------------------------------------------- */

        dsdrv6_av_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
        daxpy_(&n, &sigma, &workd[ipntr[2] - 1], &c__1, &workd[ipntr[1] - 1], &c__1);
        dgttrs_("N", &n, &c__1, adl, ad, adu, adu2, ipiv, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" Error with _gttrs in _SDRV6. \n");
            printf(" \n");
            goto EXIT;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call DSAUPD again. */
        /* --------------------------------------- */

        goto L10;
    }
    else if (ido == 2)
    {
        /* ------------------------------------------ */
        /*             Perform  y <--- M*x.           */
        /* Need matrix vector multiplication routine  */
        /* here that takes workd(ipntr(1)) as input   */
        /* and returns the result to workd(ipntr(2)). */
        /* ------------------------------------------ */

        dsdrv6_mv_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call DSAUPD again. */
        /* --------------------------------------- */

        goto L10;
    }

    /* --------------------------------------- */
    /* Either we have convergence, or there is */
    /* an error.                               */
    /* --------------------------------------- */

    if (info < 0)
    {
        /* ------------------------ */
        /* Error message, check the */
        /* documentation in DSAUPD  */
        /* ------------------------ */

        printf(" \n");
        printf(" Error with _saupd info = %d\n", info);
        printf(" Check the documentation of _saupd. \n");
        printf(" \n");

        goto EXIT;
    }

    /* ----------------------------------------- */
    /* No fatal errors occurred.                 */
    /* Post-Process using DSEUPD.                */
    /*                                           */
    /* Computed eigenvalues may be extracted.    */
    /*                                           */
    /* Eigenvectors may also be computed now if  */
    /* desired.  (indicated by rvec = .true.)    */
    /* ----------------------------------------- */

    rvec = true;

    dseupd_(&rvec, "A", select, d, v, &n, &sigma, bmat, &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &ierr);

    /* -------------------------------------------- */
    /* Eigenvalues are returned in the first column */
    /* of the two dimensional array D and the       */
    /* corresponding eigenvectors are returned in   */
    /* the first NEV columns of the two dimensional */
    /* array V if requested.  Otherwise, an         */
    /* orthogonal basis for the invariant subspace  */
    /* corresponding to the eigenvalues in D is     */
    /* returned in V.                               */
    /* -------------------------------------------- */

    if (ierr != 0)
    {
        /* ---------------------------------- */
        /* Error condition:                   */
        /* Check the documentation of DSEUPD. */
        /* ---------------------------------- */

        printf(" \n");
        printf(" Error with _seupd info = %d\n", ierr);
        printf(" Check the documentation of _seupd \n");
        printf(" \n");

        goto EXIT;
    }

    /* ------------------------- */
    /* Compute the residual norm */
    /*                           */
    /*   ||  A*x - lambda*x ||   */
    /*                           */
    /* for the NCONV accurately  */
    /* computed eigenvalues and  */
    /* eigenvectors.  (iparam(5) */
    /* indicates how many are    */
    /* accurate to the requested */
    /* tolerance)                */
    /* ------------------------- */

    nconv = iparam[4];
    for (j = 1; j <= nconv; ++j)
    {
        int k = (j - 1) * n;

        dsdrv6_av_(n, &v[k], ax);
        dsdrv6_mv_(n, &v[k], mx);
        d__1 = -d[j - 1];
        daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
        d[j + 24] = dnrm2_(&n, ax, &c__1);
        d[j + 24] /= (d__1 = d[j - 1], abs(d__1));

    }

    dmout_(&nconv, &c__2, d, &c__25, &c_n6, "Ritz values and relative residuals");

    /* ---------------------------------------- */
    /* Print additional convergence information */
    /* ---------------------------------------- */

    if (info == 1)
    {
        printf(" \n");
        printf(" Maximum number of iterations reached.\n");
        printf(" \n");
    }
    else if (info == 3)
    {
        printf(" \n");
        printf(" No shifts could be applied during implicit\n");
        printf(" Arnoldi update try increasing NCV.\n");
        printf(" \n");
    }

    printf(" \n");
    printf(" _SDRV6 \n");
    printf(" ====== \n");
    printf(" \n");
    printf(" Size of the matrix is %d\n", n);
    printf(" The number of Ritz values requested is %d\n", nev);
    printf(" The number of Arnoldi vectors generated (NCV) is %d\n", ncv);
    printf(" What portion of the spectrum: %s\n", which);
    printf(" The number of converged Ritz values is %d\n", nconv);
    printf(" The number of Implicit Arnoldi update iterations taken is %d\n", iparam[2]);
    printf(" The number of OP*x is %d\n", iparam[8]);
    printf(" The convergence criterion is %e\n", tol);
    printf(" \n");

EXIT:

    free(ad);
    free(adl);
    free(adu);
    free(adu2);
    free(temp);

    free(ax);
    free(mx);
    free(resid);
    free(v);
    free(workl);
    free(workd);

    /* ------------------------- */
    /* Done with program dsdrv6. */
    /* ------------------------- */

    return ierr;
}

/* ------------------------------------------------------------------------ */
/*     Matrix vector subroutine */
/*     where the matrix used is the 1 dimensional mass matrix */
/*     arising from using the piecewise linear finite element */
/*     on the interval [0,1]. */

int dsdrv6_mv_(const int n, double *v, double *w)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    double h;
    int j;

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    w[1] = v[1] * 4.0 + v[2];
    i__1 = n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = v[j - 1] + v[j] * 4.0 + v[j + 1];
    }
    j = n;
    w[j] = v[j - 1] + v[j] * 4.0;

    /*     Scale the vector w by h. */

    h = 1.0 / ((double) (n + 1) * 6.);
    dscal_(&n, &h, &w[1], &c__1);
    return 0;
} /* mv_ */

/* ------------------------------------------------------------------------ */
/*     Matrix vector subroutine */
/*     where the matrix is the stiffness matrix obtained from the */
/*     finite element discretization of the 1-dimensional discrete Laplacian */
/*     on the interval [0,1] with zero Dirichlet boundary condition */
/*     using piecewise linear elements. */

int dsdrv6_av_(const int n, double *v, double *w)
{
    /* System generated locals */
    int i__1;
    double d__1;

    /* Local variables */
    double h;
    int j;

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    w[1] = v[1] * 2. - v[2];
    i__1 = n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = -v[j - 1] + v[j] * 2. - v[j + 1];
    }
    j = n;
    w[j] = -v[j - 1] + v[j] * 2.0;

    /*     Scale the vector w by (1/h). */

    h = 1.0 / (double) (n + 1);
    d__1 = 1.0 / h;
    dscal_(&n, &d__1, &w[1], &c__1);
    return 0;
} /* av_ */

