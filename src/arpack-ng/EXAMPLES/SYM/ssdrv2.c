/* EXAMPLES\SYM\ssdrv2.f -- translated by f2c (version 20100827). */

#include <stdlib.h>
#include "arpack.h"

int ssdrv2_av_(const int nx, float* v, float* w);

extern int smout_(const int, const int, const float*, const int, const int, const char*);

static int c__1 = 1;

/**
 * \BeginDoc
 *
 *     Program to illustrate the idea of reverse communication
 *     in shift and invert mode for a standard symmetric eigenvalue
 *     problem.  The following program uses the two LAPACK subroutines
 *     sgttrf.f and sgttrs.f to factor and solve a tridiagonal system of
 *     equations.
 *
 *     We implement example two of ex-sym.doc in DOCUMENTS directory
 *
 * \Example-2
 *     ... Suppose we want to solve A*x = lambda*x in shift-invert mode,
 *         where A is derived from the central difference discretization
 *         of the 1-dimensional Laplacian on [0,1]  with zero Dirichlet
 *         boundary condition.
 *     ... OP = (inv[A - sigma*I]) and  B = I.
 *     ... Use mode 3 of SSAUPD.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Routines called:
 *     ssaupd  ARPACK reverse communication interface routine.
 *     sseupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     sgttrf  LAPACK tridiagonal factorization routine.
 *     sgttrs  LAPACK tridiagonal solve routine.
 *     saxpy   saxpy   Level 1 BLAS that computes y <- alpha*x+y.
 *     snrm2   Level 1 BLAS that computes the norm of a vector.
 *     av      Matrix vector multiplication routine that computes A*x.
 *
 * \EndLib
 */
int main()
{
    /* System generated locals */
    float r__1;

    /* Local variables */
    float d[50] /* (2 * MAXNCV) */;
    float h2, sigma;

    int j;
    int ierr, nconv;
    int ishfts, maxitr, mode;
    int ipiv[256];
    int ipntr[11];
    int iparam[11];
    logical select[25];
    logical rvec;

    /* Define maximum dimensions for all arrays. */

    const int MAXN   = 256; /* Maximum dimension of the A allowed. */
    const int MAXNEV =  10; /* Maximum NEV allowed */
    const int MAXNCV =  25; /* Maximum NCV allowed */

    /* -------------------------------------------------- */
    /* The number N is the dimension of the matrix.  A    */
    /* standard eigenvalue problem is solved (BMAT = 'I'. */
    /* NEV is the number of eigenvalues (closest to       */
    /* SIGMA) to be approximated.  Since the shift-invert */
    /* mode is used, WHICH is set to 'LM'.  The user can  */
    /* modify NEV, NCV, SIGMA to solve problems of        */
    /* different sizes, and to get different parts of the */
    /* spectrum.  However, The following conditions must  */
    /* be satisfied:                                      */
    /*                   N <= MAXN,                       */
    /*                 NEV <= MAXNEV,                     */
    /*             NEV + 1 <= NCV <= MAXNCV               */
    /* -------------------------------------------------- */

    int n = 100;
    int nev = 4;
    int ncv = 10;
    if (n > MAXN)
    {
        printf(" ERROR with _SDRV2: N is greater than MAXN \n");
        return 0;
    }
    else if (nev > MAXNEV)
    {
        printf(" ERROR with _SDRV2: NEV is greater than MAXNEV \n");
        return 0;
    }
    else if (ncv > MAXNCV)
    {
        printf(" ERROR with _SDRV2: NCV is greater than MAXNCV \n");
        return 0;
    }

    char* bmat = "I";
    char* which = "LM";
    sigma = 0.0f;

    /* ------------------------------------------------ */
    /* The work array WORKL is used in SSAUPD as        */
    /* workspace.  Its dimension LWORKL is set as       */
    /* illustrated below.  The parameter TOL determines */
    /* the stopping criterion.  If TOL<=0, machine      */
    /* precision is used.  The variable IDO is used for */
    /* reverse communication and is initially set to 0. */
    /* Setting INFO=0 indicates that a random vector is */
    /* generated in SSAUPD to start the Arnoldi         */
    /* iteration.                                       */
    /* ------------------------------------------------ */

    int lworkl = ncv * (ncv + 8);
    float tol = 0.0f;
    int ido = 0;
    int info = 0;

    float* ad = (float*)malloc(n * sizeof(float));
    float* adl = (float*)malloc(n * sizeof(float));
    float* adu = (float*)malloc(n * sizeof(float));
    float* adu2 = (float*)malloc(n * sizeof(float));

    float* ax = (float*)malloc(n * sizeof(float));
    float* resid = (float*)malloc(n * sizeof(float));
    float* v = (float*)malloc(n * ncv * sizeof(float));
    float* workl = (float*)malloc(lworkl * sizeof(float));
    float* workd = (float*)malloc(3 * n * sizeof(float));

    /* ------------------------------------------------- */
    /* This program uses exact shifts with respect to    */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 3 of SSAUPD is used     */
    /* (IPARAM(7) = 3).  All these options may be        */
    /* changed by the user. For details, see the         */
    /* documentation in SSAUPD.                          */
    /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode = 3;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

    /* --------------------------------------------------- */
    /* Call LAPACK routine to factor (A-SIGMA*I), where A  */
    /* is the 1-d Laplacian.                               */
    /* --------------------------------------------------- */

    h2 = 1.0f / (float) ((n + 1) * (n + 1));
    for (j = 1; j <= n; ++j)
    {
        ad[j - 1] = 2.0f / h2 - sigma;
        adl[j - 1] = -1.0f / h2;
    }
    scopy_(&n, adl, &c__1, adu, &c__1);
    sgttrf_(&n, adl, ad, adu, adu2, ipiv, &ierr);
    if (ierr != 0)
    {
        printf(" \n");
        printf(" Error with _gttrf in SDRV2.\n");
        printf(" \n");
        return 0;
    }

    /* ----------------------------------------- */
    /* M A I N   L O O P (Reverse communication) */
    /* ----------------------------------------- */

L10:

    /* ------------------------------------------- */
    /* Repeatedly call the routine SSAUPD and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    ssaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);

    if (ido == -1 || ido == 1)
    {
        /* -------------------------------------- */
        /* Perform y <-- OP*x = inv[A-sigma*I]*x. */
        /* The user only need the linear system   */
        /* solver here that takes workd(ipntr(1)) */
        /* as input, and returns the result to    */
        /* workd(ipntr(2)).                       */
        /* -------------------------------------- */

        scopy_(&n, &workd[ipntr[0] - 1], &c__1, &workd[ipntr[1] - 1], &c__1);

        sgttrs_("N", &n, &c__1, adl, ad, adu, adu2, ipiv, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" Error with _gttrs in _SDRV2. \n");
            printf(" \n");
            goto EXIT;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call SSAUPD again. */
        /* --------------------------------------- */

        goto L10;
    }

    /* -------------------------------------- */
    /* Either we have convergence or there is */
    /* an error.                              */
    /* -------------------------------------- */

    if (info < 0)
    {
        /* -------------------------- */
        /* Error message.  Check the  */
        /* documentation in SSAUPD    */
        /* -------------------------- */

        printf(" \n");
        printf(" Error with _saupd info = %d\n", info);
        printf(" Check the documentation of _saupd.\n");
        printf(" \n");

        goto EXIT;
    }

    /* ----------------------------------------- */
    /* No fatal errors occurred.                 */
    /* Post-Process using SSEUPD.                */
    /*                                           */
    /* Computed eigenvalues may be extracted.    */
    /*                                           */
    /* Eigenvectors may also be computed now if  */
    /* desired.  (indicated by rvec = .true.)    */
    /* ----------------------------------------- */

    rvec = TRUE_;

    sseupd_(&rvec, "A", select, d, v, &n, &sigma, bmat, &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &ierr);

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
        /* Check the documentation of SSEUPD. */
        /* ---------------------------------- */

        printf(" \n");
        printf(" Error with _seupd info = %d\n", ierr);
        printf(" Check the documentation of _seupd \n");
        printf(" \n");

        goto EXIT;
    }

    nconv = iparam[4];
    for (j = 1; j <= nconv; ++j)
    {
        int k = (j - 1) * n;

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

        ssdrv2_av_(n, &v[k], ax);
        r__1 = -d[j - 1];
        saxpy_(&n, &r__1, &v[k], &c__1, ax, &c__1);
        d[j + 24] = snrm2_(&n, ax, &c__1);
        d[j + 24] /= (r__1 = d[j - 1], dabs(r__1));

    }

    /* ----------------------------- */
    /* Display computed residuals    */
    /* ----------------------------- */

    smout_(nconv, 2, d, MAXNCV, -6, "Ritz values and relative residuals");
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
    printf(" _SDRV2 \n");
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

    free(ax);
    free(resid);
    free(v);
    free(workl);
    free(workd);

    /* ------------------------- */
    /* Done with program ssdrv2. */
    /* ------------------------- */

    return ierr;
}

/* ------------------------------------------------------------------------ */
/*     Matrix vector subroutine */
/*     where the matrix is the 1 dimensional discrete Laplacian on */
/*     the interval [0,1] with zero Dirichlet boundary condition. */

int ssdrv2_av_(const int n, float *v, float *w)
{
    /* System generated locals */
    int i__1;
    float r__1;

    /* Local variables */
    int j;
    float h2;

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    w[1] = v[1] * 2.0f - v[2];
    i__1 = n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = -v[j - 1] + v[j] * 2.0f - v[j + 1];
    }
    j = n;
    w[j] = -v[j - 1] + v[j] * 2.0f;

    /*     Scale the vector w by (1 / h^2). */

    h2 = 1.0f / (float) ((n + 1) * (n + 1));
    r__1 = 1.0f / h2;
    sscal_(&n, &r__1, &w[1], &c__1);
    return 0;
} /* av_ */

