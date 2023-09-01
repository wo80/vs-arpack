/* EXAMPLES\SYM\ssdrv4.f -- translated by f2c (version 20100827). */

#include <stdlib.h>
#include "arpack.h"

int ssdrv4_av_(const int nx, float* v, float* w);
int ssdrv4_mv_(const int n, float* v, float* w);

extern int smout_(const int, const int, const float*, const int, const int, const char*);

static int i_one = 1;

/**
 * \BeginDoc
 *
 *     Program to illustrate the idea of reverse communication
 *     in shift and invert mode for a generalized symmetric eigenvalue
 *     problem.  The following program uses the two LAPACK subroutines
 *     sgttrf.f and sgttrs to factor and solve a tridiagonal system of
 *     equations.
 *
 *     We implement example four of ex-sym.doc in DOCUMENTS directory
 *
 * \Example-4
 *     ... Suppose we want to solve A*x = lambda*M*x in inverse mode,
 *         where A and M are obtained from the finite element discretrization
 *         of the 1-dimensional discrete Laplacian
 *                             d^2u / dx^2
 *         on the interval [0,1] with zero Dirichlet boundary condition
 *         using piecewise linear elements.
 *
 *     ... OP = (inv[A - sigma*M])*M  and  B = M.
 *
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
 *     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
 *     scopy   Level 1 BLAS that copies one vector to another.
 *     sscal   Level 1 BLAS that scales a vector by a scalar.
 *     snrm2   Level 1 BLAS that computes the norm of a vector.
 *     av      Matrix vector multiplication routine that computes A*x.
 *     mv      Matrix vector multiplication routine that computes M*x.
 *
 * \EndLib
 */
int main()
{
    /* System generated locals */
    float r__1;

    /* Local variables */
    float d[50] /* (2 * MAXNCV) */;
    float r1, r2, h, sigma;

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
    /* generalized eigenvalue problem is solved (BMAT =   */
    /* 'G'.) NEV is the number of eigenvalues (closest to */
    /* the shift SIGMA) to be approximated.  Since the    */
    /* shift-invert mode is used, WHICH is set to 'LM'.   */
    /* The user can modify NEV, NCV, SIGMA to solve       */
    /* problems of different sizes, and to get different  */
    /* parts of the spectrum. However, The following      */
    /* conditions must be satisfied:                      */
    /*                   N <= MAXN,                       */
    /*                 NEV <= MAXNEV,                     */
    /*             NEV + 1 <= NCV <= MAXNCV               */
    /* -------------------------------------------------- */

    int n = 100;
    int nev = 4;
    int ncv = 10;
    if (n > MAXN)
    {
        printf(" ERROR with _SDRV4: N is greater than MAXN \n");
        return 0;
    }
    else if (nev > MAXNEV)
    {
        printf(" ERROR with _SDRV4: NEV is greater than MAXNEV \n");
        return 0;
    }
    else if (ncv > MAXNCV)
    {
        printf(" ERROR with _SDRV4: NCV is greater than MAXNCV \n");
        return 0;
    }
    char* bmat = "G";
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

    float* resid = (float*)malloc(n * sizeof(float));
    float* v = (float*)malloc(n * ncv * sizeof(float));
    float* workl = (float*)malloc(lworkl * sizeof(float));
    float* workd = (float*)malloc(3 * n * sizeof(float));

    /* ------------------------------------------------- */
    /* This program uses exact shifts with respect to    */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 3 specified in the      */
    /* documentation of SSAUPD is used (IPARAM(7) = 3).  */
    /* All these options may be changed by the user.     */
    /* For details, see the documentation in SSAUPD.     */
    /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode = 3;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

    /* ----------------------------------------------------- */
    /* Call LAPACK routine to factor the tridiagonal matrix  */
    /* (A-SIGMA*M).  The matrix A is the 1-d discrete        */
    /* Laplacian. The matrix M is the associated mass matrix */
    /* arising from using piecewise linear finite elements   */
    /* on the interval [0, 1].                               */
    /* ----------------------------------------------------- */

    h = 1.0f / (float) (n + 1);
    r1 = h * .66666666666666663f;
    r2 = h * .16666666666666666f;
    for (j = 1; j <= n; ++j)
    {
        ad[j - 1] = 2.0f / h - sigma * r1;
        adl[j - 1] = -1.0f / h - sigma * r2;
    }
    scopy_(&n, adl, &i_one, adu, &i_one);
    sgttrf_(&n, adl, ad, adu, adu2, ipiv, &ierr);
    if (ierr != 0)
    {
        printf(" Error with _gttrf in _SDRV4.\n");
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

    if (ido == -1)
    {
        /* ------------------------------------------ */
        /* Perform  y <--- OP*x = inv[A-SIGMA*M]*M*x  */
        /* to force the starting vector into the      */
        /* range of OP.  The user should supply       */
        /* his/her own matrix vector multiplication   */
        /* routine and a linear system solver here.   */
        /* The matrix vector multiplication routine   */
        /* takes workd(ipntr(1)) as the input vector. */
        /* The final result is returned to            */
        /* workd(ipntr(2)).                           */
        /* ------------------------------------------ */

        ssdrv4_mv_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        sgttrs_("N", &n, &i_one, adl, ad, adu, adu2, ipiv, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" Error with _gttrs in _SDRV4. \n");
            printf(" \n");
            goto EXIT;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call SSAUPD again. */
        /* --------------------------------------- */

        goto L10;
    }
    else if (ido == 1)
    {
        /* --------------------------------------- */
        /* Perform y <-- OP*x = inv[A-sigma*M]*M*x */
        /* M*x has been saved in workd(ipntr(3)).  */
        /* the user only needs the linear system   */
        /* solver here that takes workd(ipntr(3)   */
        /* as input, and returns the result to     */
        /* workd(ipntr(2)).                        */
        /* --------------------------------------- */

        scopy_(&n, &workd[ipntr[2] - 1], &i_one, &workd[ipntr[1] - 1], &i_one);
        sgttrs_("N", &n, &i_one, adl, ad, adu, adu2, ipiv, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" Error with _gttrs in _SDRV4.\n");
            printf(" \n");
            goto EXIT;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call SSAUPD again. */
        /* --------------------------------------- */

        goto L10;
    }
    else if (ido == 2)
    {
        /* --------------------------------------- */
        /*          Perform  y <--- M*x            */
        /* Need the matrix vector multiplication   */
        /* routine here that takes workd(ipntr(1)) */
        /* as the input and returns the result to  */
        /* workd(ipntr(2)).                        */
        /* --------------------------------------- */

        ssdrv4_mv_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call SSAUPD again. */
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
        /* documentation in SSAUPD. */
        /* ------------------------ */

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

        ssdrv4_av_(n, &v[k], workd);
        ssdrv4_mv_(n, &v[k], &workd[n]);
        r__1 = -d[j - 1];
        saxpy_(&n, &r__1, &workd[n], &i_one, workd, &i_one);
        d[j + 24] = snrm2_(&n, workd, &i_one);
        d[j + 24] /= (r__1 = d[j - 1], dabs(r__1));

    }

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
    printf(" _SDRV4 \n");
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

    free(resid);
    free(v);
    free(workl);
    free(workd);

    /* ------------------------- */
    /* Done with program ssdrv4. */
    /* ------------------------- */

    return ierr;
}

/* ------------------------------------------------------------------------ */
/*     matrix vector subroutine */
/*     The matrix used is the 1 dimensional mass matrix */
/*     on the interval [0,1]. */

int ssdrv4_mv_(const int n, float *v, float *w)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    float h;
    int j;

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    w[1] = v[1] * 4.0f + v[2];
    i__1 = n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = v[j - 1] + v[j] * 4.0f + v[j + 1];
    }
    j = n;
    w[j] = v[j - 1] + v[j] * 4.0f;

    /*     Scale the vector w by h. */

    h = 1.0f / ((float) (n + 1) * 6.0f);
    sscal_(&n, &h, &w[1], &i_one);
    return 0;
} /* mv_ */

/* ------------------------------------------------------------------------ */
/*     matrix vector subroutine */
/*     where the matrix is the finite element discretization of the */
/*     1 dimensional discrete Laplacian on [0,1] with zero Dirichlet */
/*     boundary condition using piecewise linear elements. */

int ssdrv4_av_(const int n, float *v, float *w)
{
    /* System generated locals */
    int i__1;
    float r__1;

    /* Local variables */
    float h;
    int j;

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

    /*     Scale the vector w by (1/h) */

    h = 1.0f / (float) (n + 1);
    r__1 = 1.0f / h;
    sscal_(&n, &r__1, &w[1], &i_one);
    return 0;
} /* av_ */

