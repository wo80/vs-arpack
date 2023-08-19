/* EXAMPLES\COMPLEX\cndrv2.f -- translated by f2c (version 20100827). */

#include <stdlib.h>
#include "arpack_internal.h"

#define RHO 10.0f

/**
 * \BeginDoc
 *
 *     Simple program to illustrate the idea of reverse communication
 *     in shift-invert mode for a standard complex nonsymmetric eigenvalue
 *     problem.
 *
 *     We implement example two of ex-complex.doc in DOCUMENTS directory
 *
 * \Example-2
 *     ... Suppose we want to solve A*x = lambda*x in shift-invert mode,
 *         where A is derived from the central difference discretization
 *         of the 1-dimensional convection-diffusion operator
 *                   (d^2u/dx^2) + rho*(du/dx)
 *         on the interval [0,1] with zero Dirichlet boundary condition.
 *     ... The shift sigma is a complex number.
 *
 *     ... OP = inv[A-sigma*I] and  B = I.
 *
 *     ... Use mode 3 of CNAUPD.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Routines called:
 *     cnaupd  ARPACK reverse communication interface routine.
 *     cneupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     cgttrf  LAPACK tridiagonal factorization routine.
 *     cgttrs  LAPACK tridiagonal solve routine.
 *     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     caxpy   Level 1 BLAS that computes y <- alpha*x+y.
 *     ccopy   Level 1 BLAS that copies one vector to another.
 *     scnrm2  Level 1 BLAS that computes the norm of a vector.
 *     av      Matrix vector multiplication routine that computes A*x.
 *
 * \EndLib
 */
int cndrv2()
{
    /* System generated locals */
    int i__1, i__2;
    complex q__1;

    /* Local variables */
    complex d[25];
    complex workev[50];
    float rd[75] /* (3 * MAXNCV) */;
    float rwork[256];

    float h, s, h2, s1, s3;
    complex sigma, s2;

    int j;
    int ierr, nconv;
    int ishfts, maxitr, mode;
    int ipiv[256];
    int ipntr[14];
    int iparam[11];
    bool select[25];
    bool rvec;

    /* Define maximum dimensions for all arrays. */

    const int MAXN   = 256; /* Maximum dimension of the A allowed. */
    const int MAXNEV =  10; /* Maximum NEV allowed */
    const int MAXNCV =  25; /* Maximum NCV allowed */

    /* ------------------------------------------------ */
    /* The number N is the dimension of the matrix.  A  */
    /* standard eigenvalue problem is solved (BMAT =    */
    /* 'I').  NEV is the number of eigenvalues (closest */
    /* to the shift SIGMA) to be approximated.  Since   */
    /* the shift-invert mode is used, WHICH is set to   */
    /* 'LM'.  The user can modify NEV, NCV, SIGMA to    */
    /* solve problems of different sizes, and to get    */
    /* different parts of the spectrum.  However, The   */
    /* following conditions must be satisfied:          */
    /*                 N <= MAXN,                       */
    /*               NEV <= MAXNEV,                     */
    /*           NEV + 2 <= NCV <= MAXNCV               */
    /* ------------------------------------------------ */

    int n = 100;
    int nev = 4;
    int ncv = 20;
    if (n > 256)
    {
        printf(" ERROR with _NDRV2: N is greater than MAXN \n");
        return 0;
    }
    else if (nev > 10)
    {
        printf(" ERROR with _NDRV2: NEV is greater than MAXNEV \n");
        return 0;
    }
    else if (ncv > 25)
    {
        printf(" ERROR with _NDRV2: NCV is greater than MAXNCV \n");
        return 0;
    }
    char* bmat = "I";
    char* which = "LM";
    sigma.r = 0.0f, sigma.i = 0.0f;

    /* -------------------------------------------------- */
    /* Construct C = A - SIGMA*I, factor C in complex     */
    /* arithmetic (using LAPACK subroutine cgttrf). The   */
    /* matrix A is chosen to be the tridiagonal matrix    */
    /* derived from standard central difference of the    */
    /* 1-d convection diffusion operator - u``+ rho*u` on */
    /* the interval [0, 1] with zero Dirichlet boundary   */
    /* condition.                                         */
    /* -------------------------------------------------- */

    complex* du = (complex*)malloc(n * sizeof(complex));
    complex* dd = (complex*)malloc(n * sizeof(complex));
    complex* dl = (complex*)malloc(n * sizeof(complex));
    complex* du2 = (complex*)malloc(n * sizeof(complex));

    h = 1.0f / (float)(n + 1);
    h2 = h * h;
    s = RHO / 2.0f;

    s1 = -1.0f / h2 - s / h;
    s2.r = 2.0f / h2 - sigma.r, s2.i = -sigma.i;
    s3 = -1.0f / h2 + s / h;

    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        i__2 = j - 1;
        dl[i__2].r = s1, dl[i__2].i = 0.0f;
        dd[i__2].r = s2.r, dd[i__2].i = s2.i;
        du[i__2].r = s3, du[i__2].i = 0.0f;
    }
    i__1 = n - 1;
    dd[i__1].r = s2.r, dd[i__1].i = s2.i;

    cgttrf_(&n, dl, dd, du, du2, ipiv, &ierr);
    if (ierr != 0)
    {
        printf(" \n");
        printf(" ERROR with _gttrf in _NDRV2.\n");
        printf(" \n");
        return 0;
    }

    /* --------------------------------------------------- */
    /* The work array WORKL is used in CNAUPD as           */
    /* workspace.  Its dimension LWORKL is set as          */
    /* illustrated below.  The parameter TOL determines    */
    /* the stopping criterion. If TOL<=0, machine          */
    /* precision is used.  The variable IDO is used for    */
    /* reverse communication, and is initially set to 0.   */
    /* Setting INFO=0 indicates that a random vector is    */
    /* generated in CNAUPD to start the Arnoldi iteration. */
    /* --------------------------------------------------- */

    /* Computing 2nd power */
    int lworkl = ncv * ncv * 3 + ncv * 5;
    float tol = 0.0f;
    int ido = 0;
    int info = 0;

    complex* ax = (complex*)malloc(n * sizeof(complex));
    complex* resid = (complex*)malloc(n * sizeof(complex));
    complex* v = (complex*)malloc(n * ncv * sizeof(complex));
    complex* workl = (complex*)malloc(lworkl * sizeof(complex));
    complex* workd = (complex*)malloc(3 * n * sizeof(complex));

    /* ------------------------------------------------- */
    /* This program uses exact shifts with respect to    */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed. Mode 3 of CNAUPD is used      */
    /* (IPARAM(7) = 3).  All these options can be        */
    /* changed by the user. For details see the          */
    /* documentation in CNAUPD.                          */
    /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode = 3;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

    /* ----------------------------------------- */
    /* M A I N   L O O P (Reverse communication) */
    /* ----------------------------------------- */

L20:

    /* ------------------------------------------- */
    /* Repeatedly call the routine CNAUPD and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    cnaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, rwork, &info);

    if (ido == -1 || ido == 1)
    {
        /* ----------------------------------------- */
        /* Perform  y <--- OP*x = inv[A-SIGMA*I]*x   */
        /* The user should supply his/her own linear */
        /* system solver here that takes             */
        /* workd(ipntr(1)) as the input, and returns */
        /* the result to workd(ipntr(2)).            */
        /* ----------------------------------------- */

        ccopy_(&n, &workd[ipntr[0] - 1], &c__1, &workd[ipntr[1] - 1], &c__1);

        cgttrs_("N", &n, &c__1, dl, dd, du, du2, ipiv, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" ERROR with _gttrs in _NDRV2.\n");
            printf(" \n");
            goto EXIT;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call CNAUPD again. */
        /* --------------------------------------- */

        goto L20;
    }

    /* --------------------------------------- */
    /* Either we have convergence, or there is */
    /* an error.                               */
    /* --------------------------------------- */

    if (info < 0)
    {
        /* ------------------------ */
        /* Error message, check the */
        /* documentation in CNAUPD  */
        /* ------------------------ */

        printf(" \n");
        printf(" Error with _naupd info = %d\n", info);
        printf(" Check the documentation of _naupd.\n");
        printf(" \n");

        ierr = info;
        goto EXIT;
    }

    /* ----------------------------------------- */
    /* No fatal errors occurred.                 */
    /* Post-Process using CNEUPD.                */
    /*                                           */
    /* Computed eigenvalues may be extracted.    */
    /*                                           */
    /* Eigenvectors may also be computed now if  */
    /* desired.  (indicated by rvec = .true.)    */
    /* ----------------------------------------- */

    rvec = true;

    cneupd_(&rvec, "A", select, d, v, &n, &sigma, workev, bmat, &n,which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);

    /* -------------------------------------------- */
    /* Eigenvalues are returned in the one          */
    /* dimensional array D.  The corresponding      */
    /* eigenvectors are returned in the first NCONV */
    /* (=IPARAM(5)) columns of the two dimensional  */
    /* array V if requested.  Otherwise, an         */
    /* orthogonal basis for the invariant subspace  */
    /* corresponding to the eigenvalues in D is     */
    /* returned in V.                               */
    /* -------------------------------------------- */

    if (ierr != 0)
    {
        /* ---------------------------------- */
        /* Error condition:                   */
        /* Check the documentation of CNEUPD. */
        /* ---------------------------------- */

        printf(" \n");
        printf(" Error with _neupd info = %d\n", ierr);
        printf(" Check the documentation of _neupd. \n");
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

        cndrv2_av_(n, &v[k], ax);
        i__2 = j - 1;
        q__1.r = -d[i__2].r, q__1.i = -d[i__2].i;
        caxpy_(&n, &q__1, &v[k], &c__1, ax, &c__1);
        i__2 = j - 1;
        rd[j - 1] = d[i__2].r;
        rd[j + 24] = d[i__2].i;
        rd[j + 49] = scnrm2_(&n, ax, &c__1);
        rd[j + 49] /= slapy2_(&rd[j - 1], &rd[j + 24]);

    }

    /* --------------------------- */
    /* Display computed residuals. */
    /* --------------------------- */

    smout_(&nconv, &c__3, rd, &c__25, &c_n6, "Ritz values (Real, Imag) and relative residuals");

    /* ----------------------------------------- */
    /* Print additional convergence information. */
    /* ----------------------------------------- */

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
    printf("_NDRV2 \n");
    printf("====== \n");
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

    free(du);
    free(dd);
    free(dl);
    free(du2);

    free(ax);
    free(resid);
    free(v);
    free(workl);
    free(workd);

    /* ------------------------- */
    /* Done with program cndrv2. */
    /* ------------------------- */

    return ierr;
}

/**
 * Matrix vector multiplication subroutine.
 */
int cndrv2_av_(const int n, complex *v, complex *w)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    int j;
    float h, s, h2, dd, dl, du;

    /* Function Body */
    h = 1.0f / (float) (n + 1);
    h2 = h * h;
    s = RHO / 2.0f;
    dd = 2.0f / h2;
    dl = -1.0f / h2 - s / h;
    du = -1.0f / h2 + s / h;

    w[0].r = dd * v[0].r + du * v[1].r;
    w[0].i = dd * v[0].i + du * v[1].i;

    i__1 = n - 1;
    for (j = 1; j < i__1; ++j)
    {
        i__2 = j - 1;
        i__3 = j + 1;
        w[j].r = (dl * v[i__2].r + dd * v[j].r) + du * v[i__3].r;
        w[j].i = (dl * v[i__2].i + dd * v[j].i) + du * v[i__3].i;
    }
    i__1 = n - 1;
    i__2 = n - 2;
    w[i__1].r = dl * v[i__2].r + dd * v[i__1].r;
    w[i__1].i = dl * v[i__2].i + dd * v[i__1].i;
    return 0;
} /* av_ */
