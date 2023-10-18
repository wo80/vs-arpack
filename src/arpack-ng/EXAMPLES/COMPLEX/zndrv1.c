/* EXAMPLES\COMPLEX\zndrv1.f -- translated by f2c (version 20100827). */

#include <stdlib.h>
#include <stdio.h>
#include "arpack.h"
#include "lapack.h"

int zndrv1_av_(const int nx, a_dcomplex* v, a_dcomplex* w);
int zndrv1_tv_(const int nx, a_dcomplex* x, a_dcomplex* y);

extern int dmout_(const int, const int, const double*, const int, const int, const char*);

static int i_one = 1;

/**
 * \BeginDoc
 *
 *     Example program to illustrate the idea of reverse communication
 *     for a standard complex nonsymmetric eigenvalue problem.
 *
 *     We implement example one of ex-complex.doc in DOCUMENTS directory
 *
 * \Example-1
 *     ... Suppose we want to solve A*x = lambda*x in regular mode,
 *         where A is obtained from the standard central difference
 *         discretization of the convection-diffusion operator
 *                 (Laplacian u) + rho*(du / dx)
 *         on the unit squre [0,1]x[0,1] with zero Dirichlet boundary
 *         condition.
 *
 *     ... OP = A  and  B = I.
 *
 *     ... Assume "call av (nx,x,y)" computes y = A*x
 *
 *     ... Use mode 1 of ZNAUPD .
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Routines called
 *     znaupd   ARPACK reverse communication interface routine.
 *     zneupd   ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     dznrm2   Level 1 BLAS that computes the norm of a complex vector.
 *     zaxpy    Level 1 BLAS that computes y <- alpha*x+y.
 *     av      Matrix vector multiplication routine that computes A*x.
 *     tv      Matrix vector multiplication routine that computes T*x,
 *             where T is a tridiagonal matrix.  It is used in routine
 *             av.
 *
 * \EndLib
 */
int main()
{
    /* System generated locals */
    int i__1, i__2;
    a_dcomplex z__1;

    /* Local variables */
    a_dcomplex workev[90];
    a_dcomplex d[30];
    double rd[90] /* (3 * MAXNCV) */;
    double rwork[30];

    a_dcomplex sigma;

    int j;
    int ierr, nconv;
    int ishfts, maxitr, mode;
    int ipntr[14];
    int iparam[11];
    logical select[30];
    logical rvec;

    /* Define maximum dimensions for all arrays. */

    const int MAXN   = 256; /* Maximum dimension of the A allowed. */
    const int MAXNEV =  12; /* Maximum NEV allowed */
    const int MAXNCV =  30; /* Maximum NCV allowed */

    /* ------------------------------------------------ */
    /* The number NX is the number of interior points   */
    /* in the discretization of the 2-dimensional       */
    /* convection-diffusion operator on the unit        */
    /* square with zero Dirichlet boundary condition.   */
    /* The number N(=NX*NX) is the dimension of the     */
    /* matrix.  A standard eigenvalue problem is        */
    /* solved (BMAT = 'I').  NEV is the number of       */
    /* eigenvalues to be approximated.  The user can    */
    /* modify NX, NEV, NCV, WHICH to solve problems of  */
    /* different sizes, and to get different parts of   */
    /* the spectrum.  However, The following            */
    /* conditions must be satisfied:                    */
    /*                   N <= MAXN                      */
    /*                 NEV <= MAXNEV                    */
    /*           NEV + 2 <= NCV <= MAXNCV               */
    /* ------------------------------------------------ */

    int nx = 10;
    int n = nx * nx;
    int nev = 4;
    int ncv = 20;
    if (n > MAXN)
    {
        printf(" ERROR with _NDRV1: N is greater than MAXN \n");
        return 0;
    }
    else if (nev > MAXNEV)
    {
        printf(" ERROR with _NDRV1: NEV is greater than MAXNEV \n");
        return 0;
    }
    else if (nev > MAXNCV)
    {
        printf(" ERROR with _NDRV1: NCV is greater than MAXNCV \n");
        return 0;
    }
    char* bmat = "I";
    char* which = "LM";

    /* ------------------------------------------------- */
    /* The work array WORKL is used in ZNAUPD  as         */
    /* workspace.  Its dimension LWORKL is set as        */
    /* illustrated below.  The parameter TOL determines  */
    /* the stopping criterion. If TOL<=0, machine        */
    /* precision is used.  The variable IDO is used for  */
    /* reverse communication, and is initially set to 0. */
    /* Setting INFO=0 indicates that a random vector is  */
    /* generated to start the ARNOLDI iteration.         */
    /* ------------------------------------------------- */

    /* Computing 2nd power */
    int lworkl = ncv * ncv * 3 + ncv * 5;
    double tol = 0.0;
    int ido = 0;
    int info = 0;

    a_dcomplex* ax = (a_dcomplex*)malloc(n * sizeof(a_dcomplex));
    a_dcomplex* resid = (a_dcomplex*)malloc(n * sizeof(a_dcomplex));
    a_dcomplex* v = (a_dcomplex*)malloc(n * ncv * sizeof(a_dcomplex));
    a_dcomplex* workl = (a_dcomplex*)malloc(lworkl * sizeof(a_dcomplex));
    a_dcomplex* workd = (a_dcomplex*)malloc(3 * n * sizeof(a_dcomplex));

    /* ------------------------------------------------- */
    /* This program uses exact shift with respect to     */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 1 of ZNAUPD  is used     */
    /* (IPARAM(7) = 1). All these options can be changed */
    /* by the user. For details see the documentation in */
    /* ZNAUPD .                                           */
    /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode = 1;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

    /* ----------------------------------------- */
    /* M A I N   L O O P (Reverse communication) */
    /* ----------------------------------------- */

L10:

    /* ------------------------------------------- */
    /* Repeatedly call the routine ZNAUPD  and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    znaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, rwork, &info);

    if (ido == -1 || ido == 1)
    {
        /* ----------------------------------------- */
        /* Perform matrix vector multiplication      */
        /*                y <--- OP*x                */
        /* The user should supply his/her own        */
        /* matrix vector multiplication routine here */
        /* that takes workd(ipntr(1)) as the input   */
        /* vector, and return the matrix vector      */
        /* product to workd(ipntr(2)).               */
        /* ----------------------------------------- */

        zndrv1_av_(nx, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call ZNAUPD  again. */
        /* --------------------------------------- */

        goto L10;
    }

    /* -------------------------------------- */
    /* Either we have convergence or there is */
    /* an error.                              */
    /* -------------------------------------- */

    if (info < 0)
    {
        /* ------------------------ */
        /* Error message, check the */
        /* documentation in ZNAUPD   */
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
    /* Post-Process using ZNEUPD .                */
    /*                                           */
    /* Computed eigenvalues may be extracted.    */
    /*                                           */
    /* Eigenvectors may also be computed now if  */
    /* desired.  (indicated by rvec = .true.)    */
    /* ----------------------------------------- */

    rvec = TRUE_;

    zneupd_(&rvec, "A", select, d, v, &n, &sigma, workev, bmat, &n,which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);

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
        /* Check the documentation of ZNEUPD . */
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

        zndrv1_av_(nx, &v[k], ax);
        i__2 = j - 1;
        z__1.r = -d[i__2].r, z__1.i = -d[i__2].i;
        zaxpy_(&n, &z__1, &v[k], &i_one, ax, &i_one);
        i__2 = j - 1;
        rd[j - 1] = d[i__2].r;
        rd[j + 29] = d[i__2].i;
        rd[j + 59] = dznrm2_(&n, ax, &i_one);
        rd[j + 59] /= dlapy2_(&rd[j - 1], &rd[j + 29]);

    }

    /* --------------------------- */
    /* Display computed residuals. */
    /* --------------------------- */

    dmout_(nconv, 3, rd, MAXNCV, -6, "Ritz values (Real, Imag) and relative residuals");

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
    printf("_NDRV1\n");
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

    free(ax);
    free(resid);
    free(v);
    free(workl);
    free(workd);

    /* ------------------------- */
    /* Done with program zndrv1 . */
    /* ------------------------- */

    return ierr;
}

/** Matrix vector subroutine.
 *
 * The matrix used is the convection-diffusion operator
 * discretized using centered difference.
 * 
 * Computes w <--- OP*v, where OP is the nx*nx by nx*nx block
 * tridiagonal matrix
 *
 *              | T -I          |
 *              |-I  T -I       |
 *         OP = |   -I  T       |
 *              |        ...  -I|
 *              |           -I T|
 *
 * derived from the standard central difference  discretization
 * of the 2-dimensional convection-diffusion operator
 *              (Laplacian u) + rho*(du/dx)
 * on the unit squqre with zero boundary condition.
 *
 * The subroutine TV is called to computed y<---T*x.
 */
int zndrv1_av_(const int nx, a_dcomplex *v, a_dcomplex *w)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    int j;
    int lo;
    double h2;
    a_dcomplex z;

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    h2 = 1.0 / (double)((nx + 1) * (nx + 1));

    z.r = -1.0 / h2;
    z.i = -0.0;

    zndrv1_tv_(nx, &v[1], &w[1]);
    zaxpy_(&nx, &z, &v[nx + 1], &i_one, &w[1], &i_one);

    i__1 = nx - 1;
    for (j = 2; j <= i__1; ++j)
    {
        lo = (j - 1) * nx;
        zndrv1_tv_(nx, &v[lo + 1], &w[lo + 1]);
        zaxpy_(&nx, &z, &v[lo - nx + 1], &i_one, &w[lo + 1], &i_one);
        zaxpy_(&nx, &z, &v[lo + nx + 1], &i_one, &w[lo + 1], &i_one);
    }

    lo = (nx - 1) * nx;
    zndrv1_tv_(nx, &v[lo + 1], &w[lo + 1]);
    zaxpy_(&nx, &z, &v[lo - nx + 1], &i_one, &w[lo + 1], &i_one);

    return 0;
} /* av_ */

/**
 * Compute the matrix vector multiplication y<---T*x
 * where T is a nx by nx tridiagonal matrix with DD on the
 * diagonal, DL on the subdiagonal, and DU on the superdiagonal
 */
int zndrv1_tv_(const int nx, a_dcomplex *x, a_dcomplex *y)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    int j;
    double h, h2, dd, dl, du;

    /* Function Body */
    h = 1.0 / (double)(nx + 1);
    h2 = h * h;
    dd = 4.0 / h2;
    dl = -1.0 / h2 - 50.0 / h;
    du = -1.0 / h2 + 50.0 / h;

    y[0].r = dd * x[0].r + du * x[1].r;
    y[0].i = dd * x[0].i + du * x[1].i;

    i__1 = nx - 1;
    for (j = 1; j < i__1; ++j)
    {
        i__2 = j - 1;
        i__3 = j + 1;
        y[j].r = (dl * x[i__2].r + dd * x[j].r) + du * x[i__3].r;
        y[j].i = (dl * x[i__2].i + dd * x[j].i) + du * x[i__3].i;
    }
    i__1 = nx - 1;
    i__2 = nx - 2;
    y[i__1].r = dl * x[i__2].r + dd * x[i__1].r;
    y[i__1].i = dl * x[i__2].i + dd * x[i__1].i;
    return 0;
} /* tv_ */
