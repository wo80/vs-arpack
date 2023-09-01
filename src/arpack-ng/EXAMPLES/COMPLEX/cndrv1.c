/* EXAMPLES\COMPLEX\cndrv1.f -- translated by f2c (version 20100827). */

#include <stdlib.h>
#include "arpack.h"

int cndrv1_av_(const int nx, a_fcomplex* v, a_fcomplex* w);
int cndrv1_tv_(const int nx, a_fcomplex* x, a_fcomplex* y);

extern int smout_(const int, const int, const float*, const int, const int, const char*);

static int c__1 = 1;

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
 *     ... Use mode 1 of CNAUPD.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Routines called
 *     cnaupd  ARPACK reverse communication interface routine.
 *     cneupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     scnrm2  Level 1 BLAS that computes the norm of a complex vector.
 *     caxpy   Level 1 BLAS that computes y <- alpha*x+y.
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
    a_fcomplex q__1;

    /* Local variables */
    a_fcomplex workev[90];
    a_fcomplex d[30];
    float rd[90] /* (3 * MAXNCV) */;
    float rwork[30];

    a_fcomplex sigma;

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
    if (n > 256)
    {
        printf(" ERROR with _NDRV1: N is greater than MAXN \n");
        return 0;
    }
    else if (nev > 12)
    {
        printf(" ERROR with _NDRV1: NEV is greater than MAXNEV \n");
        return 0;
    }
    else if (ncv > 30)
    {
        printf(" ERROR with _NDRV1: NCV is greater than MAXNCV \n");
        return 0;
    }
    char* bmat = "I";
    char* which = "LM";

    /* ------------------------------------------------- */
    /* The work array WORKL is used in CNAUPD as         */
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
    float tol = 0.0f;
    int ido = 0;
    int info = 0;

    a_fcomplex* ax = (a_fcomplex*)malloc(n * sizeof(a_fcomplex));
    a_fcomplex* resid = (a_fcomplex*)malloc(n * sizeof(a_fcomplex));
    a_fcomplex* v = (a_fcomplex*)malloc(n * ncv * sizeof(a_fcomplex));
    a_fcomplex* workl = (a_fcomplex*)malloc(lworkl * sizeof(a_fcomplex));
    a_fcomplex* workd = (a_fcomplex*)malloc(3 * n * sizeof(a_fcomplex));

    /* ------------------------------------------------- */
    /* This program uses exact shift with respect to     */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 1 of CNAUPD is used     */
    /* (IPARAM(7) = 1). All these options can be changed */
    /* by the user. For details see the documentation in */
    /* CNAUPD.                                           */
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
    /* Repeatedly call the routine CNAUPD and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    cnaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, rwork, &info);

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

        cndrv1_av_(nx, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call CNAUPD again. */
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

    rvec = TRUE_;

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

        cndrv1_av_(nx, &v[k], ax);
        i__2 = j - 1;
        q__1.r = -d[i__2].r, q__1.i = -d[i__2].i;
        caxpy_(&n, &q__1, &v[k], &c__1, ax, &c__1);
        i__2 = j - 1;
        rd[j - 1] = d[i__2].r;
        rd[j + 29] = d[i__2].i;
        rd[j + 59] = scnrm2_(&n, ax, &c__1);
        rd[j + 59] /= slapy2_(&rd[j - 1], &rd[j + 29]);

    }

    /* --------------------------- */
    /* Display computed residuals. */
    /* --------------------------- */

    smout_(nconv, 3, rd, 30, -6, "Ritz values (Real, Imag) and relative residuals");

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
    /* Done with program cndrv1. */
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
int cndrv1_av_(const int nx, a_fcomplex *v, a_fcomplex *w)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    int j;
    int lo;
    float h2;
    a_fcomplex z;

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    h2 = 1.0f / (float)((nx + 1) * (nx + 1));

    z.r = -1.0f / h2;
    z.i = -0.0f;

    cndrv1_tv_(nx, &v[1], &w[1]);
    caxpy_(&nx, &z, &v[nx + 1], &c__1, &w[1], &c__1);

    i__1 = nx - 1;
    for (j = 2; j <= i__1; ++j)
    {
        lo = (j - 1) * nx;
        cndrv1_tv_(nx, &v[lo + 1], &w[lo + 1]);
        caxpy_(&nx, &z, &v[lo - nx + 1], &c__1, &w[lo + 1], &c__1);
        caxpy_(&nx, &z, &v[lo + nx + 1], &c__1, &w[lo + 1], &c__1);
    }

    lo = (nx - 1) * nx;
    cndrv1_tv_(nx, &v[lo + 1], &w[lo + 1]);
    caxpy_(&nx, &z, &v[lo - nx + 1], &c__1, &w[lo + 1], &c__1);

    return 0;
} /* av_ */

/**
 * Compute the matrix vector multiplication y<---T*x
 * where T is a nx by nx tridiagonal matrix with DD on the
 * diagonal, DL on the subdiagonal, and DU on the superdiagonal
 */
int cndrv1_tv_(const int nx, a_fcomplex *x, a_fcomplex *y)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    int j;
    float h, h2, dd, dl, du;

    /* Function Body */
    h = 1.0f / (float)(nx + 1);
    h2 = h * h;
    dd = 4.0f / h2;
    dl = -1.0f / h2 - 50.0f / h;
    du = -1.0f / h2 + 50.0f / h;

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
