/* EXAMPLES\SYM\ssdrv1.f -- translated by f2c (version 20100827). */

#include <stdlib.h>
#include "arpack.h"

int ssdrv1_av_(const int nx, float* v, float* w);

extern int smout_(const int, const int, const float*, const int, const int, const char*);

static int c__1 = 1;
static float s_m1 = -1.0f;

/**
 * \BeginDoc
 *
 *     Simple program to illustrate the idea of reverse communication
 *     in regular mode for a standard symmetric eigenvalue problem.
 *
 *     We implement example one of ex-sym.doc in SRC directory
 *
 * \Example-1
 *     ... Suppose we want to solve A*x = lambda*x in regular mode,
 *         where A is derived from the central difference discretization
 *         of the 2-dimensional Laplacian on the unit square [0,1]x[0,1]
 *         with zero Dirichlet boundary condition.
 *
 *     ... OP = A  and  B = I.
 *
 *     ... Assume "call av (n,x,y)" computes y = A*x.
 *
 *     ... Use mode 1 of SSAUPD.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Routines called:
 *     ssaupd  ARPACK reverse communication interface routine.
 *     sseupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     snrm2   Level 1 BLAS that computes the norm of a vector.
 *     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
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
    float r__1;

    /* Local variables */
    float d[50] /* (2 * MAXNCV) */;
    float sigma;

    int j;
    int ierr, nconv;
    int ishfts, maxitr, mode;
    int ipntr[11];
    int iparam[11];
    logical select[25];
    logical rvec;

    /* Define maximum dimensions for all arrays. */

    const int MAXN   = 256; /* Maximum dimension of the A allowed. */
    const int MAXNEV =  10; /* Maximum NEV allowed */
    const int MAXNCV =  25; /* Maximum NCV allowed */

    /* -------------------------------------------------- */
    /* The number NX is the number of interior points     */
    /* in the discretization of the 2-dimensional         */
    /* Laplacian on the unit square with zero Dirichlet   */
    /* boundary condition.  The number N(=NX*NX) is the   */
    /* dimension of the matrix.  A standard eigenvalue    */
    /* problem is solved (BMAT = 'I'). NEV is the number  */
    /* of eigenvalues to be approximated.  The user can   */
    /* modify NEV, NCV, WHICH to solve problems of        */
    /* different sizes, and to get different parts of the */
    /* spectrum.  However, The following conditions must  */
    /* be satisfied:                                      */
    /*                   N <= MAXN,                       */
    /*                 NEV <= MAXNEV,                     */
    /*             NEV + 1 <= NCV <= MAXNCV               */
    /* -------------------------------------------------- */

    int nx = 10;
    int n = nx * nx;
    int nev = 4;
    int ncv = 10;
    if (n > 256)
    {
        printf(" ERROR with _SDRV1: N is greater than MAXN \n");
        return 0;
    }
    else if (nev > 10)
    {
        printf(" ERROR with _SDRV1: NEV is greater than MAXNEV \n");
        return 0;
    }
    else if (ncv > 25)
    {
        printf(" ERROR with _SDRV1: NCV is greater than MAXNCV \n");
        return 0;
    }
    char* bmat = "I";
    char* which = "SM";

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
    int info = 0;
    int ido = 0;

    float* ax = (float*)malloc(n * sizeof(float));
    float* resid = (float*)malloc(n * sizeof(float));
    float* v = (float*)malloc(n * ncv * sizeof(float));
    float* workl = (float*)malloc(lworkl * sizeof(float));
    float* workd = (float*)malloc(3 * n * sizeof(float));

    /* ------------------------------------------------- */
    /* This program uses exact shifts with respect to    */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 1 of SSAUPD is used     */
    /* (IPARAM(7) = 1).  All these options may be        */
    /* changed by the user. For details, see the         */
    /* documentation in SSAUPD.                          */
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
    /* Repeatedly call the routine SSAUPD and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    ssaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);

    if (ido == -1 || ido == 1)
    {
        /* ------------------------------------ */
        /* Perform matrix vector multiplication */
        /*              y <--- OP*x             */
        /* The user should supply his/her own   */
        /* matrix vector multiplication routine */
        /* here that takes workd(ipntr(1)) as   */
        /* the input, and return the result to  */
        /* workd(ipntr(2)).                     */
        /* ------------------------------------ */

        ssdrv1_av_(nx, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

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
        /* ------------------------ */
        /* Error message. Check the */
        /* documentation in SSAUPD. */
        /* ------------------------ */

        printf(" \n");
        printf(" Error with _saupd info = %d\n", info);
        printf(" Check the documentation of _saupd.\n");
        printf(" \n");

        ierr = info;
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
        printf(" Check the documentation of _seupd. \n");
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

        ssdrv1_av_(nx, &v[k], ax);
        r__1 = -d[j - 1];
        saxpy_(&n, &r__1, &v[k], &c__1, ax, &c__1);
        d[j + 24] = snrm2_(&n, ax, &c__1);
        d[j + 24] /= (r__1 = d[j - 1], dabs(r__1));

    }

    /* ----------------------------- */
    /* Display computed residuals    */
    /* ----------------------------- */

    smout_(nconv, 2, d, 25, -6, "Ritz values and relative residuals");
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
    printf(" _SDRV1 \n");
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

    free(ax);
    free(resid);
    free(v);
    free(workl);
    free(workd);

    /* ------------------------- */
    /* Done with program ssdrv1. */
    /* ------------------------- */

    return ierr;
}

/* ------------------------------------------------------------------ */
/*     matrix vector subroutine */

/*     The matrix used is the 2 dimensional discrete Laplacian on unit */
/*     square with zero Dirichlet boundary condition. */

/*     Computes w <--- OP*v, where OP is the nx*nx by nx*nx block */
/*     tridiagonal matrix */

/*                  | T -I          | */
/*                  |-I  T -I       | */
/*             OP = |   -I  T       | */
/*                  |        ...  -I| */
/*                  |           -I T| */

/*     The subroutine TV is called to computed y<---T*x. */

int ssdrv1_av_(const int nx, float *v, float *w)
{
    /* System generated locals */
    int i__1;
    float r__1;

    /* Local variables */
    int j;
    float h2;
    int n2, lo;

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    ssdrv1_tv_(nx, &v[1], &w[1]);
    saxpy_(&nx, &s_m1, &v[nx + 1], &c__1, &w[1], &c__1);

    i__1 = nx - 1;
    for (j = 2; j <= i__1; ++j)
    {
        lo = (j - 1) * nx;
        ssdrv1_tv_(nx, &v[lo + 1], &w[lo + 1]);
        saxpy_(&nx, &s_m1, &v[lo - nx + 1], &c__1, &w[lo + 1], &c__1);
        saxpy_(&nx, &s_m1, &v[lo + nx + 1], &c__1, &w[lo + 1], &c__1);
    }

    lo = (nx - 1) * nx;
    ssdrv1_tv_(nx, &v[lo + 1], &w[lo + 1]);
    saxpy_(&nx, &s_m1, &v[lo - nx + 1], &c__1, &w[lo + 1], &c__1);

    /*     Scale the vector w by (1/h^2), where h is the mesh size */

    n2 = nx * nx;
    h2 = 1.0f / (float) ((nx + 1) * (nx + 1));
    r__1 = 1.0f / h2;
    sscal_(&n2, &r__1, &w[1], &c__1);
    return 0;
} /* av_ */

/* ------------------------------------------------------------------- */
int ssdrv1_tv_(const int nx, float *x, float *y)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    int j;
    float dd, dl, du;

    /*     Compute the matrix vector multiplication y<---T*x */
    /*     where T is a nx by nx tridiagonal matrix with DD on the */
    /*     diagonal, DL on the subdiagonal, and DU on the superdiagonal. */

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    dd = 4.0f;
    dl = -1.f;
    du = -1.f;

    y[1] = dd * x[1] + du * x[2];
    i__1 = nx - 1;
    for (j = 2; j <= i__1; ++j)
    {
        y[j] = dl * x[j - 1] + dd * x[j] + du * x[j + 1];
    }
    y[nx] = dl * x[nx - 1] + dd * x[nx];
    return 0;
} /* tv_ */

