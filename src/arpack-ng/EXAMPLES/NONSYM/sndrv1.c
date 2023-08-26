/* EXAMPLES\NONSYM\sndrv1.f -- translated by f2c (version 20100827). */

#include <stdlib.h>
#include "arpack_internal.h"

int sndrv1_av_(const int nx, float* v, float* w);
int sndrv1_tv_(const int nx, float* x, float* y);

/**
 * \BeginDoc
 *
 *     Example program to illustrate the idea of reverse communication
 *     for a standard nonsymmetric eigenvalue problem.
 *
 *     We implement example one of ex-nonsym.doc in DOCUMENTS directory
 *
 * \Example-1
 *     ... Suppose we want to solve A*x = lambda*x in regular mode,
 *         where A is obtained from the standard central difference
 *         discretization of the convection-diffusion operator
 *                 (Laplacian u) + rho*(du / dx)
 *         on the unit square [0,1]x[0,1] with zero Dirichlet boundary
 *         condition.
 *
 *     ... OP = A  and  B = I.
 *
 *     ... Assume "call av (nx,x,y)" computes y = A*x.c
 *
 *     ... Use mode 1 of SNAUPD.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Routines called:
 *     snaupd  ARPACK reverse communication interface routine.
 *     sneupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
 *     snrm2   Level 1 BLAS that computes the norm of a vector.
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
    float d[90] /* (3 * MAXNCV) */;
    float workev[90];

    float sigmar, sigmai;

    int j;
    int ierr, nconv;
    int ishfts, maxitr, mode;

    int ipntr[14];
    int iparam[11];
    logical select[30];
    logical first, rvec;

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
    char* which = "SM";

    /* --------------------------------------------------- */
    /* The work array WORKL is used in SNAUPD as           */
    /* workspace.  Its dimension LWORKL is set as          */
    /* illustrated below.  The parameter TOL determines    */
    /* the stopping criterion. If TOL<=0, machine          */
    /* precision is used.  The variable IDO is used for    */
    /* reverse communication, and is initially set to 0.   */
    /* Setting INFO=0 indicates that a random vector is    */
    /* generated in SNAUPD to start the Arnoldi iteration. */
    /* --------------------------------------------------- */

    int lworkl = ncv * ncv * 3 + ncv * 6;
    float tol = 0.0f;
    int ido = 0;
    int info = 0;

    float* ax = (float*)malloc(n * sizeof(float));
    float* resid = (float*)malloc(n * sizeof(float));
    float* v = (float*)malloc(n * ncv * sizeof(float));
    float* workl = (float*)malloc(lworkl * sizeof(float));
    float* workd = (float*)malloc(3 * n * sizeof(float));

    /* ------------------------------------------------- */
    /* This program uses exact shifts with respect to    */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 1 of SNAUPD is used     */
    /* (IPARAM(7) = 1). All these options can be changed */
    /* by the user. For details see the documentation in */
    /* SNAUPD.                                           */
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
    /* Repeatedly call the routine SNAUPD and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    snaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);

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

        sndrv1_av_(nx, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call SNAUPD again. */
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
        /* documentation in SNAUPD. */
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
    /* Post-Process using SNEUPD.                */
    /*                                           */
    /* Computed eigenvalues may be extracted.    */
    /*                                           */
    /* Eigenvectors may also be computed now if  */
    /* desired.  (indicated by rvec = .true.)    */
    /* ----------------------------------------- */

    rvec = TRUE_;

    sneupd_(&rvec, "A", select, d, &d[30], v, &n, &sigmar, &sigmai, workev, bmat, &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &ierr);

    /* --------------------------------------------- */
    /* The real part of the eigenvalue is returned   */
    /* in the first column of the two dimensional    */
    /* array D, and the imaginary part is returned   */
    /* in the second column of D.  The corresponding */
    /* eigenvectors are returned in the first NEV    */
    /* columns of the two dimensional array V if     */
    /* requested.  Otherwise, an orthogonal basis    */
    /* for the invariant subspace corresponding to   */
    /* the eigenvalues in D is returned in V.        */
    /* --------------------------------------------- */

    if (ierr != 0)
    {
        /* ---------------------------------- */
        /* Error condition:                   */
        /* Check the documentation of SNEUPD. */
        /* ---------------------------------- */

        printf(" \n");
        printf(" Error with _neupd info = %d\n", ierr);
        printf(" Check the documentation of _neupd. \n");
        printf(" \n");

        goto EXIT;
    }

    first = TRUE_;
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

        if (d[j + 29] == 0.0f)
        {
            /* ------------------ */
            /* Ritz value is real */
            /* ------------------ */

            sndrv1_av_(nx, &v[k], ax);
            r__1 = -d[j - 1];
            saxpy_(&n, &r__1, &v[k], &c__1, ax, &c__1);
            d[j + 59] = snrm2_(&n, ax, &c__1);
            d[j + 59] /= (r__1 = d[j - 1], dabs(r__1));
        }
        else if (first)
        {
            /* ---------------------- */
            /* Ritz value is complex. */
            /* Residual of one Ritz   */
            /* value of the conjugate */
            /* pair is computed.      */
            /* ---------------------- */

            sndrv1_av_(nx, &v[k], ax);
            r__1 = -d[j - 1];
            saxpy_(&n, &r__1, &v[k], &c__1, ax, &c__1);
            saxpy_(&n, &d[j + 29], &v[j * n], &c__1, ax, &c__1);
            d[j + 59] = snrm2_(&n, ax, &c__1);
            sndrv1_av_(nx, &v[j * n], ax);
            r__1 = -d[j + 29];
            saxpy_(&n, &r__1, &v[k], &c__1, ax, &c__1);
            r__1 = -d[j - 1];
            saxpy_(&n, &r__1, &v[j * n], &c__1, ax, &c__1);
            r__1 = snrm2_(&n, ax, &c__1);
            d[j + 59] = slapy2_(&d[j + 59], &r__1);
            d[j + 59] /= slapy2_(&d[j - 1], &d[j + 29]);
            d[j + 60] = d[j + 59];
            first = FALSE_;
        }
        else
        {
            first = TRUE_;
        }
    }

    /* --------------------------- */
    /* Display computed residuals. */
    /* --------------------------- */

    smout_(&nconv, &c__3, d, &c__30, &c_n6, "Ritz values (Real,Imag) and relative residuals");

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
    printf(" _NDRV1 \n");
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
    /* Done with program sndrv1. */
    /* ------------------------- */

    return ierr;
}

/* ========================================================================== */

/*     matrix vector subroutine */

/*     The matrix used is the 2 dimensional convection-diffusion */
/*     operator discretized using central difference. */

int sndrv1_av_(const int nx, float *v, float *w)
{
    /* System generated locals */
    int i__1;
    float r__1;

    /* Local variables */
    int j;
    float h2;
    int lo;

    /*     Computes w <--- OP*v, where OP is the nx*nx by nx*nx block */
    /*     tridiagonal matrix */

    /*                  | T -I          | */
    /*                  |-I  T -I       | */
    /*             OP = |   -I  T       | */
    /*                  |        ...  -I| */
    /*                  |           -I T| */

    /*     derived from the standard central difference discretization */
    /*     of the 2 dimensional convection-diffusion operator */
    /*     (Laplacian u) + rho*(du/dx) on a unit square with zero boundary */
    /*     condition. */

    /*     When rho*h/2 <= 1, the discrete convection-diffusion operator */
    /*     has real eigenvalues.  When rho*h/2 > 1, it has COMPLEX */
    /*     eigenvalues. */

    /*     The subroutine TV is called to compute y<---T*x. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    h2 = 1.0f / (float) ((nx + 1) * (nx + 1));

    sndrv1_tv_(nx, &v[1], &w[1]);
    r__1 = -1.0f / h2;
    saxpy_(&nx, &r__1, &v[nx + 1], &c__1, &w[1], &c__1);

    i__1 = nx - 1;
    for (j = 2; j <= i__1; ++j)
    {
        lo = (j - 1) * nx;
        sndrv1_tv_(nx, &v[lo + 1], &w[lo + 1]);
        r__1 = -1.0f / h2;
        saxpy_(&nx, &r__1, &v[lo - nx + 1], &c__1, &w[lo + 1], &c__1);
        r__1 = -1.0f / h2;
        saxpy_(&nx, &r__1, &v[lo + nx + 1], &c__1, &w[lo + 1], &c__1);
    }

    lo = (nx - 1) * nx;
    sndrv1_tv_(nx, &v[lo + 1], &w[lo + 1]);
    r__1 = -1.0f / h2;
    saxpy_(&nx, &r__1, &v[lo - nx + 1], &c__1, &w[lo + 1], &c__1);

    return 0;
} /* av_ */

/* ========================================================================= */
int sndrv1_tv_(const int nx, float *x, float *y)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    float h;
    int j;
    float h2, dd, dl, du;

    /*     Compute the matrix vector multiplication y<---T*x */
    /*     where T is a nx by nx tridiagonal matrix with DD on the */
    /*     diagonal, DL on the subdiagonal, and DU on the superdiagonal. */

    /*     When rho*h/2 <= 1, the discrete convection-diffusion operator */
    /*     has real eigenvalues.  When rho*h/2 > 1, it has COMPLEX */
    /*     eigenvalues. */

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    h = 1.0f / (float) (nx + 1);
    h2 = h * h;
    dd = 4.0f / h2;
    dl = -1.0f / h2 - 0.0f / h;
    du = -1.0f / h2 + 0.0f / h;

    y[1] = dd * x[1] + du * x[2];
    i__1 = nx - 1;
    for (j = 2; j <= i__1; ++j)
    {
        y[j] = dl * x[j - 1] + dd * x[j] + du * x[j + 1];
    }
    y[nx] = dl * x[nx - 1] + dd * x[nx];
    return 0;
} /* tv_ */

