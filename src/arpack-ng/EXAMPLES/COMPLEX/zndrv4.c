/* EXAMPLES\COMPLEX\zndrv4.f -- translated by f2c (version 20100827). */

#include <stdlib.h>
#include <stdio.h>
#include "arpack.h"
#include "lapack.h"

int zndrv4_av_(const int n, a_dcomplex* v, a_dcomplex* w);
int zndrv4_mv_(const int n, a_dcomplex* v, a_dcomplex* w);

extern int dmout_(const int, const int, const double*, const int, const int, const char*);

static int i_one = 1;

#define RHO 10.0

/**
 * \BeginDoc
 *
 *     Simple program to illustrate the idea of reverse communication
 *     in shift and invert mode for a generalized complex nonsymmetric
 *     eigenvalue problem.
 *
 *     We implement example four of ex-complex.doc in DOCUMENTS directory
 *
 * \Example-4
 *     ... Suppose we want to solve A*x = lambda*B*x in shift-invert mode,
 *         where A and B are derived from a finite element discretization
 *         of a 1-dimensional convection-diffusion operator
 *                         (d^2u/dx^2) + rho*(du/dx)
 *         on the interval [0,1] with zero boundary condition using
 *         piecewise linear elements.
 *
 *     ... where the shift sigma is a complex number.
 *
 *     ... OP = inv[A-SIGMA*M]*M  and  B = M.
 *
 *     ... Use mode 3 of ZNAUPD .
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Routines called:
 *     znaupd   ARPACK reverse communication interface routine.
 *     zneupd   ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     zgttrf   LAPACK tridiagonal factorization routine.
 *     zgttrs   LAPACK tridiagonal solve routine.
 *     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     zaxpy    Level 1 BLAS that computes y <- alpha*x+y.
 *     zcopy    Level 1 BLAS that copies one vector to another.
 *     dznrm2   Level 1 BLAS that computes the norm of a complex vector.
 *     av      Matrix vector multiplication routine that computes A*x.
 *     mv      Matrix vector multiplication routine that computes M*x.
 *
 * \EndLib
 */
int main()
{
    /* System generated locals */
    int i__1, i__2;
    a_dcomplex z__1;

    /* Local variables */
    a_dcomplex d[25];
    a_dcomplex workev[50];
    double rd[75] /* (3 * MAXNCV) */;
    double rwork[256];

    double h, s;
    a_dcomplex sigma, s1, s2, s3;

    int j;
    int ierr, nconv;
    int ishfts, maxitr, mode;
    int ipiv[256];
    int ipntr[14];
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
    /* 'G').  NEV is the number of eigenvalues (closest   */
    /* to SIGMAR) to be approximated.  Since the          */
    /* shift-invert mode is used,  WHICH is set to 'LM'.  */
    /* The user can modify NEV, NCV, SIGMA to solve       */
    /* problems of different sizes, and to get different  */
    /* parts of the spectrum.  However, The following     */
    /* conditions must be satisfied:                      */
    /*                     N <= MAXN,                     */
    /*                   NEV <= MAXNEV,                   */
    /*               NEV + 2 <= NCV <= MAXNCV             */
    /* -------------------------------------------------- */

    int n = 100;
    int nev = 4;
    int ncv = 20;
    if (n > MAXN)
    {
        printf(" ERROR with _NDRV4: N is greater than MAXN \n");
        return 0;
    }
    else if (nev > MAXNEV)
    {
        printf(" ERROR with _NDRV4: NEV is greater than MAXNEV \n");
        return 0;
    }
    else if (ncv > MAXNCV)
    {
        printf(" ERROR with _NDRV4: NCV is greater than MAXNCV \n");
        return 0;
    }
    char* bmat = "G";
    char* which = "LM";
    sigma.r = 1.0, sigma.i = 0.0;

    /* ------------------------------------------------ */
    /* Construct C = A - SIGMA*M in COMPLEX arithmetic. */
    /* Factor C in COMPLEX arithmetic (using LAPACK     */
    /* subroutine zgttrf ). The matrix A is chosen to be */
    /* the tridiagonal matrix derived from the standard */
    /* central difference discretization of the 1-d     */
    /* convection-diffusion operator u``+ rho*u` on the */
    /* interval [0, 1] with zero Dirichlet boundary     */
    /* condition.  The matrix M is chosen to be the     */
    /* symmetric tridiagonal matrix with 4.0 on the     */
    /* diagonal and 1.0 on the off-diagonals.           */
    /* ------------------------------------------------ */

    a_dcomplex* du = (a_dcomplex*)malloc(n * sizeof(a_dcomplex));
    a_dcomplex* dd = (a_dcomplex*)malloc(n * sizeof(a_dcomplex));
    a_dcomplex* dl = (a_dcomplex*)malloc(n * sizeof(a_dcomplex));
    a_dcomplex* du2 = (a_dcomplex*)malloc(n * sizeof(a_dcomplex));

    h = 1.0 / (double)(n + 1);
    s = RHO / 2.0;

    s1.r = (-1.0 / h - s) - sigma.r * h / 6.0;
    s1.i = -sigma.i * h / 6.0;

    s2.r = 2.0 / h - (sigma.r * 4.0 * h / 6.0);
    s2.i = -sigma.i * 4.0 * h / 6.0;

    s3.r = (-1.0 / h + s) - (sigma.r * h / 6.0);
    s3.i = -sigma.i * h / 6.0;

    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        i__2 = j - 1;
        dl[i__2].r = s1.r, dl[i__2].i = s1.i;
        dd[i__2].r = s2.r, dd[i__2].i = s2.i;
        du[i__2].r = s3.r, du[i__2].i = s3.i;
    }
    i__1 = n - 1;
    dd[i__1].r = s2.r, dd[i__1].i = s2.i;

    zgttrf_(&n, dl, dd, du, du2, ipiv, &ierr);
    if (ierr != 0)
    {
        printf(" \n");
        printf(" ERROR with _gttrf in _NDRV4.\n");
        printf(" \n");
        return 0;
    }

    /* --------------------------------------------------- */
    /* The work array WORKL is used in ZNAUPD  as           */
    /* workspace.  Its dimension LWORKL is set as          */
    /* illustrated below.  The parameter TOL determines    */
    /* the stopping criterion. If TOL<=0, machine          */
    /* precision is used.  The variable IDO is used for    */
    /* reverse communication, and is initially set to 0.   */
    /* Setting INFO=0 indicates that a random vector is    */
    /* generated in ZNAUPD  to start the Arnoldi iteration. */
    /* --------------------------------------------------- */

    /* Computing 2nd power */
    int lworkl = ncv * ncv * 3 + ncv * 5;
    double tol = 0.0;
    int ido = 0;
    int info = 0;

    a_dcomplex* ax = (a_dcomplex*)malloc(n * sizeof(a_dcomplex));
    a_dcomplex* mx = (a_dcomplex*)malloc(n * sizeof(a_dcomplex));
    a_dcomplex* resid = (a_dcomplex*)malloc(n * sizeof(a_dcomplex));
    a_dcomplex* v = (a_dcomplex*)malloc(n * ncv * sizeof(a_dcomplex));
    a_dcomplex* workl = (a_dcomplex*)malloc(lworkl * sizeof(a_dcomplex));
    a_dcomplex* workd = (a_dcomplex*)malloc(3 * n * sizeof(a_dcomplex));

    /* ------------------------------------------------- */
    /* This program uses exact shifts with respect to    */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed. Mode 3 of ZNAUPD  is used      */
    /* (IPARAM(7) = 3).  All these options can be        */
    /* changed by the user. For details see the          */
    /* documentation in ZNAUPD .                          */
    /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode = 3;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

    /* ---------------------------------------- */
    /* M A I N   L O O P(Reverse communication) */
    /* ---------------------------------------- */

L20:

    /* ------------------------------------------- */
    /* Repeatedly call the routine ZNAUPD  and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    znaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, rwork, &info);

    if (ido == -1)
    {
        /* ----------------------------------------- */
        /* Perform  y <--- OP*x = inv[A-SIGMA*M]*M*x */
        /* to force starting vector into the range   */
        /* of OP.   The user should supply his/her   */
        /* own matrix vector multiplication routine  */
        /* and a linear system solver.  The matrix   */
        /* vector multiplication routine should take */
        /* workd(ipntr(1)) as the input. The final   */
        /* result should be returned to              */
        /* workd(ipntr(2)).                          */
        /* ----------------------------------------- */

        zndrv4_mv_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
        zgttrs_("N", &n, &i_one, dl, dd, du, du2, ipiv, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" ERROR with _gttrs in _NDRV4.\n");
            printf(" \n");
            goto EXIT;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call ZNAUPD  again. */
        /* --------------------------------------- */

        goto L20;
    }
    else if (ido == 1)
    {
        /* --------------------------------------- */
        /* Perform y <-- OP*x = inv[A-sigma*M]*M*x */
        /* M*x has been saved in workd(ipntr(3)).  */
        /* The user only need the linear system    */
        /* solver here that takes workd(ipntr(3))  */
        /* as input, and returns the result to     */
        /* workd(ipntr(2)).                        */
        /* --------------------------------------- */

        zcopy_(&n, &workd[ipntr[2] - 1], &i_one, &workd[ipntr[1] - 1], &i_one);
        zgttrs_("N", &n, &i_one, dl, dd, du, du2, ipiv, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" ERROR with _gttrs in _NDRV4.\n");
            printf(" \n");
            goto EXIT;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call ZNAUPD  again. */
        /* --------------------------------------- */

        goto L20;
    }
    else if (ido == 2)
    {
        /* ------------------------------------------- */
        /*          Perform  y <--- M*x                */
        /* Need matrix vector multiplication routine   */
        /* here that takes workd(ipntr(1)) as input    */
        /* and returns the result to workd(ipntr(2)).  */
        /* ------------------------------------------- */

        zndrv4_mv_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call ZNAUPD  again. */
        /* --------------------------------------- */

        goto L20;
    }

    /* --------------------------------------- */
    /* Either we have convergence, or there is */
    /* an error.                               */
    /* --------------------------------------- */

    if (info < 0)
    {
        /* -------------------------- */
        /*  Error message, check the  */
        /*  documentation in ZNAUPD    */
        /* -------------------------- */

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

        zndrv4_av_(n, &v[k], ax);
        zndrv4_mv_(n, &v[k], mx);
        i__2 = j - 1;
        z__1.r = -d[i__2].r, z__1.i = -d[i__2].i;
        zaxpy_(&n, &z__1, mx, &i_one, ax, &i_one);
        i__2 = j - 1;
        rd[j - 1] = d[i__2].r;
        rd[j + 24] = d[i__2].i;
        rd[j + 49] = dznrm2_(&n, ax, &i_one);
        rd[j + 49] /= dlapy2_(&rd[j - 1], &rd[j + 24]);

    }

    /* --------------------------- */
    /* Display computed residuals. */
    /* --------------------------- */

    dmout_(nconv, 3, rd, MAXNCV, -6, "Ritz values (Real, Imag) and direct residuals");

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
    printf("_NDRV4 \n");
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
    free(mx);
    free(resid);
    free(v);
    free(workl);
    free(workd);

    return ierr;
}

/** Matrix vector multiplication subroutine.
 *
 * Compute the matrix vector multiplication y<---M*x
 * where M is a n by n symmetric tridiagonal matrix with 4 on the
 * diagonal, 1 on the subdiagonal and superdiagonal.
 */
int zndrv4_mv_(const int n, a_dcomplex *v, a_dcomplex *w)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    a_dcomplex h;
    int j;

    /* Function Body */
    w[0].r = (v[0].r * 4.0 + v[1].r) / 6.0;
    w[0].i = (v[0].i * 4.0 + v[1].i) / 6.0;
    i__1 = n - 1;
    for (j = 1; j < i__1; ++j)
    {
        i__2 = j - 1;
        i__3 = j + 1;
        w[j].r = (v[i__2].r + v[j].r * 4.0 + v[i__3].r) / 6.0;
        w[j].i = (v[i__2].i + v[j].i * 4.0 + v[i__3].i) / 6.0;
    }
    i__1 = n - 1;
    i__2 = n - 2;
    w[i__1].r = (v[i__2].r + v[i__1].r * 4.0) / 6.0;
    w[i__1].i = (v[i__2].i + v[i__1].i * 4.0) / 6.0;

    h.r = 1.0 / (double)(n + 1);
    h.i = 0.0;
    zscal_(&n, &h, w, &i_one);
    return 0;
} /* mv_ */

int zndrv4_av_(const int n, a_dcomplex *v, a_dcomplex *w)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    int j;
    double h, s, dd, dl, du;

    /* Function Body */
    h = 1.0 / (double)(n + 1);
    s = RHO / 2.0;
    dd = 2.0 / h;
    dl = -1.0 / h - s;
    du = -1.0 / h + s;

    w[0].r = dd * v[0].r + du * v[1].r;
    w[0].i = dd * v[0].i + du * v[1].i;

    i__1 = n - 1;
    for (j = 1; j < i__1; ++j)
    {
        i__2 = j - 1;
        i__3 = j + 1;
        w[j].r = dl * v[i__2].r + dd * v[j].r + du * v[i__3].r;
        w[j].i = dl * v[i__2].i + dd * v[j].i + du * v[i__3].i;
    }
    i__1 = n - 1;
    i__2 = n - 2;
    w[i__1].r = dl * v[i__2].r + dd * v[i__1].r;
    w[i__1].i = dl * v[i__2].i + dd * v[i__1].i;
    return 0;
} /* av_ */
