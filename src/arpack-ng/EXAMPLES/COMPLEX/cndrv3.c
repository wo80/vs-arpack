/* EXAMPLES\COMPLEX\cndrv3.f -- translated by f2c (version 20100827). */

#include <stdlib.h>
#include "arpack.h"

int cndrv3_av_(const int n, a_fcomplex* v, a_fcomplex* w);
int cndrv3_mv_(const int n, a_fcomplex* v, a_fcomplex* w);

extern int smout_(const int, const int, const float*, const int, const int, const char*);

static int i_one = 1;

/**
 * \BeginDoc
 *
 *     Simple program to illustrate the idea of reverse communication
 *     in inverse mode for a generalized complex nonsymmetric eigenvalue
 *     problem.
 *
 *     We implement example three of ex-complex.doc in DOCUMENTS directory
 *
 * \Example-3
 *     ... Suppose we want to solve A*x = lambda*B*x in regular mode,
 *         where A and B are derived from the finite element discretization
 *         of the 1-dimensional convection-diffusion operator
 *                   (d^2u/dx^2) + rho*(du/dx)
 *         on the interval [0,1] with zero boundary condition using
 *         piecewise linear elements.
 *
 *     ... OP = inv[M]*A  and  B = M.
 *
 *     ... Use mode 2 of CNAUPD.
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
 *     scnrm2  Level 1 BLAS that computes the norm of a vector.
 *     av      Matrix vector multiplication routine that computes A*x.
 *     mv      Matrix vector multiplication routine that computes M*x.
 *
 * \EndLib
 */
int main()
{
    /* System generated locals */
    int i__1, i__2;
    a_fcomplex q__1;

    /* Local variables */
    a_fcomplex d[25];
    a_fcomplex workev[50];
    float rd[75] /* (3 * MAXNCV) */;
    float rwork[256];

    a_fcomplex sigma, h;

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
    /* 'G').  NEV is the number of eigenvalues to be      */
    /* approximated.  The user can modify NEV, NCV, WHICH */
    /* to solve problems of different sizes, and to get   */
    /* different parts of the spectrum.  However, The     */
    /* following conditions must be satisfied:            */
    /*                    N <= MAXN,                      */
    /*                  NEV <= MAXNEV,                    */
    /*              NEV + 2 <= NCV <= MAXNCV              */
    /* -------------------------------------------------- */

    int n = 100;
    int nev = 4;
    int ncv = 20;
    if (n > MAXN)
    {
        printf(" ERROR with _NDRV3: N is greater than MAXN \n");
        return 0;
    }
    else if (nev > MAXNEV)
    {
        printf(" ERROR with _NDRV3: NEV is greater than MAXNEV \n");
        return 0;
    }
    else if (ncv > MAXNCV)
    {
        printf(" ERROR with _NDRV3: NCV is greater than MAXNCV \n");
        return 0;
    }
    char* bmat = "G";
    char* which = "LM";
    sigma.r = 0.0f, sigma.i = 0.0f;

    /* --------------------------------------------------- */
    /* The matrix M is chosen to be the symmetric tri-     */
    /* diagonal matrix with 4 on the diagonal and 1 on the */
    /* off diagonals. It is factored by LAPACK subroutine  */
    /* cgttrf.                                             */
    /* --------------------------------------------------- */

    a_fcomplex* du = (a_fcomplex*)malloc(n * sizeof(a_fcomplex));
    a_fcomplex* dd = (a_fcomplex*)malloc(n * sizeof(a_fcomplex));
    a_fcomplex* dl = (a_fcomplex*)malloc(n * sizeof(a_fcomplex));
    a_fcomplex* du2 = (a_fcomplex*)malloc(n * sizeof(a_fcomplex));

    h.r = 1.0f / (float)(n + 1);
    h.i = 0.0f;
    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        i__2 = j - 1;
        dl[i__2].r = h.r * 1.0f, dl[i__2].i = h.i * 1.0f;
        dd[i__2].r = h.r * 4.0f, dd[i__2].i = h.i * 4.0f;
        du[i__2].r = h.r * 1.0f, du[i__2].i = h.i * 1.0f;
    }
    i__1 = n - 1;
    dd[i__1].r = h.r * 4.0f, dd[i__1].i = h.i * 4.0f;

    cgttrf_(&n, dl, dd, du, du2, ipiv, &ierr);
    if (ierr != 0)
    {
        printf(" \n");
        printf(" ERROR with _gttrf. \n");
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

    a_fcomplex* ax = (a_fcomplex*)malloc(n * sizeof(a_fcomplex));
    a_fcomplex* mx = (a_fcomplex*)malloc(n * sizeof(a_fcomplex));
    a_fcomplex* resid = (a_fcomplex*)malloc(n * sizeof(a_fcomplex));
    a_fcomplex* v = (a_fcomplex*)malloc(n * ncv * sizeof(a_fcomplex));
    a_fcomplex* workl = (a_fcomplex*)malloc(lworkl * sizeof(a_fcomplex));
    a_fcomplex* workd = (a_fcomplex*)malloc(3 * n * sizeof(a_fcomplex));

    /* ------------------------------------------------- */
    /* This program uses exact shifts with respect to    */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 2 of CNAUPD is used     */
    /* (IPARAM(7) = 2).  All these options can be        */
    /* changed by the user. For details, see the         */
    /* documentation in CNAUPD.                          */
    /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode = 2;

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
        /* -------------------------------------- */
        /* Perform  y <--- OP*x = inv[M]*A*x      */
        /* The user should supply his/her own     */
        /* matrix vector routine and a linear     */
        /* system solver.  The matrix-vector      */
        /* subroutine should take workd(ipntr(1)) */
        /* as input, and the final result should  */
        /* be returned to workd(ipntr(2)).        */
        /* -------------------------------------- */

        cndrv3_av_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
        cgttrs_("N", &n, &i_one, dl, dd, du, du2, ipiv, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" ERROR with _gttrs. \n");
            printf(" \n");
            goto EXIT;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call CNAUPD again. */
        /* --------------------------------------- */

        goto L10;
    }
    else if (ido == 2)
    {
        /* ----------------------------------- */
        /*        Perform  y <--- M*x          */
        /* The matrix vector multiplication    */
        /* routine should take workd(ipntr(1)) */
        /* as input and return the result to   */
        /* workd(ipntr(2)).                    */
        /* ----------------------------------- */

        cndrv3_mv_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call CNAUPD again. */
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
        /* Error message. Check the */
        /* documentation in CNAUPD. */
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
        printf(" Check the documentation of _neupd\n");
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
        /*  ||  A*x - lambda*M*x ||  */
        /*                           */
        /* for the NCONV accurately  */
        /* computed eigenvalues and  */
        /* eigenvectors.  (iparam(5) */
        /* indicates how many are    */
        /* accurate to the requested */
        /* tolerance)                */
        /* ------------------------- */

        cndrv3_av_(n, &v[k], ax);
        cndrv3_mv_(n, &v[k], mx);
        i__2 = j - 1;
        q__1.r = -d[i__2].r, q__1.i = -d[i__2].i;
        caxpy_(&n, &q__1, mx, &i_one, ax, &i_one);
        i__2 = j - 1;
        rd[j - 1] = d[i__2].r;
        rd[j + 24] = d[i__2].i;
        rd[j + 49] = scnrm2_(&n, ax, &i_one);
        rd[j + 49] /= slapy2_(&rd[j - 1], &rd[j + 24]);

    }

    /* --------------------------- */
    /* Display computed residuals. */
    /* --------------------------- */

    smout_(nconv, 3, rd, MAXNCV, -6, "Ritz values (Real, Imag) and relative residuals");

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
    printf("_NDRV3 \n");
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
 * Compute the matrix vector multiplication y<---A*x
 * where A is the stiffness matrix formed by using piecewise linear
 * elements on [0,1].
 */
int cndrv3_av_(const int n, a_fcomplex *v, a_fcomplex *w)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    int j;
    float h, s, dd, dl, du;

    /* Function Body */
    h = 1.0f / (float)(n + 1);
    s = 10.0f / 2.0f;
    dd = 2.0f / h;
    dl = -1.0f / h - s;
    du = -1.0f / h + s;

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

/**
 * Compute the matrix vector multiplication y<---M*x
 * where M is the mass matrix formed by using piecewise linear elements
 * on [0,1].
 */
int cndrv3_mv_(const int n, a_fcomplex *v, a_fcomplex *w)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    a_fcomplex h;
    int j;

    /* Function Body */
    w[0].r = v[0].r * 4.0f + v[1].r * 1.0f;
    w[0].i = v[0].i * 4.0f + v[1].i * 1.0f;

    i__1 = n - 1;
    for (j = 1; j < i__1; ++j)
    {
        i__2 = j - 1;
        i__3 = j + 1;
        w[j].r = (v[i__2].r * 1.0f + v[j].r * 4.0f) + v[i__3].r * 1.0f;
        w[j].i = (v[i__2].i * 1.0f + v[j].i * 4.0f) + v[i__3].i * 1.0f;
    }
    i__1 = n - 1;
    i__2 = n - 2;
    w[i__1].r = v[i__2].r * 1.0f + v[i__1].r * 4.0f;
    w[i__1].i = v[i__2].i * 1.0f + v[i__1].i * 4.0f;

    h.r = 1.0f / (float)(n + 1);
    h.i = 0.0f;
    cscal_(&n, &h, w, &i_one);
    return 0;
} /* mv_ */
