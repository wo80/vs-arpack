/* EXAMPLES\COMPLEX\cndrv4.f -- translated by f2c (version 20100827). */

#include "arpack.h"

struct
{
    complex rho;
} convct_;

#define convct_1 convct_

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
 *     scnrm2  Level 1 BLAS that computes the norm of a complex vector.
 *     av      Matrix vector multiplication routine that computes A*x.
 *     mv      Matrix vector multiplication routine that computes M*x.
 *
 * \EndLib
 */
int cndrv4()
{
    /* System generated locals */
    int32_t i__1, i__2;
    complex q__1, q__2, q__3, q__4, q__5, q__6;

    void c_div(complex *, complex *, complex *);
    double r_imag(complex *);

    /* Local variables */
    complex d[25], h;
    int32_t j, n;
    complex s, s1, s2, s3, dd[256], dl[256];
    float rd[75]	/* was [25][3] */;
    complex ax[256], du[256];
    complex mx[256], du2[256];
    int32_t ido, ncv, nev;
    float tol;
    char* bmat;
    int32_t mode, info;
    bool rvec;
    int32_t ierr, ipiv[256];
    complex sigma;
    char* which;
    int32_t nconv;
    complex *v	/* was [256][25] */;
    complex *resid;
    complex *workd;
    complex *workl;
    int32_t ipntr[14];
    float rwork[256];
    int32_t iparam[11];
    bool select[25];
    int32_t ishfts, maxitr;
    int32_t lworkl;
    complex workev[50];

    resid = (complex*)malloc(256 * sizeof(complex));
    v = (complex*)malloc(6400 * sizeof(complex));
    workl = (complex*)malloc(2000 * sizeof(complex));
    workd = (complex*)malloc(768 * sizeof(complex));

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

    n = 100;
    nev = 4;
    ncv = 20;
    if (n > 256)
    {
        printf(" ERROR with _NDRV4: N is greater than MAXN \n");
        return 0;
    }
    else if (nev > 10)
    {
        printf(" ERROR with _NDRV4: NEV is greater than MAXNEV \n");
        return 0;
    }
    else if (ncv > 25)
    {
        printf(" ERROR with _NDRV4: NCV is greater than MAXNCV \n");
        return 0;
    }
    bmat = "G";
    which = "LM";
    sigma.r = 1.f, sigma.i = 0.f;

    /* ------------------------------------------------ */
    /* Construct C = A - SIGMA*M in COMPLEX arithmetic. */
    /* Factor C in COMPLEX arithmetic (using LAPACK     */
    /* subroutine cgttrf). The matrix A is chosen to be */
    /* the tridiagonal matrix derived from the standard */
    /* central difference discretization of the 1-d     */
    /* convection-diffusion operator u``+ rho*u` on the */
    /* interval [0, 1] with zero Dirichlet boundary     */
    /* condition.  The matrix M is chosen to be the     */
    /* symmetric tridiagonal matrix with 4.0 on the     */
    /* diagonal and 1.0 on the off-diagonals.           */
    /* ------------------------------------------------ */

    convct_1.rho.r = 10.f, convct_1.rho.i = 0.f;
    i__1 = n + 1;
    q__2.r = (float) i__1, q__2.i = 0.f;
    c_div(&q__1, &c_one, &q__2);
    h.r = q__1.r, h.i = q__1.i;
    c_div(&q__1, &convct_1.rho, &c_b3);
    s.r = q__1.r, s.i = q__1.i;

    q__4.r = -1.f, q__4.i = -0.f;
    c_div(&q__3, &q__4, &h);
    q__2.r = q__3.r - s.r, q__2.i = q__3.i - s.i;
    q__6.r = sigma.r * h.r - sigma.i * h.i, q__6.i = sigma.r * h.i +
             sigma.i * h.r;
    c_div(&q__5, &q__6, &c_b5);
    q__1.r = q__2.r - q__5.r, q__1.i = q__2.i - q__5.i;
    s1.r = q__1.r, s1.i = q__1.i;
    c_div(&q__2, &c_b3, &h);
    q__5.r = sigma.r * 4.f - sigma.i * 0.f, q__5.i = sigma.i * 4.f + sigma.r *
             0.f;
    q__4.r = q__5.r * h.r - q__5.i * h.i, q__4.i = q__5.r * h.i +
             q__5.i * h.r;
    c_div(&q__3, &q__4, &c_b5);
    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
    s2.r = q__1.r, s2.i = q__1.i;
    q__4.r = -1.f, q__4.i = -0.f;
    c_div(&q__3, &q__4, &h);
    q__2.r = q__3.r + s.r, q__2.i = q__3.i + s.i;
    q__6.r = sigma.r * h.r - sigma.i * h.i, q__6.i = sigma.r * h.i +
             sigma.i * h.r;
    c_div(&q__5, &q__6, &c_b5);
    q__1.r = q__2.r - q__5.r, q__1.i = q__2.i - q__5.i;
    s3.r = q__1.r, s3.i = q__1.i;

    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        i__2 = j - 1;
        dl[i__2].r = s1.r, dl[i__2].i = s1.i;
        i__2 = j - 1;
        dd[i__2].r = s2.r, dd[i__2].i = s2.i;
        i__2 = j - 1;
        du[i__2].r = s3.r, du[i__2].i = s3.i;
        /* L10: */
    }
    i__1 = n - 1;
    dd[i__1].r = s2.r, dd[i__1].i = s2.i;

    cgttrf_(&n, dl, dd, du, du2, ipiv, &ierr);
    if (ierr != 0)
    {
        printf(" \n");
        printf(" ERROR with _gttrf in _NDRV4.\n");
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
    i__1 = ncv;
    lworkl = i__1 * i__1 * 3 + ncv * 5;
    tol = 0.f;
    ido = 0;
    info = 0;

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

    /* ---------------------------------------- */
    /* M A I N   L O O P(Reverse communication) */
    /* ---------------------------------------- */

L20:

    /* ------------------------------------------- */
    /* Repeatedly call the routine CNAUPD and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    cnaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, rwork, &info);

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

        cndrv4_mv_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
        cgttrs_("N", &n, &c__1, dl, dd, du, du2, ipiv, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" ERROR with _gttrs in _NDRV4.\n");
            printf(" \n");
            return ierr;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call CNAUPD again. */
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

        ccopy_(&n, &workd[ipntr[2] - 1], &c__1, &workd[ipntr[1] - 1], &c__1);
        cgttrs_("N", &n, &c__1, dl, dd, du, du2, ipiv, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" ERROR with _gttrs in _NDRV4.\n");
            printf(" \n");
            return ierr;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call CNAUPD again. */
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

        cndrv4_mv_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

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

        /* -------------------------- */
        /*  Error message, check the  */
        /*  documentation in CNAUPD   */
        /* -------------------------- */

        printf(" \n");
        printf(" Error with _naupd info = %d\n", info);
        printf(" Check the documentation of _naupd.\n");
        printf(" \n");

    }
    else
    {

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

        cneupd_(&rvec, "A", select, d, v, &c__256, &sigma, workev, bmat, &n,which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);

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

        }
        else
        {

            nconv = iparam[4];
            i__1 = nconv;
            for (j = 1; j <= i__1; ++j)
            {

                cndrv4_av_(&n, &v[(j << 8) - 256], ax);
                cndrv4_mv_(&n, &v[(j << 8) - 256], mx);
                i__2 = j - 1;
                q__1.r = -d[i__2].r, q__1.i = -d[i__2].i;
                caxpy_(&n, &q__1, mx, &c__1, ax, &c__1);
                i__2 = j - 1;
                rd[j - 1] = d[i__2].r;
                rd[j + 24] = r_imag(&d[j - 1]);
                rd[j + 49] = scnrm2_(&n, ax, &c__1);
                rd[j + 49] /= slapy2_(&rd[j - 1], &rd[j + 24]);
                /* L80: */
            }

            /* --------------------------- */
            /* Display computed residuals. */
            /* --------------------------- */

            smout_(&nconv, &c__3, rd, &c__25, &c_n6, "Ritz values (Real, Imag) and direct residuals");

        }

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

    }

    free(resid);
    free(v);
    free(workl);
    free(workd);

    return 0;
} /* MAIN__ */

/* ========================================================================== */

/*     matrix vector multiplication subroutine */

int cndrv4_mv_(int32_t *n, complex *v, complex *w)
{
    /* System generated locals */
    int32_t i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2, q__3, q__4, q__5, q__6;

    /* Builtin functions */
    void c_div(complex *, complex *, complex *);

    /* Local variables */
    complex h;
    int32_t j;

    /*     Compute the matrix vector multiplication y<---M*x */
    /*     where M is a n by n symmetric tridiagonal matrix with 4 on the */
    /*     diagonal, 1 on the subdiagonal and superdiagonal. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    q__3.r = v[1].r * 4.f - v[1].i * 0.f, q__3.i = v[1].i * 4.f + v[1].r *
             0.f;
    q__4.r = v[2].r * 1.f - v[2].i * 0.f, q__4.i = v[2].i * 1.f + v[2].r *
             0.f;
    q__2.r = q__3.r + q__4.r, q__2.i = q__3.i + q__4.i;
    c_div(&q__1, &q__2, &c_b5);
    w[1].r = q__1.r, w[1].i = q__1.i;
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        i__2 = j;
        i__3 = j - 1;
        q__4.r = v[i__3].r * 1.f - v[i__3].i * 0.f, q__4.i = v[i__3].i * 1.f
                 + v[i__3].r * 0.f;
        i__4 = j;
        q__5.r = v[i__4].r * 4.f - v[i__4].i * 0.f, q__5.i = v[i__4].i * 4.f
                 + v[i__4].r * 0.f;
        q__3.r = q__4.r + q__5.r, q__3.i = q__4.i + q__5.i;
        i__5 = j + 1;
        q__6.r = v[i__5].r * 1.f - v[i__5].i * 0.f, q__6.i = v[i__5].i * 1.f
                 + v[i__5].r * 0.f;
        q__2.r = q__3.r + q__6.r, q__2.i = q__3.i + q__6.i;
        c_div(&q__1, &q__2, &c_b5);
        w[i__2].r = q__1.r, w[i__2].i = q__1.i;
        /* L40: */
    }
    i__1 = *n;
    i__2 = *n - 1;
    q__3.r = v[i__2].r * 1.f - v[i__2].i * 0.f, q__3.i = v[i__2].i * 1.f + v[
                 i__2].r * 0.f;
    i__3 = *n;
    q__4.r = v[i__3].r * 4.f - v[i__3].i * 0.f, q__4.i = v[i__3].i * 4.f + v[
                 i__3].r * 0.f;
    q__2.r = q__3.r + q__4.r, q__2.i = q__3.i + q__4.i;
    c_div(&q__1, &q__2, &c_b5);
    w[i__1].r = q__1.r, w[i__1].i = q__1.i;

    i__1 = *n + 1;
    q__2.r = (float) i__1, q__2.i = 0.f;
    c_div(&q__1, &c_one, &q__2);
    h.r = q__1.r, h.i = q__1.i;
    cscal_(n, &h, &w[1], &c__1);
    return 0;
} /* mv_ */

/* ------------------------------------------------------------------ */
int cndrv4_av_(int32_t *n, complex *v, complex *w)
{
    /* System generated locals */
    int32_t i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2, q__3, q__4, q__5;

    /* Builtin functions */
    void c_div(complex *, complex *, complex *);

    /* Local variables */
    complex h;
    int32_t j;
    complex s, dd, dl, du;

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    i__1 = *n + 1;
    q__2.r = (float) i__1, q__2.i = 0.f;
    c_div(&q__1, &c_one, &q__2);
    h.r = q__1.r, h.i = q__1.i;
    c_div(&q__1, &convct_1.rho, &c_b3);
    s.r = q__1.r, s.i = q__1.i;
    c_div(&q__1, &c_b3, &h);
    dd.r = q__1.r, dd.i = q__1.i;
    q__3.r = -1.f, q__3.i = -0.f;
    c_div(&q__2, &q__3, &h);
    q__1.r = q__2.r - s.r, q__1.i = q__2.i - s.i;
    dl.r = q__1.r, dl.i = q__1.i;
    q__3.r = -1.f, q__3.i = -0.f;
    c_div(&q__2, &q__3, &h);
    q__1.r = q__2.r + s.r, q__1.i = q__2.i + s.i;
    du.r = q__1.r, du.i = q__1.i;

    q__2.r = dd.r * v[1].r - dd.i * v[1].i, q__2.i = dd.r * v[1].i + dd.i * v[
                 1].r;
    q__3.r = du.r * v[2].r - du.i * v[2].i, q__3.i = du.r * v[2].i + du.i * v[
                 2].r;
    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
    w[1].r = q__1.r, w[1].i = q__1.i;
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        i__2 = j;
        i__3 = j - 1;
        q__3.r = dl.r * v[i__3].r - dl.i * v[i__3].i, q__3.i = dl.r * v[i__3]
                 .i + dl.i * v[i__3].r;
        i__4 = j;
        q__4.r = dd.r * v[i__4].r - dd.i * v[i__4].i, q__4.i = dd.r * v[i__4]
                 .i + dd.i * v[i__4].r;
        q__2.r = q__3.r + q__4.r, q__2.i = q__3.i + q__4.i;
        i__5 = j + 1;
        q__5.r = du.r * v[i__5].r - du.i * v[i__5].i, q__5.i = du.r * v[i__5]
                 .i + du.i * v[i__5].r;
        q__1.r = q__2.r + q__5.r, q__1.i = q__2.i + q__5.i;
        w[i__2].r = q__1.r, w[i__2].i = q__1.i;
        /* L40: */
    }
    i__1 = *n;
    i__2 = *n - 1;
    q__2.r = dl.r * v[i__2].r - dl.i * v[i__2].i, q__2.i = dl.r * v[i__2].i +
             dl.i * v[i__2].r;
    i__3 = *n;
    q__3.r = dd.r * v[i__3].r - dd.i * v[i__3].i, q__3.i = dd.r * v[i__3].i +
             dd.i * v[i__3].r;
    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
    w[i__1].r = q__1.r, w[i__1].i = q__1.i;
    return 0;
} /* av_ */

