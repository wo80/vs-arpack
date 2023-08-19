/* EXAMPLES\NONSYM\dndrv6.f -- translated by f2c (version 20100827). */

#include <stdlib.h>
#include "arpack_internal.h"

/**
 * \BeginDoc
 *
 *     Simple program to illustrate the idea of reverse communication
 *     in shift-invert mode for a generalized nonsymmetric eigenvalue problem.
 *
 *     We implement example six of ex-nonsym.doc in DOCUMENTS directory
 *
 * \Example-6
 *
 *     ... Suppose we want to solve A*x = lambda*B*x in shift-invert mode
 *         The matrix A is the tridiagonal matrix with 2 on the diagonal,
 *         -2 on the subdiagonal and 3 on the superdiagonal.  The matrix M
 *         is the tridiagonal matrix with 4 on the diagonal and 1 on the
 *         off-diagonals.
 *     ... The shift sigma is a complex number (sigmar, sigmai).
 *     ... OP = Imaginary_Part{inv[A-(SIGMAR,SIGMAI)*M]*M  and  B = M.
 *     ... Use mode 4 of DNAUPD.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Routines called:
 *     dnaupd  ARPACK reverse communication interface routine.
 *     dneupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     zgttrf  LAPACK complex matrix factorization routine.
 *     zgttrs  LAPACK complex linear system solve routine.
 *     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
 *     ddot    Level 1 BLAS that computes the dot product of two vectors.
 *     dnrm2   Level 1 BLAS that computes the norm of a vector.
 *     av      Matrix vector subroutine that computes A*x.
 *     mv      Matrix vector subroutine that computes M*x.
 *
 * \EndLib
 */
int dndrv6()
{
    /* System generated locals */
    int i__1, i__2, i__3;
    double d__1, d__2;
    zomplex z__1;

    /* Local variables */
    double d[75]; /* (3 * MAXNCV) */
    double workev[75];

    double denr, deni;
    double numr, numi;
    double sigmar, sigmai;

    zomplex c1, c2, c3;

    int j;
    int ierr, nconv;
    int ishfts, maxitr, mode;

    int ipiv[256];
    int ipntr[14];
    int iparam[11];
    bool select[25];
    bool first, rvec;

    /* Define maximum dimensions for all arrays. */

    const int MAXN   = 256; /* Maximum dimension of the A allowed. */
    const int MAXNEV =  10; /* Maximum NEV allowed */
    const int MAXNCV =  25; /* Maximum NCV allowed */

    /* -------------------------------------------------- */
    /* The number N is the dimension of the matrix.  A    */
    /* generalized eigenvalue problem is solved (BMAT =   */
    /* 'G').  NEV is the number of eigenvalues (closest   */
    /* to the shift (SIGMAR,SIGMAI)) to be approximated.  */
    /* Since the shift-invert mode is used, WHICH is set  */
    /* to 'LM'.  The user can modify NEV, NCV, SIGMA to   */
    /* solve problems of different sizes, and to get      */
    /* different parts of the spectrum. However, The      */
    /* following conditions must be satisfied:            */
    /*                     N <= MAXN,                     */
    /*                   NEV <= MAXNEV,                   */
    /*               NEV + 2 <= NCV <= MAXNCV             */
    /* -------------------------------------------------- */

    int n = 100;
    int nev = 4;
    int ncv = 20;
    if (n > 256)
    {
        printf(" ERROR with _NDRV6: N is greater than MAXN \n");
        return 0;
    }
    else if (nev > 10)
    {
        printf(" ERROR with _NDRV6: NEV is greater than MAXNEV \n");
        return 0;
    }
    else if (ncv > 25)
    {
        printf(" ERROR with _NDRV6: NCV is greater than MAXNCV \n");
        return 0;
    }
    char* bmat = "G";
    char* which = "LM";
    sigmar = 0.4;
    sigmai = 0.6;

    /* -------------------------------------------------- */
    /* Construct C = A - (SIGMAR,SIGMAI)*M in complex     */
    /* arithmetic, and factor C in complex arithmetic     */
    /* (using LAPACK subroutine zgttrf). The matrix A is  */
    /* chosen to be the tridiagonal matrix with -2 on the */
    /* subdiagonal, 2 on the diagonal and 3 on the        */
    /* superdiagonal. The matrix M is chosen to be the    */
    /* symmetric tridiagonal matrix with 4 on the         */
    /* diagonal and 1 on the off-diagonals.               */
    /* -------------------------------------------------- */

    zomplex* cdd = (zomplex*)malloc(n * sizeof(zomplex));
    zomplex* cdl = (zomplex*)malloc(n * sizeof(zomplex));
    zomplex* cdu = (zomplex*)malloc(n * sizeof(zomplex));
    zomplex* cdu2 = (zomplex*)malloc(n * sizeof(zomplex));
    zomplex* ctemp = (zomplex*)malloc(n * sizeof(zomplex));

    d__1 = -2.0 - sigmar;
    d__2 = -sigmai;
    z__1.r = d__1, z__1.i = d__2;
    c1.r = z__1.r, c1.i = z__1.i;
    d__1 = 2. - sigmar * 4.0;
    d__2 = sigmai * -4.0;
    z__1.r = d__1, z__1.i = d__2;
    c2.r = z__1.r, c2.i = z__1.i;
    d__1 = 3.0 - sigmar;
    d__2 = -sigmai;
    z__1.r = d__1, z__1.i = d__2;
    c3.r = z__1.r, c3.i = z__1.i;

    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        i__2 = j - 1;
        cdl[i__2].r = c1.r, cdl[i__2].i = c1.i;
        cdd[i__2].r = c2.r, cdd[i__2].i = c2.i;
        cdu[i__2].r = c3.r, cdu[i__2].i = c3.i;
    }
    cdd[i__1].r = c2.r, cdd[i__1].i = c2.i;

    zgttrf_(&n, cdl, cdd, cdu, cdu2, ipiv, &ierr);
    if (ierr != 0)
    {
        printf(" \n");
        printf(" ERROR with _gttrf in _NDRV6.\n");
        printf(" \n");
        return 0;
    }

    /* --------------------------------------------------- */
    /* The work array WORKL is used in DNAUPD as           */
    /* workspace.  Its dimension LWORKL is set as          */
    /* illustrated below.  The parameter TOL determines    */
    /* the stopping criterion. If TOL<=0, machine          */
    /* precision is used.  The variable IDO is used for    */
    /* reverse communication, and is initially set to 0.   */
    /* Setting INFO=0 indicates that a random vector is    */
    /* generated in DNAUPD to start the Arnoldi iteration. */
    /* --------------------------------------------------- */

    int lworkl = ncv * ncv * 3 + ncv * 6;
    double tol = 0.0;
    int ido = 0;
    int info = 0;

    double* ax = (double*)malloc(n * sizeof(double));
    double* mx = (double*)malloc(n * sizeof(double));
    double* resid = (double*)malloc(n * sizeof(double));
    double* v = (double*)malloc(n * ncv * sizeof(double));
    double* workl = (double*)malloc(lworkl * sizeof(double));
    double* workd = (double*)malloc(3 * n * sizeof(double));

    /* ------------------------------------------------- */
    /* This program uses exact shift with respect to     */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 3 of DNAUPD is used     */
    /* (IPARAM(7) = 3).  All these options can be        */
    /* changed by the user. For details, see the         */
    /* documentation in DNAUPD.                          */
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
    /* Repeatedly call the routine DNAUPD and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    dnaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);

    if (ido == -1)
    {
        /* ---------------------------------------------------------- */
        /*                           Perform                          */
        /* y <--- OP*x = Imaginary_Part{inv[A-(SIGMAR,SIGMAI)*M]*M*x} */
        /* to force starting vector into the range of OP. The user    */
        /* should supply his/her own matrix vector multiplication     */
        /* routine and a complex linear system solver.  The matrix    */
        /* vector multiplication routine should take workd(ipntr(1))  */
        /* as the input. The final result (a real vector) should be   */
        /* returned to workd(ipntr(2)).                               */
        /* ---------------------------------------------------------- */

        dndrv6_mv_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
        for (j = 1; j <= n; ++j)
        {
            i__2 = j - 1;
            i__3 = ipntr[1] + j - 2;
            z__1.r = workd[i__3], z__1.i = 0.0;
            ctemp[i__2].r = z__1.r, ctemp[i__2].i = z__1.i;

        }

        zgttrs_("N", &n, &c__1, cdl, cdd, cdu, cdu2, ipiv, ctemp, &c__256, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" ERROR with _gttrs in _NDRV6.\n");
            printf(" \n");
            goto EXIT;
        }
        for (j = 1; j <= n; ++j)
        {
            workd[ipntr[1] + j - 2] = ctemp[j - 1].i;

        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call DNAUPD again. */
        /* --------------------------------------- */

        goto L20;
    }
    else if (ido == 1)
    {
        /* ---------------------------------------------------------- */
        /*                          Perform                           */
        /* y <--- OP*x = Imaginary_Part{inv[A-(SIGMAR,SIGMAI)*M]*M*x} */
        /* M*x has been saved in workd(ipntr(3)). The user only need  */
        /* the complex linear system solver here that takes           */
        /* complex[workd(ipntr(3))] as input, and returns the result  */
        /* to workd(ipntr(2)).                                        */
        /* ---------------------------------------------------------- */

        for (j = 1; j <= n; ++j)
        {
            i__2 = j - 1;
            i__3 = ipntr[2] + j - 2;
            z__1.r = workd[i__3], z__1.i = 0.0;
            ctemp[i__2].r = z__1.r, ctemp[i__2].i = z__1.i;

        }
        zgttrs_("N", &n, &c__1, cdl, cdd, cdu, cdu2, ipiv, ctemp, &c__256, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" ERROR with _gttrs in _NDRV6.\n");
            printf(" \n");
            goto EXIT;
        }
        for (j = 1; j <= n; ++j)
        {
            workd[ipntr[1] + j - 2] = ctemp[j - 1].i;

        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call DNAUPD again. */
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

        dndrv6_mv_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call DNAUPD again. */
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
        /* documentation in DNAUPD. */
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
    /* Post-Process using DNEUPD.                */
    /*                                           */
    /* Computed eigenvalues may be extracted.    */
    /*                                           */
    /* Eigenvectors may also be computed now if  */
    /* desired.  (indicated by rvec = .true.)    */
    /* ----------------------------------------- */

    rvec = true;
    dneupd_(&rvec, "A", select, d, &d[25], v, &n, &sigmar, &sigmai, workev, bmat, &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &ierr);

    /* --------------------------------------------- */
    /* The real part of the eigenvalue is returned   */
    /* in the first column of the two dimensional    */
    /* array D, and the IMAGINARY part is returned   */
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
        /* Check the documentation of DNEUPD. */
        /* ---------------------------------- */

        printf(" \n");
        printf(" Error with _neupd info = %d\n", ierr);
        printf(" Check the documentation of _neupd. \n");
        printf(" \n");

        goto EXIT;
    }

    first = true;
    nconv = iparam[4];
    for (j = 1; j <= nconv; ++j)
    {
        int k = (j - 1) * n;

        /* ----------------------------------- */
        /* Use Rayleigh Quotient to recover    */
        /* eigenvalues of the original problem.*/
        /* ----------------------------------- */

        if (d[j + 24] == 0.0)
        {
            /* -------------------------- */
            /* Eigenvalue is real.        */
            /* Compute d = x'(Ax)/x'(Mx). */
            /* -------------------------- */

            dndrv6_av_(n, &v[k], ax);
            numr = ddot_(&n, &v[k], &c__1, ax, &c__1);
            dndrv6_mv_(n, &v[k], ax);
            denr = ddot_(&n, &v[k], &c__1, ax, &c__1);
            d[j - 1] = numr / denr;
        }
        else if (first)
        {
            /* ---------------------- */
            /* Eigenvalue is complex. */
            /* Compute the first one  */
            /* of the conjugate pair. */
            /* ---------------------- */

            /* -------------- */
            /* Compute x'(Ax) */
            /* -------------- */

            dndrv6_av_(n, &v[k], ax);
            numr = ddot_(&n, &v[k], &c__1, ax, &c__1);
            numi = ddot_(&n, &v[j * n], &c__1, ax, &c__1);
            dndrv6_av_(n, &v[j * n], ax);
            numr += ddot_(&n, &v[j * n], &c__1, ax, &c__1);
            numi = -numi + ddot_(&n, &v[k], &c__1, ax, &c__1);

            /* -------------- */
            /* Compute x'(Mx) */
            /* -------------- */

            dndrv6_mv_(n, &v[k], ax);
            denr = ddot_(&n, &v[k], &c__1, ax, &c__1);
            deni = ddot_(&n, &v[j * n], &c__1, ax, &c__1);
            dndrv6_mv_(n, &v[j * n], ax);
            denr += ddot_(&n, &v[j * n], &c__1, ax, &c__1);
            deni = -deni + ddot_(&n, &v[k], &c__1, ax, &c__1);

            /* -------------- */
            /* d=x'(Ax)/x'(Mx)*/
            /* -------------- */

            d[j - 1] = (numr * denr + numi * deni) / dlapy2_(&denr, &deni);
            d[j + 24] = (numi * denr - numr * deni) / dlapy2_(&denr,&deni);
            first = false;
        }
        else
        {
            /* ---------------------------- */
            /* Get the second eigenvalue of */
            /* the conjugate pair by taking */
            /* the conjugate of the last    */
            /* eigenvalue computed.         */
            /* ---------------------------- */

            d[j - 1] = d[j - 2];
            d[j + 24] = -d[j + 23];
            first = true;
        }
    }

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

    first = true;
    nconv = iparam[4];
    for (j = 1; j <= nconv; ++j)
    {
        int k = (j - 1) * n;

        if (d[j + 24] == 0.0)
        {
            /* ------------------ */
            /* Ritz value is real */
            /* ------------------ */

            dndrv6_av_(n, &v[k], ax);
            dndrv6_mv_(n, &v[k], mx);
            d__1 = -d[j - 1];
            daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
            d[j + 49] = dnrm2_(&n, ax, &c__1);
            d[j + 49] /= (d__1 = d[j - 1], abs(d__1));
        }
        else if (first)
        {
            /* ---------------------- */
            /* Ritz value is complex  */
            /* Residual of one Ritz   */
            /* value of the conjugate */
            /* pair is computed.      */
            /* ---------------------- */

            dndrv6_av_(n, &v[k], ax);
            dndrv6_mv_(n, &v[k], mx);
            d__1 = -d[j - 1];
            daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
            dndrv6_mv_(n, &v[j * n], mx);
            daxpy_(&n, &d[j + 24], mx, &c__1, ax, &c__1);
            d[j + 49] = dnrm2_(&n, ax, &c__1);
            dndrv6_av_(n, &v[j * n], ax);
            dndrv6_mv_(n, &v[j * n], mx);
            d__1 = -d[j - 1];
            daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
            dndrv6_mv_(n, &v[k], mx);
            d__1 = -d[j + 24];
            daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
            d__1 = dnrm2_(&n, ax, &c__1);
            d[j + 49] = dlapy2_(&d[j + 49], &d__1);
            d[j + 49] /= dlapy2_(&d[j - 1], &d[j + 24]);
            d[j + 50] = d[j + 49];
            first = false;
        }
        else
        {
            first = true;
        }
    }

    /* --------------------------- */
    /* Display computed residuals. */
    /* --------------------------- */

    dmout_(&nconv, &c__3, d, &c__25, &c_n6, "Ritz values (Real,Imag) and relative residuals");

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
    printf(" _NDRV6 \n");
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

    free(cdd);
    free(cdl);
    free(cdu);
    free(cdu2);
    free(ctemp);

    free(ax);
    free(mx);
    free(resid);
    free(v);
    free(workl);
    free(workd);

    /* ------------------------- */
    /* Done with program dndrv6. */
    /* ------------------------- */

    return ierr;
}

/* ========================================================================== */

/*     matrix vector multiplication subroutine */

int dndrv6_mv_(const int n, double *v, double *w)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    int j;

    /*     Compute the matrix vector multiplication y<---M*x */
    /*     where M is a n by n symmetric tridiagonal matrix with 4 on the */
    /*     diagonal, 1 on the subdiagonal and superdiagonal. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    w[1] = v[1] * 4.0 + v[2] * 1.0;
    i__1 = n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = v[j - 1] * 1. + v[j] * 4.0 + v[j + 1] * 1.0;
    }
    w[n] = v[n - 1] * 1. + v[n] * 4.0;
    return 0;
} /* mv_ */

/* ------------------------------------------------------------------ */
int dndrv6_av_(const int n, double *v, double *w)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    int j;

    /*     Compute the matrix vector multiplication y<---A*x */
    /*     where M is a n by n symmetric tridiagonal matrix with 2 on the */
    /*     diagonal, -2 on the subdiagonal and 3 on the superdiagonal. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    w[1] = v[1] * 2. + v[2] * 3.0;
    i__1 = n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = v[j - 1] * -2. + v[j] * 2. + v[j + 1] * 3.0;
    }
    w[n] = v[n - 1] * -2. + v[n] * 2.0;
    return 0;
} /* av_ */

