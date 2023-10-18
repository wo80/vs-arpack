/* EXAMPLES\NONSYM\sndrv5.f -- translated by f2c (version 20100827). */

#include <stdlib.h>
#include <stdio.h>
#include "arpack.h"
#include "lapack.h"

int sndrv5_av_(const int n, float* v, float* w);
int sndrv5_mv_(const int n, float* v, float* w);

extern int smout_(const int, const int, const float*, const int, const int, const char*);

static int i_one = 1;

/**
 * \BeginDoc
 *
 *     Simple program to illustrate the idea of reverse communication
 *     in shift-invert mode for a generalized nonsymmetric eigenvalue problem.
 *
 *     We implement example five of ex-nonsym.doc in DOCUMENTS directory
 *
 * \Example-5
 *
 *     ... Suppose we want to solve A*x = lambda*B*x in shift-invert mode
 *         The matrix A is the tridiagonal matrix with 2 on the diagonal,
 *         -2 on the subdiagonal and 3 on the superdiagonal.  The matrix M
 *         is the tridiagonal matrix with 4 on the diagonal and 1 on the
 *         off-diagonals.
 *     ... The shift sigma is a complex number (sigmar, sigmai).
 *     ... OP = Real_Part{inv[A-(SIGMAR,SIGMAI)*M]*M and  B = M.
 *     ... Use mode 3 of SNAUPD.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Routines called:
 *     snaupd  ARPACK reverse communication interface routine.
 *     sneupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     cgttrf  LAPACK complex matrix factorization routine.
 *     cgttrs  LAPACK complex linear system solve routine.
 *     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
 *     sdot    Level 1 BLAS that computes the dot product of two vectors.
 *     snrm2   Level 1 BLAS that computes the norm of a vector
 *     av      Matrix vector subroutine that computes A*x.
 *     mv      Matrix vector subroutine that computes M*x.
 *
 * \EndLib
 */
int main()
{
    /* System generated locals */
    int i__1, i__2, i__3;
    float r__1, r__2;
    a_fcomplex q__1;

    /* Local variables */
    float d[75]; /* (3 * MAXNCV) */
    float workev[75];

    float denr, deni;
    float numr, numi;
    float sigmar, sigmai;

    a_fcomplex c1, c2, c3;

    int j;
    int ierr, nconv;
    int ishfts, maxitr, mode;

    int ipiv[256];
    int ipntr[14];
    int iparam[11];
    logical select[25];
    logical first, rvec;

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
    /* to 'LM'.  The user can modify NEV, NCV, SIGMAR,    */
    /* SIGMAI to solve problems of different sizes, and   */
    /* to get different parts of the spectrum. However,   */
    /* The following conditions must be satisfied:        */
    /*                     N <= MAXN,                     */
    /*                   NEV <= MAXNEV,                   */
    /*               NEV + 2 <= NCV <= MAXNCV             */
    /* -------------------------------------------------- */

    int n = 100;
    int nev = 4;
    int ncv = 20;
    if (n > MAXN)
    {
        printf(" ERROR with _NDRV5: N is greater than MAXN \n");
        return 0;
    }
    else if (nev > MAXNEV)
    {
        printf(" ERROR with _NDRV5: NEV is greater than MAXNEV \n");
        return 0;
    }
    else if (ncv > MAXNCV)
    {
        printf(" ERROR with _NDRV5: NCV is greater than MAXNCV \n");
        return 0;
    }
    char* bmat = "G";
    char* which = "LM";
    sigmar = 0.4f;
    sigmai = 0.6f;

    /* ------------------------------------------------- */
    /* Construct C = A - (SIGMAR,SIGMAI)*M in complex    */
    /* arithmetic, and factor C in complex arithmetic    */
    /* (using LAPACK subroutine cgttrf). The matrix A is */
    /* chosen to be the tridiagonal matrix with -2 on    */
    /* the subdiagonal, 2 on the diagonal and 3 on the   */
    /* superdiagonal. The matrix M is chosen to the      */
    /* symmetric tridiagonal matrix with 4 on the        */
    /* diagonal and 1 on the off-diagonals.              */
    /* ------------------------------------------------- */

    a_fcomplex* cdd = (a_fcomplex*)malloc(n * sizeof(a_fcomplex));
    a_fcomplex* cdl = (a_fcomplex*)malloc(n * sizeof(a_fcomplex));
    a_fcomplex* cdu = (a_fcomplex*)malloc(n * sizeof(a_fcomplex));
    a_fcomplex* cdu2 = (a_fcomplex*)malloc(n * sizeof(a_fcomplex));
    a_fcomplex* ctemp = (a_fcomplex*)malloc(n * sizeof(a_fcomplex));

    r__1 = -2.0f - sigmar;
    r__2 = -sigmai;
    q__1.r = r__1, q__1.i = r__2;
    c1.r = q__1.r, c1.i = q__1.i;
    r__1 = 2.0f - sigmar * 4.0f;
    r__2 = sigmai * -4.0f;
    q__1.r = r__1, q__1.i = r__2;
    c2.r = q__1.r, c2.i = q__1.i;
    r__1 = 3.0f - sigmar;
    r__2 = -sigmai;
    q__1.r = r__1, q__1.i = r__2;
    c3.r = q__1.r, c3.i = q__1.i;

    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        i__2 = j - 1;
        cdl[i__2].r = c1.r, cdl[i__2].i = c1.i;
        cdd[i__2].r = c2.r, cdd[i__2].i = c2.i;
        cdu[i__2].r = c3.r, cdu[i__2].i = c3.i;
    }
    cdd[i__1].r = c2.r, cdd[i__1].i = c2.i;

    cgttrf_(&n, cdl, cdd, cdu, cdu2, ipiv, &ierr);
    if (ierr != 0)
    {
        printf(" \n");
        printf(" ERROR with _gttrf in _NDRV5.\n");
        printf(" \n");
        return 0;
    }

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
    float* mx = (float*)malloc(n * sizeof(float));
    float* resid = (float*)malloc(n * sizeof(float));
    float* v = (float*)malloc(n * ncv * sizeof(float));
    float* workl = (float*)malloc(lworkl * sizeof(float));
    float* workd = (float*)malloc(3 * n * sizeof(float));

    /* ------------------------------------------------- */
    /* This program uses exact shift with respect to     */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 3 of SNAUPD is used     */
    /* (IPARAM(7) = 3).  All these options can be        */
    /* changed by the user. For details, see the         */
    /* documentation in SNAUPD.                          */
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
    /* Repeatedly call the routine SNAUPD and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    snaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);

    if (ido == -1)
    {
        /* ----------------------------------------------------- */
        /*                       Perform                         */
        /* y <--- OP*x = Real_Part{inv[A-(SIGMAR,SIGMAI)*M]*M*x} */
        /* to force starting vector into the range of OP. The    */
        /* user should supply his/her own matrix vector          */
        /* multiplication routine and a complex linear system    */
        /* solver.  The matrix vector multiplication routine     */
        /* should take workd(ipntr(1)) as the input. The final   */
        /* result (a real vector) should be returned to          */
        /* workd(ipntr(2)).                                      */
        /* ----------------------------------------------------- */

        sndrv5_mv_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
        for (j = 1; j <= n; ++j)
        {
            i__2 = j - 1;
            i__3 = ipntr[1] + j - 2;
            q__1.r = workd[i__3], q__1.i = 0.0f;
            ctemp[i__2].r = q__1.r, ctemp[i__2].i = q__1.i;

        }

        cgttrs_("N", &n, &i_one, cdl, cdd, cdu, cdu2, ipiv, ctemp, &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" ERROR with _gttrs in _NDRV5.\n");
            printf(" \n");
            goto EXIT;
        }
        for (j = 1; j <= n; ++j)
        {
            i__2 = j - 1;
            workd[ipntr[1] + j - 2] = ctemp[i__2].r;

        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call SNAUPD again. */
        /* --------------------------------------- */

        goto L20;
    }
    else if (ido == 1)
    {
        /* ----------------------------------------------------- */
        /*                        Perform                        */
        /* y <--- OP*x = Real_Part{inv[A-(SIGMAR,SIGMAI)*M]*M*x} */
        /* M*x has been saved in workd(ipntr(3)). The user only  */
        /* needs the complex linear system solver here that      */
        /* takes complex[workd(ipntr(3))] as input, and returns  */
        /* the result to workd(ipntr(2)).                        */
        /* ----------------------------------------------------- */

        for (j = 1; j <= n; ++j)
        {
            i__2 = j - 1;
            i__3 = ipntr[2] + j - 2;
            q__1.r = workd[i__3], q__1.i = 0.0f;
            ctemp[i__2].r = q__1.r, ctemp[i__2].i = q__1.i;

        }
        cgttrs_("N", &n, &i_one, cdl, cdd, cdu, cdu2, ipiv, ctemp, &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" ERROR with _gttrs in _NDRV5.\n");
            printf(" \n");
            goto EXIT;
        }
        for (j = 1; j <= n; ++j)
        {
            i__2 = j - 1;
            workd[ipntr[1] + j - 2] = ctemp[i__2].r;

        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call SNAUPD again. */
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

        sndrv5_mv_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call SNAUPD again. */
        /* --------------------------------------- */

        goto L20;
    }

    /* ---------------------------------------- */
    /* Either we have convergence, or there is  */
    /* an error.                                */
    /* ---------------------------------------- */

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
    sneupd_(&rvec, "A", select, d, &d[MAXNCV], v, &n, &sigmar, &sigmai, workev, bmat, &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &ierr);

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
        /* Check the documentation of SNEUPD. */
        /* ---------------------------------- */

        printf(" \n");
        printf(" Error with _neupd = \n");
        printf("%d", ierr);
        printf(" Check the documentation of _neupd. \n");
        printf(" \n");

        goto EXIT;
    }

    first = TRUE_;
    nconv = iparam[4];
    for (j = 1; j <= nconv; ++j)
    {
        int k = (j - 1) * n;

        /* ----------------------------------- */
        /* Use Rayleigh Quotient to recover    */
        /* eigenvalues of the original problem.*/
        /* ----------------------------------- */

        if (d[j + 24] == 0.0f)
        {
            /* ------------------------- */
            /*    Eigenvalue is real.    */
            /* Compute d = x'(Ax)/x'(Mx).*/
            /* ------------------------- */

            sndrv5_av_(n, &v[k], ax);
            numr = sdot_(&n, &v[k], &i_one, ax, &i_one);
            sndrv5_mv_(n, &v[k], ax);
            denr = sdot_(&n, &v[k], &i_one, ax, &i_one);
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
            sndrv5_av_(n, &v[k], ax);
            numr = sdot_(&n, &v[k], &i_one, ax, &i_one);
            numi = sdot_(&n, &v[j * n], &i_one, ax, &i_one);
            sndrv5_av_(n, &v[j * n], ax);
            numr += sdot_(&n, &v[j * n], &i_one, ax, &i_one);
            numi = -numi + sdot_(&n, &v[k], &i_one, ax, &i_one);

            /* -------------- */
            /* Compute x'(Mx) */
            /* -------------- */

            sndrv5_mv_(n, &v[k], ax);
            denr = sdot_(&n, &v[k], &i_one, ax, &i_one);
            deni = sdot_(&n, &v[j * n], &i_one, ax, &i_one);
            sndrv5_mv_(n, &v[j * n], ax);
            denr += sdot_(&n, &v[j * n], &i_one, ax, &i_one);
            deni = -deni + sdot_(&n, &v[k], &i_one, ax, &i_one);

            /* -------------- */
            /* d=x'(Ax)/x'(Mx)*/
            /* -------------- */

            d[j - 1] = (numr * denr + numi * deni) / slapy2_(&denr, &deni);
            d[j + 24] = (numi * denr - numr * deni) / slapy2_(&denr,&deni);
            first = FALSE_;
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
            first = TRUE_;
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

    first = TRUE_;
    for (j = 1; j <= nconv; ++j)
    {
        int k = (j - 1) * n;

        if (d[j + 24] == 0.0f)
        {
            /* ------------------ */
            /* Ritz value is real */
            /* ------------------ */

            sndrv5_av_(n, &v[k], ax);
            sndrv5_mv_(n, &v[k], mx);
            r__1 = -d[j - 1];
            saxpy_(&n, &r__1, mx, &i_one, ax, &i_one);
            d[j + 49] = snrm2_(&n, ax, &i_one);
            d[j + 49] /= (r__1 = d[j - 1], dabs(r__1));
        }
        else if (first)
        {
            /* ---------------------- */
            /* Ritz value is complex  */
            /* Residual of one Ritz   */
            /* value of the conjugate */
            /* pair is computed.      */
            /* ---------------------- */

            sndrv5_av_(n, &v[k], ax);
            sndrv5_mv_(n, &v[k], mx);
            r__1 = -d[j - 1];
            saxpy_(&n, &r__1, mx, &i_one, ax, &i_one);
            sndrv5_mv_(n, &v[j * n], mx);
            saxpy_(&n, &d[j + 24], mx, &i_one, ax, &i_one);
            d[j + 49] = snrm2_(&n, ax, &i_one);
            sndrv5_av_(n, &v[j * n], ax);
            sndrv5_mv_(n, &v[j * n], mx);
            r__1 = -d[j - 1];
            saxpy_(&n, &r__1, mx, &i_one, ax, &i_one);
            sndrv5_mv_(n, &v[k], mx);
            r__1 = -d[j + 24];
            saxpy_(&n, &r__1, mx, &i_one, ax, &i_one);
            r__1 = snrm2_(&n, ax, &i_one);
            d[j + 49] = slapy2_(&d[j + 49], &r__1);
            d[j + 49] /= slapy2_(&d[j - 1], &d[j + 24]);
            d[j + 50] = d[j + 49];
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

    smout_(nconv, 3, d, MAXNCV, -6, "Ritz values (Real,Imag) and relative residuals");

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
    printf(" _NDRV5 \n");
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
    /* Done with program sndrv5. */
    /* ------------------------- */

    return ierr;
}

/* ========================================================================== */

/*     matrix vector multiplication subroutine */

int sndrv5_mv_(const int n, float *v, float *w)
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
    w[1] = v[1] * 4.0f + v[2] * 1.0f;
    i__1 = n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = v[j - 1] * 1.0f + v[j] * 4.0f + v[j + 1] * 1.0f;
    }
    w[n] = v[n - 1] * 1.0f + v[n] * 4.0f;
    return 0;
} /* mv_ */

/* ------------------------------------------------------------------ */
int sndrv5_av_(const int n, float *v, float *w)
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
    w[1] = v[1] * 2.0f + v[2] * 3.0f;
    i__1 = n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = v[j - 1] * -2.0f + v[j] * 2.0f + v[j + 1] * 3.0f;
    }
    w[n] = v[n - 1] * -2.0f + v[n] * 2.0f;
    return 0;
} /* av_ */

