/* EXAMPLES\NONSYM\sndrv2.f -- translated by f2c (version 20100827). */

#include <stdlib.h>
#include "arpack.h"

struct
{
    float rho;
} convct_;

#define convct_1 convct_

/**
 * \BeginDoc
 *
 *     Simple program to illustrate the idea of reverse communication
 *     in shift-invert mode for a standard nonsymmetric eigenvalue problem.
 *
 *     We implement example two of ex-nonsym.doc in DOCUMENTS directory
 *
 * \Example-2
 *     ... Suppose we want to solve A*x = lambda*x in shift-invert mode,
 *         where A is derived from the centered difference discretization
 *         of the 1-dimensional convection-diffusion operator
 *                          (d^2u / dx^2) + rho*(du/dx)
 *         on the interval [0,1] with zero Dirichlet boundary condition.
 *
 *     ... The shift sigma is a real number.
 *
 *     ... OP = inv[A-sigma*I] and  B = I.
 *
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
 *     sgttrf  LAPACK tridiagonal factorization routine.
 *     sgttrs  LAPACK tridiagonal solve routine.
 *     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
 *     scopy   Level 1 BLAS that copies one vector to another.
 *     sdot    Level 1 BLAS that computes the dot product of two vectors.
 *     snrm2   Level 1 BLAS that computes the norm of a vector.
 *     av      Matrix vector multiplication routine that computes A*x.
 *
 * \EndLib
 */
int sndrv2()
{
    /* System generated locals */
    int i__1;
    float r__1;

    /* Local variables */
    float d[75]	/* was [25][3] */, h;
    int j, n;
    float s, s1, s2, s3, dd[256], dl[256];
    float ax[256], du[256], du2[256];
    int ido, ncv, nev;
    float tol;
    char* bmat;
    int mode, info;
    bool rvec;
    int ierr, ipiv[256];
    char* which;
    int nconv;
    float *v	/* was [256][25] */;
    float *resid;
    float *workd;
    float *workl;
    bool first;
    int ipntr[14];
    int iparam[11];
    float sigmai;
    bool select[25];
    float sigmar;
    int ishfts, maxitr;
    int lworkl;
    float workev[75];

    resid = (float*)malloc(256 * sizeof(float));
    v = (float*)malloc(6400 * sizeof(float));
    workl = (float*)malloc(2025 * sizeof(float));
    workd = (float*)malloc(768 * sizeof(float));

    /* Define maximum dimensions for all arrays. */

    const int MAXN   = 256; /* Maximum dimension of the A allowed. */
    const int MAXNEV =  10; /* Maximum NEV allowed */
    const int MAXNCV =  25; /* Maximum NCV allowed */

    /* ------------------------------------------------ */
    /* The number N is the dimension of the matrix.  A  */
    /* standard eigenvalue problem is solved (BMAT =    */
    /* 'I').  NEV is the number of eigenvalues (closest */
    /* to the shift SIGMAR) to be approximated.  Since  */
    /* the shift-invert mode is used, WHICH is set to   */
    /* 'LM'.  The user can modify NEV, NCV, SIGMAR to   */
    /* solve problems of different sizes, and to get    */
    /* different parts of the spectrum.  However, The   */
    /* following conditions must be satisfied:          */
    /*                 N <= MAXN,                       */
    /*               NEV <= MAXNEV,                     */
    /*           NEV + 2 <= NCV <= MAXNCV               */
    /* ------------------------------------------------ */

    n = 100;
    nev = 4;
    ncv = 20;
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
    bmat = "I";
    which = "LM";
    sigmar = 1.f;
    sigmai = 0.f;

    /* -------------------------------------------------- */
    /* Construct C = A - SIGMA*I in real arithmetic, and  */
    /* factor C in real arithmetic using LAPACK           */
    /* subroutine sgttrf. The matrix A is chosen to be    */
    /* the tridiagonal matrix derived from standard       */
    /* central difference of the 1-d convection diffusion */
    /* operator u" + rho*u' on the interval [0, 1] with   */
    /* zero Dirichlet boundary condition.                 */
    /* -------------------------------------------------- */

    convct_1.rho = 10.f;
    h = 1.f / (float) (n + 1);
    s = convct_1.rho * h / 2.f;

    s1 = -1.f - s;
    s2 = 2.f - sigmar;
    s3 = s - 1.f;

    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        dl[j - 1] = s1;
        dd[j - 1] = s2;
        du[j - 1] = s3;
    }
    dd[n - 1] = s2;

    sgttrf_(&n, dl, dd, du, du2, ipiv, &ierr);
    if (ierr != 0)
    {
        printf(" \n");
        printf(" ERROR with _gttrf in _NDRV2.\n");
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

    /* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 * 3 + ncv * 6;
    tol = 0.f;
    ido = 0;
    info = 0;

    /* ------------------------------------------------- */
    /* This program uses exact shifts with respect to    */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed. Mode 3 of SNAUPD is used      */
    /* (IPARAM(7) = 3).  All these options can be        */
    /* changed by the user. For details see the          */
    /* documentation in SNAUPD.                          */
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
    /* Repeatedly call the routine SNAUPD and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    snaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &info);

    if (ido == -1 || ido == 1)
    {
        /* ----------------------------------------- */
        /* Perform  y <--- OP*x = inv[A-SIGMA*I]*x   */
        /* The user should supply his/her own linear */
        /* system solver here that takes             */
        /* workd(ipntr(1)) as the input, and returns */
        /* the result to workd(ipntr(2)).            */
        /* ----------------------------------------- */

        scopy_(&n, &workd[ipntr[0] - 1], &c__1, &workd[ipntr[1] - 1], &c__1);

        sgttrs_("N", &n, &c__1, dl, dd, du, du2, ipiv, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" ERROR with _gttrs in _NDRV2.\n");
            printf(" \n");
            return ierr;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call SNAUPD again. */
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
        /* documentation in SNAUPD. */
        /* ------------------------ */

        printf(" \n");
        printf(" Error with _naupd info = %d\n", info);
        printf(" Check the documentation in _naupd.\n");
        printf(" \n");
    }
    else
    {
        /* ----------------------------------------- */
        /* No fatal errors occurred.                 */
        /* Post-Process using SNEUPD.                */
        /*                                           */
        /* Computed eigenvalues may be extracted.    */
        /*                                           */
        /* Eigenvectors may also be computed now if  */
        /* desired.  (indicated by rvec = .true.)    */
        /* ----------------------------------------- */

        rvec = true;

        sneupd_(&rvec, "A", select, d, &d[25], v, &c__256, &sigmar, &sigmai, workev, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &ierr);

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

        }
        else
        {
            first = true;
            nconv = iparam[4];
            i__1 = nconv;
            for (j = 1; j <= i__1; ++j)
            {
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

                if (d[j + 24] == 0.f)
                {
                    /* ------------------ */
                    /* Ritz value is real */
                    /* ------------------ */

                    sndrv2_av_(&n, &v[(j << 8) - 256], ax);
                    r__1 = -d[j - 1];
                    saxpy_(&n, &r__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
                    d[j + 49] = snrm2_(&n, ax, &c__1);
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

                    sndrv2_av_(&n, &v[(j << 8) - 256], ax);
                    r__1 = -d[j - 1];
                    saxpy_(&n, &r__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
                    saxpy_(&n, &d[j + 24], &v[(j + 1 << 8) - 256], &c__1, ax, &c__1);
                    d[j + 49] = snrm2_(&n, ax, &c__1);
                    sndrv2_av_(&n, &v[(j + 1 << 8) - 256], ax);
                    r__1 = -d[j + 24];
                    saxpy_(&n, &r__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
                    r__1 = -d[j - 1];
                    saxpy_(&n, &r__1, &v[(j + 1 << 8) - 256], &c__1, ax, &c__1);
                    r__1 = snrm2_(&n, ax, &c__1);
                    d[j + 49] = slapy2_(&d[j + 49], &r__1);
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

            smout_(&nconv, &c__3, d, &c__25, &c_n6, "Ritz values (Real,Imag) and relative residuals");

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
        printf(" _NDRV2 \n");
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
    }

    free(resid);
    free(v);
    free(workl);
    free(workd);

    /* ------------------------- */
    /* Done with program sndrv2. */
    /* ------------------------- */

    return 0;
} /* MAIN__ */

/* ------------------------------------------------------------------- */

/*     matrix vector multiplication subroutine */

int sndrv2_av_(int *n, float *v, float *w)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    float h;
    int j;
    float s, dd, dl, du;

    /*     Compute the matrix vector multiplication y<---A*x */
    /*     where A is a n by n nonsymmetric tridiagonal matrix derived from */
    /*     the central difference discretization of the 1-dimensional */
    /*     convection diffusion operator on the interval [0,1] with */
    /*     zero Dirichlet boundary condition. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    h = 1.f / (float) (*n + 1);
    s = convct_1.rho * h / 2.f;
    dd = 2.f;
    dl = -1.f - s;
    du = s - 1.f;

    w[1] = dd * v[1] + du * v[2];
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = dl * v[j - 1] + dd * v[j] + du * v[j + 1];
    }
    w[*n] = dl * v[*n - 1] + dd * v[*n];
    return 0;
} /* av_ */

