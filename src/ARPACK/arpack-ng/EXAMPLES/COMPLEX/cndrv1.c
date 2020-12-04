/* EXAMPLES\COMPLEX\cndrv1.f -- translated by f2c (version 20100827). */

#include "arpack.h"

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
int cndrv1()
{
    /* System generated locals */
    int32_t i__1, i__2;
    complex q__1;

    double r_imag(complex *);

    /* Local variables */
    complex d[30];
    int32_t j, n;
    float rd[90]	/* was [30][3] */;
    complex ax[256];
    int32_t nx, ido, ncv, nev;
    float tol;
    char* bmat;
    int32_t mode, info;
    bool rvec;
    int32_t ierr;
    complex sigma;
    char* which;
    int32_t nconv;
    complex *v	/* was [256][30] */;
    complex *resid;
    complex *workd;
    complex *workl;
    int32_t ipntr[14];
    float rwork[30];
    int32_t iparam[11];
    bool select[30];
    int32_t ishfts, maxitr, lworkl;
    complex workev[90];

    resid = (complex*)malloc(256 * sizeof(complex));
    v = (complex*)malloc(7680 * sizeof(complex));
    workl = (complex*)malloc(2850 * sizeof(complex));
    workd = (complex*)malloc(768 * sizeof(complex));

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

    nx = 10;
    n = nx * nx;
    nev = 4;
    ncv = 20;
    if (n > 256) {
	printf(" ERROR with _NDRV1: N is greater than MAXN \n");
	return 0;
    } else if (nev > 12) {
	printf(" ERROR with _NDRV1: NEV is greater than MAXNEV \n");
	return 0;
    } else if (ncv > 30) {
	printf(" ERROR with _NDRV1: NCV is greater than MAXNCV \n");
	return 0;
    }
    bmat = "I";
    which = "LM";

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
    i__1 = ncv;
    lworkl = i__1 * i__1 * 3 + ncv * 5;
    tol = 0.f;
    ido = 0;
    info = 0;

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

    cnaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, rwork, &info);

    if (ido == -1 || ido == 1) {

           /* ----------------------------------------- */
           /* Perform matrix vector multiplication      */
           /*                y <--- OP*x                */
           /* The user should supply his/her own        */
           /* matrix vector multiplication routine here */
           /* that takes workd(ipntr(1)) as the input   */
           /* vector, and return the matrix vector      */
           /* product to workd(ipntr(2)).               */
           /* ----------------------------------------- */

	cndrv1_av_(&nx, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

           /* --------------------------------------- */
           /* L O O P   B A C K to call CNAUPD again. */
           /* --------------------------------------- */

	goto L10;
    }

     /* -------------------------------------- */
     /* Either we have convergence or there is */
     /* an error.                              */
     /* -------------------------------------- */

    if (info < 0) {

        /* ------------------------ */
        /* Error message, check the */
        /* documentation in CNAUPD  */
        /* ------------------------ */

	printf(" \n");
	printf(" Error with _naupd info = %d\n", info);
	printf(" Check the documentation of _naupd\n");
	printf(" \n");

    } else {

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

	if (ierr != 0) {

           /* ---------------------------------- */
           /* Error condition:                   */
           /* Check the documentation of CNEUPD. */
           /* ---------------------------------- */

	    printf(" \n");
	    printf(" Error with _neupd info = %d\n", ierr);
	    printf(" Check the documentation of _neupd. \n");
	    printf(" \n");

	} else {

	    nconv = iparam[4];
	    i__1 = nconv;
	    for (j = 1; j <= i__1; ++j) {

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

		cndrv1_av_(&nx, &v[(j << 8) - 256], ax);
		i__2 = j - 1;
		q__1.r = -d[i__2].r, q__1.i = -d[i__2].i;
		caxpy_(&n, &q__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
		i__2 = j - 1;
		rd[j - 1] = d[i__2].r;
		rd[j + 29] = r_imag(&d[j - 1]);
		rd[j + 59] = scnrm2_(&n, ax, &c__1);
		rd[j + 59] /= slapy2_(&rd[j - 1], &rd[j + 29]);
/* L20: */
	    }

            /* --------------------------- */
            /* Display computed residuals. */
            /* --------------------------- */

	    smout_(&nconv, &c__3, rd, &c__30, &c_n6, "Ritz values (Real, Imag) and relative residuals");
	}

        /* ----------------------------------------- */
        /* Print additional convergence information. */
        /* ----------------------------------------- */

	if (info == 1) {
	    printf(" \n");
	    printf(" Maximum number of iterations reached.\n");
	    printf(" \n");
	} else if (info == 3) {
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

    }

    free(resid);
    free(v);
    free(workl);
    free(workd);

     /* ------------------------- */
     /* Done with program cndrv1. */
     /* ------------------------- */

    return 0;
} /* MAIN__ */

/* ========================================================================== */

/*     matrix vector subroutine */

/*     The matrix used is the convection-diffusion operator */
/*     discretized using centered difference. */

int cndrv1_av_(int32_t *nx, complex *v, complex *w)
{
    /* System generated locals */
    int32_t i__1;
    complex q__1, q__2;

    /* Builtin functions */
    void c_div(complex *, complex *, complex *);

    /* Local variables */
    int32_t j;
    complex h2;
    int32_t lo;

/*     Computes w <--- OP*v, where OP is the nx*nx by nx*nx block */
/*     tridiagonal matrix */

/*                  | T -I          | */
/*                  |-I  T -I       | */
/*             OP = |   -I  T       | */
/*                  |        ...  -I| */
/*                  |           -I T| */

/*     derived from the standard central difference  discretization */
/*     of the convection-diffusion operator (Laplacian u) + rho*(du/dx) */
/*     with zero boundary condition. */

/*     The subroutine TV is called to computed y<---T*x. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    i__1 = (*nx + 1) * (*nx + 1);
    q__2.r = (float) i__1, q__2.i = 0.f;
    c_div(&q__1, &c_one, &q__2);
    h2.r = q__1.r, h2.i = q__1.i;

	cndrv1_tv_(nx, &v[1], &w[1]);
    q__2.r = -1.f, q__2.i = -0.f;
    c_div(&q__1, &q__2, &h2);
    caxpy_(nx, &q__1, &v[*nx + 1], &c__1, &w[1], &c__1);

    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j) {
	lo = (j - 1) * *nx;
	cndrv1_tv_(nx, &v[lo + 1], &w[lo + 1]);
	q__2.r = -1.f, q__2.i = -0.f;
	c_div(&q__1, &q__2, &h2);
	caxpy_(nx, &q__1, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);
	q__2.r = -1.f, q__2.i = -0.f;
	c_div(&q__1, &q__2, &h2);
	caxpy_(nx, &q__1, &v[lo + *nx + 1], &c__1, &w[lo + 1], &c__1);
/* L10: */
    }

    lo = (*nx - 1) * *nx;
	cndrv1_tv_(nx, &v[lo + 1], &w[lo + 1]);
    q__2.r = -1.f, q__2.i = -0.f;
    c_div(&q__1, &q__2, &h2);
    caxpy_(nx, &q__1, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);

    return 0;
} /* av_ */

/* ========================================================================= */
int cndrv1_tv_(int32_t *nx, complex *x, complex *y)
{
    /* System generated locals */
    int32_t i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2, q__3, q__4, q__5;

    /* Builtin functions */
    void c_div(complex *, complex *, complex *);

    /* Local variables */
    complex h;
    int32_t j;
    complex h2, dd, dl, du;

/*     Compute the matrix vector multiplication y<---T*x */
/*     where T is a nx by nx tridiagonal matrix with DD on the */
/*     diagonal, DL on the subdiagonal, and DU on the superdiagonal */

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    i__1 = *nx + 1;
    q__2.r = (float) i__1, q__2.i = 0.f;
    c_div(&q__1, &c_one, &q__2);
    h.r = q__1.r, h.i = q__1.i;
    q__1.r = h.r * h.r - h.i * h.i, q__1.i = h.r * h.i + h.i * 
	    h.r;
    h2.r = q__1.r, h2.i = q__1.i;
    c_div(&q__1, &c_b151, &h2);
    dd.r = q__1.r, dd.i = q__1.i;
    q__3.r = -1.f, q__3.i = -0.f;
    c_div(&q__2, &q__3, &h2);
    q__5.r = 50.f, q__5.i = 0.f;
    c_div(&q__4, &q__5, &h);
    q__1.r = q__2.r - q__4.r, q__1.i = q__2.i - q__4.i;
    dl.r = q__1.r, dl.i = q__1.i;
    q__3.r = -1.f, q__3.i = -0.f;
    c_div(&q__2, &q__3, &h2);
    q__5.r = 50.f, q__5.i = 0.f;
    c_div(&q__4, &q__5, &h);
    q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
    du.r = q__1.r, du.i = q__1.i;

    q__2.r = dd.r * x[1].r - dd.i * x[1].i, q__2.i = dd.r * x[1].i + dd.i * x[
	    1].r;
    q__3.r = du.r * x[2].r - du.i * x[2].i, q__3.i = du.r * x[2].i + du.i * x[
	    2].r;
    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
    y[1].r = q__1.r, y[1].i = q__1.i;
    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j - 1;
	q__3.r = dl.r * x[i__3].r - dl.i * x[i__3].i, q__3.i = dl.r * x[i__3]
		.i + dl.i * x[i__3].r;
	i__4 = j;
	q__4.r = dd.r * x[i__4].r - dd.i * x[i__4].i, q__4.i = dd.r * x[i__4]
		.i + dd.i * x[i__4].r;
	q__2.r = q__3.r + q__4.r, q__2.i = q__3.i + q__4.i;
	i__5 = j + 1;
	q__5.r = du.r * x[i__5].r - du.i * x[i__5].i, q__5.i = du.r * x[i__5]
		.i + du.i * x[i__5].r;
	q__1.r = q__2.r + q__5.r, q__1.i = q__2.i + q__5.i;
	y[i__2].r = q__1.r, y[i__2].i = q__1.i;
/* L10: */
    }
    i__1 = *nx;
    i__2 = *nx - 1;
    q__2.r = dl.r * x[i__2].r - dl.i * x[i__2].i, q__2.i = dl.r * x[i__2].i + 
	    dl.i * x[i__2].r;
    i__3 = *nx;
    q__3.r = dd.r * x[i__3].r - dd.i * x[i__3].i, q__3.i = dd.r * x[i__3].i + 
	    dd.i * x[i__3].r;
    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
    y[i__1].r = q__1.r, y[i__1].i = q__1.i;
    return 0;
} /* tv_ */

