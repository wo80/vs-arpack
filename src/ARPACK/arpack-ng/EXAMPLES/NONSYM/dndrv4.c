/* EXAMPLES\NONSYM\dndrv4.f -- translated by f2c (version 20100827). */

#include "arpack.h"

struct {
    double rho;
} convct_;

#define convct_1 convct_

int dndrv4()
{
    /* System generated locals */
    int32_t i__1;
    double d__1;

    /* Local variables */
    double d[75]	/* was [25][3] */, h;
    int32_t j, n;
    double s, v[6400]	/* was [256][25] */, s1, s2, s3, dd[256], dl[
	    256];
    double ax[256], du[256];
    double mx[256], du2[256];
    int32_t ido, ncv, nev;
    double tol;
    char* bmat;
    int32_t mode, info;
    bool rvec;
    int32_t ierr, ipiv[256];
    char* which;
    double resid[256];
    int32_t nconv;
    double workd[768];
    bool first;
    int32_t ipntr[14];
    double workl[2025];
    int32_t iparam[11];
    double sigmai;
    bool select[25];
    double sigmar;
    int32_t ishfts, maxitr;
    int32_t lworkl;
    double workev[75];

    /* Fortran I/O blocks */

/*     Simple program to illustrate the idea of reverse communication */
/*     in shift-invert mode for a generalized nonsymmetric eigenvalue */
/*     problem. */

/*     We implement example four of ex-nonsym.doc in DOCUMENTS directory */

/* \Example-4 */
/*     ... Suppose we want to solve A*x = lambda*B*x in inverse mode, */
/*         where A and B are derived from the finite element discretization */
/*         of the 1-dimensional convection-diffusion operator */
/*                           (d^2u / dx^2) + rho*(du/dx) */
/*         on the interval [0,1] with zero Dirichlet boundary condition */
/*         using linear elements. */

/*     ... The shift sigma is a real number. */

/*     ... OP = inv[A-SIGMA*M]*M  and  B = M. */

/*     ... Use mode 3 of DNAUPD. */
/**
 * \BeginLib
 *
 * \Routines called:
 *     dnaupd  ARPACK reverse communication interface routine.
 *     dneupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     dgttrf  LAPACK tridiagonal factorization routine.
 *     dgttrs  LAPACK tridiagonal linear system solve routine.
 *     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
 *     dcopy   Level 1 BLAS that copies one vector to another.
 *     ddot    Level 1 BLAS that computes the dot product of two vectors.
 *     dnrm2   Level 1 BLAS that computes the norm of a vector.
 *     av      Matrix vector multiplication routine that computes A*x.
 *     mv      Matrix vector multiplication routine that computes M*x.
 *
 * \Author
 *     Richard Lehoucq
 *     Danny Sorensen
 *     Chao Yang
 *     Dept. of Computational &
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * \SCCS Information: @(#)
 * FILE: ndrv4.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
 *
 * \Remarks
 *     1. None
 *
 * \EndLib
 */
     /* --------------------------- */
     /* Define leading dimensions   */
     /* for all arrays.             */
     /* MAXN:   Maximum dimension   */
     /*         of the A allowed.   */
     /* MAXNEV: Maximum NEV allowed */
     /* MAXNCV: Maximum NCV allowed */
     /* --------------------------- */

     /* --------------------- */
     /* Executable Statements */
     /* --------------------- */

     /* -------------------------------------------------- */
     /* The number N is the dimension of the matrix.  A    */
     /* generalized eigenvalue problem is solved (BMAT =   */
     /* 'G').  NEV is the number of eigenvalues (closest   */
     /* to SIGMAR) to be approximated.  Since the          */
     /* shift-invert mode is used,  WHICH is set to 'LM'.  */
     /* The user can modify NEV, NCV, SIGMAR to solve      */
     /* problems of different sizes, and to get different  */
     /* parts of the spectrum.  However, The following     */
     /* conditions must be satisfied:                      */
     /*                     N <= MAXN,                     */
     /*                   NEV <= MAXNEV,                   */
     /*               NEV + 2 <= NCV <= MAXNCV             */
     /* -------------------------------------------------- */

    n = 100;
    nev = 4;
    ncv = 10;
    if (n > 256) {
	printf(" ERROR with _NDRV4: N is greater than MAXN \n");
	return 0;
    } else if (nev > 10) {
	printf(" ERROR with _NDRV4: NEV is greater than MAXNEV \n");
	return 0;
    } else if (ncv > 25) {
	printf(" ERROR with _NDRV4: NCV is greater than MAXNCV \n");
	return 0;
    }
    bmat = "G";
    strcpy(which, "LM");
    sigmar = 1.;
    sigmai = 0.;

     /* ------------------------------------------------ */
     /* Construct C = A - SIGMA*M in real arithmetic,    */
     /* and factor C in real arithmetic (using LAPACK    */
     /* subroutine dgttrf). The matrix A is chosen to be */
     /* the tridiagonal matrix derived from the standard */
     /* central difference discretization of the 1-d     */
     /* convection-diffusion operator u" + rho*u' on the */
     /* interval [0, 1] with zero Dirichlet boundary     */
     /* condition.  The matrix M is the mass matrix      */
     /* formed by using piecewise linear elements on     */
     /* [0,1].                                           */
     /* ------------------------------------------------ */

    convct_1.rho = 10.;
    h = 1. / (double) (n + 1);
    s = convct_1.rho / 2.;

    s1 = -1. / h - s - sigmar * h / 6.;
    s2 = 2. / h - sigmar * 4. * h / 6.;
    s3 = -1. / h + s - sigmar * h / 6.;

    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j) {
	dl[j - 1] = s1;
	dd[j - 1] = s2;
	du[j - 1] = s3;
/* L10: */
    }
    dd[n - 1] = s2;

    dgttrf_(&n, dl, dd, du, du2, ipiv, &ierr);
    if (ierr != 0) {
	printf(" \n");
	printf(" ERROR with _gttrf in _NDRV4.\n");
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

/* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 * 3 + ncv * 6;
    tol = 0.;
    ido = 0;
    info = 0;

     /* ------------------------------------------------- */
     /* This program uses exact shifts with respect to    */
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

    dnaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, 
	    iparam, ipntr, workd, workl, &lworkl, &info, (ftnlen)1, (ftnlen)2)
	    ;

    if (ido == -1) {

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

	dndrv4_mv_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
	dgttrs_("N", &n, &c__1, dl, dd, du, du2, ipiv, &workd[ipntr[1] - 1], &
		n, &ierr, (ftnlen)1);
	if (ierr != 0) {
	    printf(" \n");
	    printf(" ERROR with _gttrs in _NDRV4.\n");
	    printf(" \n");
	    return ierr;
	}

           /* --------------------------------------- */
           /* L O O P   B A C K to call DNAUPD again. */
           /* --------------------------------------- */

	goto L20;

    } else if (ido == 1) {

           /* --------------------------------------- */
           /* Perform y <-- OP*x = inv[A-sigma*M]*M*x */
           /* M*x has been saved in workd(ipntr(3)).  */
           /* The user only need the linear system    */
           /* solver here that takes workd(ipntr(3))  */
           /* as input, and returns the result to     */
           /* workd(ipntr(2)).                        */
           /* --------------------------------------- */

	dcopy_(&n, &workd[ipntr[2] - 1], &c__1, &workd[ipntr[1] - 1], &c__1);
	dgttrs_("N", &n, &c__1, dl, dd, du, du2, ipiv, &workd[ipntr[1] - 1], &
		n, &ierr, (ftnlen)1);
	if (ierr != 0) {
	    printf(" \n");
	    printf(" ERROR with _gttrs in _NDRV4.\n");
	    printf(" \n");
	    return ierr;
	}

           /* --------------------------------------- */
           /* L O O P   B A C K to call DNAUPD again. */
           /* --------------------------------------- */

	goto L20;

    } else if (ido == 2) {

           /* ------------------------------------------- */
           /*          Perform  y <--- M*x                */
           /* Need matrix vector multiplication routine   */
           /* here that takes workd(ipntr(1)) as input    */
           /* and returns the result to workd(ipntr(2)).  */
           /* ------------------------------------------- */

	dndrv4_mv_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

           /* --------------------------------------- */
           /* L O O P   B A C K to call DNAUPD again. */
           /* --------------------------------------- */

	goto L20;

    }

     /* --------------------------------------- */
     /* Either we have convergence, or there is */
     /* an error.                               */
     /* --------------------------------------- */

    if (info < 0) {

        /* ------------------------ */
        /* Error message, check the */
        /* documentation in DNAUPD. */
        /* ------------------------ */

	printf(" \n");
	printf(" Error with _naupd info = %d\n", info);
	printf(" Check the documentation in _naupd.\n");
	printf(" \n");

    } else {

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
	dneupd_(&rvec, "A", select, d, &d[25], v, &c__256, &sigmar, &
		sigmai, workev, bmat, &n, which, &nev, &tol, resid, &ncv, v, &
		c__256, iparam, ipntr, workd, workl, &lworkl, &ierr, (ftnlen)
		1, (ftnlen)1, (ftnlen)2);

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

	if (ierr != 0) {

            /* ---------------------------------- */
            /* Error condition:                   */
            /* Check the documentation of DNEUPD. */
            /* ---------------------------------- */

	    printf(" \n");
	    printf(" Error with _neupd info = %d\n", ierr);
	    printf(" Check the documentation of _neupd. \n");
	    printf(" \n");

	} else {

	    first = true;
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

		if (d[j + 24] == 0.) {

                    /* ------------------ */
                    /* Ritz value is real */
                    /* ------------------ */

			dndrv4_av_(&n, &v[(j << 8) - 256], ax);
			dndrv4_mv_(&n, &v[(j << 8) - 256], mx);
		    d__1 = -d[j - 1];
		    daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
		    d[j + 49] = dnrm2_(&n, ax, &c__1);
		    d[j + 49] /= (d__1 = d[j - 1], abs(d__1));

		} else if (first) {

                    /* ---------------------- */
                    /* Ritz value is complex. */
                    /* Residual of one Ritz   */
                    /* value of the conjugate */
                    /* pair is computed.      */
                    /* ---------------------- */

			dndrv4_av_(&n, &v[(j << 8) - 256], ax);
			dndrv4_mv_(&n, &v[(j << 8) - 256], mx);
		    d__1 = -d[j - 1];
		    daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
			dndrv4_mv_(&n, &v[(j + 1 << 8) - 256], mx);
		    daxpy_(&n, &d[j + 24], mx, &c__1, ax, &c__1);
		    d[j + 49] = dnrm2_(&n, ax, &c__1);
			dndrv4_av_(&n, &v[(j + 1 << 8) - 256], ax);
			dndrv4_mv_(&n, &v[(j + 1 << 8) - 256], mx);
		    d__1 = -d[j - 1];
		    daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
			dndrv4_mv_(&n, &v[(j << 8) - 256], mx);
		    d__1 = -d[j + 24];
		    daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
		    d__1 = dnrm2_(&n, ax, &c__1);
		    d[j + 49] = dlapy2_(&d[j + 49], &d__1);
		    d[j + 49] /= dlapy2_(&d[j - 1], &d[j + 24]);
		    d[j + 50] = d[j + 49];
		    first = false;
		} else {
		    first = true;
		}

/* L30: */
	    }

            /* --------------------------- */
            /* Display computed residuals. */
            /* --------------------------- */

	    dmout_(&nconv, &c__3, d, &c__25, &c_n6, "Ritz values (Real,Imag) and relative residuals");

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
	printf(" _NDRV4 \n");
	printf(" ====== \n");
	printf(" \n");
	printf(" Size of the matrix is %d\n", n);
	printf(" The number of Ritz values requested is %d\n", nev);
	printf(" The number of Arnoldi vectors generated (NCV) is %d\n", ncv);
	printf(" What portion of the spectrum: %s\n", which);
	printf(" The number of converged Ritz values is %d\n", nconv);
	printf(" The number of Implicit Arnoldi update iterations taken is %d\n", iparam[2]);
	printf(" The number of OP*x is %d\n", iparam[8]);
	printf(" The convergence criterion is %f\n", tol);
	printf(" \n");

    }

     /* ------------------------- */
     /* Done with program dndrv4. */
     /* ------------------------- */

    return 0;
} /* MAIN__ */

/* ========================================================================== */

/*     matrix vector multiplication subroutine */

int dndrv4_mv_(int32_t *n, double *v, double *w)
{
    /* System generated locals */
    int32_t i__1;

    /* Local variables */
    double h;
    int32_t j;

/*     Compute the matrix vector multiplication y<---M*x */
/*     where M is mass matrix formed by using piecewise linear elements */
/*     on [0,1]. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    w[1] = (v[1] * 4. + v[2] * 1.) / 6.;
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j) {
	w[j] = (v[j - 1] * 1. + v[j] * 4. + v[j + 1] * 1.) / 6.;
/* L10: */
    }
    w[*n] = (v[*n - 1] * 1. + v[*n] * 4.) / 6.;

    h = 1. / (double) (*n + 1);
    dscal_(n, &h, &w[1], &c__1);
    return 0;
} /* mv_ */

/* ------------------------------------------------------------------ */
int dndrv4_av_(int32_t *n, double *v, double *w)
{
    /* System generated locals */
    int32_t i__1;

    /* Local variables */
    double h;
    int32_t j;
    double s, dd, dl, du;

/*     Compute the matrix vector multiplication y<---A*x */
/*     where A is obtained from the finite element discretization of the */
/*     1-dimensional convection diffusion operator */
/*                     d^u/dx^2 + rho*(du/dx) */
/*     on the interval [0,1] with zero Dirichlet boundary condition */
/*     using linear elements. */
/*     This routine is only used in residual calculation. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    h = 1. / (double) (*n + 1);
    s = convct_1.rho / 2.;
    dd = 2. / h;
    dl = -1. / h - s;
    du = -1. / h + s;

    w[1] = dd * v[1] + du * v[2];
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j) {
	w[j] = dl * v[j - 1] + dd * v[j] + du * v[j + 1];
/* L10: */
    }
    w[*n] = dl * v[*n - 1] + dd * v[*n];
    return 0;
} /* av_ */

