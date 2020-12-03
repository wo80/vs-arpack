/* EXAMPLES\NONSYM\dndrv1.f -- translated by f2c (version 20100827). */

#include <stdlib.h>
#include "arpack.h"

int dndrv1()
{
    /* System generated locals */
    int32_t i__1;
    double d__1;

    /* Local variables */
    double d[90]	/* was [30][3] */;
    int32_t j, n;
    double *v	/* was [256][30] */;
    double ax[256];
    int32_t nx, ido, ncv, nev;
    double tol;
    char* bmat;
    int32_t mode, info;
    bool rvec;
    int32_t ierr;
    char* which;
    double* resid;
    int32_t nconv;
    double* workd;
    bool first;
    int32_t ipntr[14];
    double *workl;
    int32_t iparam[11];
    double sigmai;
    bool select[30];
    double sigmar;
    int32_t ishfts, maxitr, lworkl;
    double workev[90];

    resid = (double*)malloc(256 * sizeof(double));
    v = (double*)malloc(7680 * sizeof(double));
    workl = (double*)malloc(2880 * sizeof(double));
    workd = (double*)malloc(768 * sizeof(double));

    /* Fortran I/O blocks */

/*     Example program to illustrate the idea of reverse communication */
/*     for a standard nonsymmetric eigenvalue problem. */

/*     We implement example one of ex-nonsym.doc in DOCUMENTS directory */

/* \Example-1 */
/*     ... Suppose we want to solve A*x = lambda*x in regular mode, */
/*         where A is obtained from the standard central difference */
/*         discretization of the convection-diffusion operator */
/*                 (Laplacian u) + rho*(du / dx) */
/*         on the unit square [0,1]x[0,1] with zero Dirichlet boundary */
/*         condition. */

/*     ... OP = A  and  B = I. */

/*     ... Assume "call av (nx,x,y)" computes y = A*x.c */

/*     ... Use mode 1 of DNAUPD. */
/**
 * \BeginLib
 *
 * \Routines called:
 *     dnaupd  ARPACK reverse communication interface routine.
 *     dneupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
 *     dnrm2   Level 1 BLAS that computes the norm of a vector.
 *     av      Matrix vector multiplication routine that computes A*x.
 *     tv      Matrix vector multiplication routine that computes T*x,
 *             where T is a tridiagonal matrix.  It is used in routine
 *             av.
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
 * FILE: ndrv1.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
 *
 * \Remarks
 *     1. None
 *
 * \EndLib
 */
     /* --------------------------- */
     /* Define maximum dimensions   */
     /* for all arrays.             */
     /* MAXN:   Maximum dimension   */
     /*         of the A allowed.   */
     /* MAXNEV: Maximum NEV allowed */
     /* MAXNCV: Maximum NCV allowed */
     /* --------------------------- */

     /* --------------------- */
     /* Executable Statements */
     /* --------------------- */

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
    which = "SM";

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
     /* iterations allowed.  Mode 1 of DNAUPD is used     */
     /* (IPARAM(7) = 1). All these options can be changed */
     /* by the user. For details see the documentation in */
     /* DNAUPD.                                           */
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
        /* Repeatedly call the routine DNAUPD and take */
        /* actions indicated by parameter IDO until    */
        /* either convergence is indicated or maxitr   */
        /* has been exceeded.                          */
        /* ------------------------------------------- */

    dnaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &info);

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

	dndrv1_av_(&nx, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

           /* --------------------------------------- */
           /* L O O P   B A C K to call DNAUPD again. */
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
        /* documentation in DNAUPD. */
        /* ------------------------ */

	printf(" \n");
	printf(" Error with _naupd info = %d\n", info);
	printf(" Check the documentation of _naupd\n");
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

	dneupd_(&rvec, "A", select, d, &d[30], v, &c__256, &sigmar, &sigmai, workev, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &ierr);

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

		if (d[j + 29] == 0.) {

                  /* ------------------ */
                  /* Ritz value is real */
                  /* ------------------ */

			dndrv1_av_(&nx, &v[(j << 8) - 256], ax);
		    d__1 = -d[j - 1];
		    daxpy_(&n, &d__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
		    d[j + 59] = dnrm2_(&n, ax, &c__1);
		    d[j + 59] /= (d__1 = d[j - 1], abs(d__1));

		} else if (first) {

                  /* ---------------------- */
                  /* Ritz value is complex. */
                  /* Residual of one Ritz   */
                  /* value of the conjugate */
                  /* pair is computed.      */
                  /* ---------------------- */

			dndrv1_av_(&nx, &v[(j << 8) - 256], ax);
		    d__1 = -d[j - 1];
		    daxpy_(&n, &d__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
		    daxpy_(&n, &d[j + 29], &v[(j + 1 << 8) - 256], &c__1, ax, &c__1);
		    d[j + 59] = dnrm2_(&n, ax, &c__1);
			dndrv1_av_(&nx, &v[(j + 1 << 8) - 256], ax);
		    d__1 = -d[j + 29];
		    daxpy_(&n, &d__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
		    d__1 = -d[j - 1];
		    daxpy_(&n, &d__1, &v[(j + 1 << 8) - 256], &c__1, ax, &c__1);
		    d__1 = dnrm2_(&n, ax, &c__1);
		    d[j + 59] = dlapy2_(&d[j + 59], &d__1);
		    d[j + 59] /= dlapy2_(&d[j - 1], &d[j + 29]);
		    d[j + 60] = d[j + 59];
		    first = false;
		} else {
		    first = true;
		}

/* L20: */
	    }

            /* --------------------------- */
            /* Display computed residuals. */
            /* --------------------------- */

	    dmout_(&nconv, &c__3, d, &c__30, &c_n6, "Ritz values (Real,Imag) and relative residuals");
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
	printf(" The convergence criterion is %f\n", tol);
	printf(" \n");

    }

    free(resid);
    free(v);
    free(workl);
    free(workd);

     /* ------------------------- */
     /* Done with program dndrv1. */
     /* ------------------------- */

    return 0;
} /* MAIN__ */

/* ========================================================================== */

/*     matrix vector subroutine */

/*     The matrix used is the 2 dimensional convection-diffusion */
/*     operator discretized using central difference. */

int dndrv1_av_(int32_t *nx, double *v, double *w)
{
    /* System generated locals */
    int32_t i__1;
    double d__1;

    /* Local variables */
    int32_t j;
    double h2;
    int32_t lo;

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
    h2 = 1. / (double) ((*nx + 1) * (*nx + 1));

	dndrv1_tv_(nx, &v[1], &w[1]);
    d__1 = -1. / h2;
    daxpy_(nx, &d__1, &v[*nx + 1], &c__1, &w[1], &c__1);

    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j) {
	lo = (j - 1) * *nx;
	dndrv1_tv_(nx, &v[lo + 1], &w[lo + 1]);
	d__1 = -1. / h2;
	daxpy_(nx, &d__1, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);
	d__1 = -1. / h2;
	daxpy_(nx, &d__1, &v[lo + *nx + 1], &c__1, &w[lo + 1], &c__1);
/* L10: */
    }

    lo = (*nx - 1) * *nx;
	dndrv1_tv_(nx, &v[lo + 1], &w[lo + 1]);
    d__1 = -1. / h2;
    daxpy_(nx, &d__1, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);

    return 0;
} /* av_ */

/* ========================================================================= */
int dndrv1_tv_(int32_t *nx, double *x, double *y)
{
    /* System generated locals */
    int32_t i__1;

    /* Local variables */
    double h;
    int32_t j;
    double h2, dd, dl, du;

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
    h = 1. / (double) (*nx + 1);
    h2 = h * h;
    dd = 4. / h2;
    dl = -1. / h2 - 0. / h;
    du = -1. / h2 + 0. / h;

    y[1] = dd * x[1] + du * x[2];
    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j) {
	y[j] = dl * x[j - 1] + dd * x[j] + du * x[j + 1];
/* L10: */
    }
    y[*nx] = dl * x[*nx - 1] + dd * x[*nx];
    return 0;
} /* tv_ */

