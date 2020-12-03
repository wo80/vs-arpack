/* EXAMPLES\SIMPLE\dssimp.f -- translated by f2c (version 20100827). */

#include "arpack.h"

int dssimp()
{
    /* System generated locals */
    int32_t i__1;
    double d__1;

    /* Builtin functions */

    int32_t s_wsle(cilist *), do_lio(int32_t *, int32_t *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    double d[50]	/* was [25][2] */;
    int32_t j, n;
    double v[6400]	/* was [256][25] */;
    double ax[256];
    int32_t nx, ido, ncv, nev;
    double tol;
    char* bmat;
    int32_t info;
    bool rvec;
    int32_t ierr, mode1;
    double sigma;
    char* which;
    double resid[256];
    int32_t nconv;
    double workd[768];
    int32_t ipntr[11];
    double workl[825];
    int32_t iparam[11];
    bool select[25];
    int32_t ishfts, maxitr, lworkl;

    /* Fortran I/O blocks */

/*     This example program is intended to illustrate the */
/*     simplest case of using ARPACK in considerable detail. */
/*     This code may be used to understand basic usage of ARPACK */
/*     and as a template for creating an interface to ARPACK. */

/*     This code shows how to use ARPACK to find a few eigenvalues */
/*     (lambda) and corresponding eigenvectors (x) for the standard */
/*     eigenvalue problem: */

/*                        A*x = lambda*x */

/*     where A is an n by n real symmetric matrix. */

/*     The main points illustrated here are */

/*        1) How to declare sufficient memory to find NEV */
/*           eigenvalues of largest magnitude.  Other options */
/*           are available. */

/*        2) Illustration of the reverse communication interface */
/*           needed to utilize the top level ARPACK routine DSAUPD */
/*           that computes the quantities needed to construct */
/*           the desired eigenvalues and eigenvectors(if requested). */

/*        3) How to extract the desired eigenvalues and eigenvectors */
/*           using the ARPACK routine DSEUPD. */

/*     The only thing that must be supplied in order to use this */
/*     routine on your problem is to change the array dimensions */
/*     appropriately, to specify WHICH eigenvalues you want to compute */
/*     and to supply a matrix-vector product */

/*                         w <-  Av */

/*     in place of the call to AV( ) below. */

/*     Once usage of this routine is understood, you may wish to explore */
/*     the other available options to improve convergence, to solve generalized */
/*     problems, etc.  Look at the file ex-sym.doc in DOCUMENTS directory. */
/*     This codes implements */

/* \Example-1 */
/*     ... Suppose we want to solve A*x = lambda*x in regular mode, */
/*         where A is derived from the central difference discretization */
/*         of the 2-dimensional Laplacian on the unit square with */
/*         zero Dirichlet boundary condition. */
/*     ... OP = A  and  B = I. */
/*     ... Assume "call av (n,x,y)" computes y = A*x */
/*     ... Use mode 1 of DSAUPD. */
/**
 * \BeginLib
 *
 * \Routines called:
 *     dsaupd  ARPACK reverse communication interface routine.
 *     dseupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     dnrm2   Level 1 BLAS that computes the norm of a vector.
 *     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
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
 * FILE: ssimp.F   SID: 2.6   DATE OF SID: 10/17/00   RELEASE: 2
 *
 * \Remarks
 *     1. None
 *
 * \EndLib
 */
     /* ---------------------------------------------------- */
     /* Storage Declarations:                                */
     /*                                                      */
     /* The maximum dimensions for all arrays are            */
     /* set here to accommodate a problem size of            */
     /* N .le. MAXN                                          */
     /*                                                      */
     /* NEV is the number of eigenvalues requested.          */
     /*     See specifications for ARPACK usage below.       */
     /*                                                      */
     /* NCV is the largest number of basis vectors that will */
     /*     be used in the Implicitly Restarted Arnoldi      */
     /*     Process.  Work per major iteration is            */
     /*     proportional to N*NCV*NCV.                       */
     /*                                                      */
     /* You must set:                                        */
     /*                                                      */
     /* MAXN:   Maximum dimension of the A allowed.          */
     /* MAXNEV: Maximum NEV allowed.                         */
     /* MAXNCV: Maximum NCV allowed.                         */
     /* ---------------------------------------------------- */

     /* --------------------- */
     /* Executable Statements */
     /* --------------------- */

     /* ----------------------------------------------- */
     /* The following include statement and assignments */
     /* initiate trace output from the internal         */
     /* actions of ARPACK.  See debug.doc in the        */
     /* DOCUMENTS directory for usage.  Initially, the  */
     /* most useful information will be a breakdown of  */
     /* time spent in the various stages of computation */
     /* given by setting msaupd = 1.                    */
     /* ----------------------------------------------- */

/* \SCCS Information: @(#) */
/* FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 */

     /* ------------------------------- */
     /* See debug.doc for documentation */
     /* ------------------------------- */
    debug_1.ndigit = -3;
    debug_1.logfil = 6;
    debug_1.msgets = 0;
    debug_1.msaitr = 0;
    debug_1.msapps = 0;
    debug_1.msaupd = 1;
    debug_1.msaup2 = 0;
    debug_1.mseigt = 0;
    debug_1.mseupd = 0;

     /* ----------------------------------------------- */
     /* The following sets dimensions for this problem. */
     /* ----------------------------------------------- */

    nx = 10;
    n = nx * nx;

     /* --------------------------------------------- */
     /*                                               */
     /* Specifications for ARPACK usage are set       */
     /* below:                                        */
     /*                                               */
     /*    1) NEV = 4  asks for 4 eigenvalues to be   */
     /*       computed.                               */
     /*                                               */
     /*    2) NCV = 20 sets the length of the Arnoldi */
     /*       factorization                           */
     /*                                               */
     /*    3) This is a standard problem              */
     /*         (indicated by bmat  = 'I')            */
     /*                                               */
     /*    4) Ask for the NEV eigenvalues of          */
     /*       largest magnitude                       */
     /*         (indicated by which = 'LM')           */
     /*       See documentation in DSAUPD for the     */
     /*       other options SM, LA, SA, LI, SI.       */
     /*                                               */
     /* Note: NEV and NCV must satisfy the following  */
     /* conditions:                                   */
     /*              NEV <= MAXNEV                    */
     /*          NEV + 1 <= NCV <= MAXNCV             */
     /* --------------------------------------------- */

    nev = 4;
    ncv = 20;
    bmat = "I";
    which = "LM";

    if (n > 256) {
	printf(" ERROR with _SSIMP: N is greater than MAXN \n");
	return 0;
    } else if (nev > 10) {
	printf(" ERROR with _SSIMP: NEV is greater than MAXNEV \n");
	return 0;
    } else if (ncv > 25) {
	printf(" ERROR with _SSIMP: NCV is greater than MAXNCV \n");
	return 0;
    }

     /* --------------------------------------------------- */
     /*                                                     */
     /* Specification of stopping rules and initial         */
     /* conditions before calling DSAUPD                    */
     /*                                                     */
     /* TOL  determines the stopping criterion.             */
     /*                                                     */
     /*      Expect                                         */
     /*           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) */
     /*               computed   true                       */
     /*                                                     */
     /*      If TOL .le. 0,  then TOL <- macheps            */
     /*           (machine precision) is used.              */
     /*                                                     */
     /* IDO  is the REVERSE COMMUNICATION parameter         */
     /*      used to specify actions to be taken on return  */
     /*      from DSAUPD. (See usage below.)                */
     /*                                                     */
     /*      It MUST initially be set to 0 before the first */
     /*      call to DSAUPD.                                */
     /*                                                     */
     /* INFO on entry specifies starting vector information */
     /*      and on return indicates error codes            */
     /*                                                     */
     /*      Initially, setting INFO=0 indicates that a     */
     /*      random starting vector is requested to         */
     /*      start the ARNOLDI iteration.  Setting INFO to  */
     /*      a nonzero value on the initial call is used    */
     /*      if you want to specify your own starting       */
     /*      vector (This vector must be placed in RESID.)  */
     /*                                                     */
     /* The work array WORKL is used in DSAUPD as           */
     /* workspace.  Its dimension LWORKL is set as          */
     /* illustrated below.                                  */
     /*                                                     */
     /* --------------------------------------------------- */

    lworkl = ncv * (ncv + 8);
    tol = 0.;
    info = 0;
    ido = 0;

     /* ------------------------------------------------- */
     /* Specification of Algorithm Mode:                  */
     /*                                                   */
     /* This program uses the exact shift strategy        */
     /* (indicated by setting PARAM(1) = 1).              */
     /* IPARAM(3) specifies the maximum number of Arnoldi */
     /* iterations allowed.  Mode 1 of DSAUPD is used     */
     /* (IPARAM(7) = 1). All these options can be changed */
     /* by the user. For details see the documentation in */
     /* DSAUPD.                                           */
     /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode1 = 1;

    iparam[0] = ishfts;

    iparam[2] = maxitr;

    iparam[6] = mode1;

     /* ---------------------------------------------- */
     /* M A I N   L O O P (Reverse communication loop) */
     /* ---------------------------------------------- */

L10:

        /* ------------------------------------------- */
        /* Repeatedly call the routine DSAUPD and take */
        /* actions indicated by parameter IDO until    */
        /* either convergence is indicated or maxitr   */
        /* has been exceeded.                          */
        /* ------------------------------------------- */

    dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, 
	    iparam, ipntr, workd, workl, &lworkl, &info, (ftnlen)1, (ftnlen)2)
	    ;

    if (ido == -1 || ido == 1) {

           /* ------------------------------------ */
           /* Perform matrix vector multiplication */
           /*              y <--- OP*x             */
           /* The user should supply his/her own   */
           /* matrix vector multiplication routine */
           /* here that takes workd(ipntr(1)) as   */
           /* the input, and return the result to  */
           /* workd(ipntr(2)).                     */
           /* ------------------------------------ */

	dssimp_av_(&nx, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

           /* --------------------------------------- */
           /* L O O P   B A C K to call DSAUPD again. */
           /* --------------------------------------- */

	goto L10;

    }

     /* -------------------------------------- */
     /* Either we have convergence or there is */
     /* an error.                              */
     /* -------------------------------------- */

    if (info < 0) {

        /* ------------------------ */
        /* Error message. Check the */
        /* documentation in DSAUPD. */
        /* ------------------------ */

	printf(" \n");
	printf(" Error with _saupd info = %d\n", info);
	printf(" Check documentation in _saupd \n");
	printf(" \n");

    } else {

        /* ----------------------------------------- */
        /* No fatal errors occurred.                 */
        /* Post-Process using DSEUPD.                */
        /*                                           */
        /* Computed eigenvalues may be extracted.    */
        /*                                           */
        /* Eigenvectors may be also computed now if  */
        /* desired.  (indicated by rvec = .true.)    */
        /*                                           */
        /* The routine DSEUPD now called to do this  */
        /* post processing (Other modes may require  */
        /* more complicated post processing than     */
        /* mode1.)                                   */
        /*                                           */
        /* ----------------------------------------- */

	rvec = true;

	dseupd_(&rvec, "All", select, d, v, &c__256, &sigma, bmat, &n, 
		which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, 
		workd, workl, &lworkl, &ierr, (ftnlen)3, (ftnlen)1, (ftnlen)2)
		;

         /* -------------------------------------------- */
         /* Eigenvalues are returned in the first column */
         /* of the two dimensional array D and the       */
         /* corresponding eigenvectors are returned in   */
         /* the first NCONV (=IPARAM(5)) columns of the  */
         /* two dimensional array V if requested.        */
         /* Otherwise, an orthogonal basis for the       */
         /* invariant subspace corresponding to the      */
         /* eigenvalues in D is returned in V.           */
         /* -------------------------------------------- */

	if (ierr != 0) {

            /* ---------------------------------- */
            /* Error condition:                   */
            /* Check the documentation of DSEUPD. */
            /* ---------------------------------- */

	    printf(" \n");
	    printf(" Error with _seupd info = %d\n", ierr);
	    printf(" Check the documentation of _seupd. \n");
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

		dssimp_av_(&nx, &v[(j << 8) - 256], ax);
		d__1 = -d[j - 1];
		daxpy_(&n, &d__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
		d[j + 24] = dnrm2_(&n, ax, &c__1);
		d[j + 24] /= (d__1 = d[j - 1], abs(d__1));

/* L20: */
	    }

            /* --------------------------- */
            /* Display computed residuals. */
            /* --------------------------- */

	    dmout_(&nconv, &c__2, d, &c__25, &c_n6, "Ritz values and relative residuals");
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
	printf(" _SSIMP \n");
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
     /* Done with program dssimp. */
     /* ------------------------- */

    return 0;
} /* MAIN__ */

/* ------------------------------------------------------------------ */
/*     matrix vector subroutine */

/*     The matrix used is the 2 dimensional discrete Laplacian on unit */
/*     square with zero Dirichlet boundary condition. */

/*     Computes w <--- OP*v, where OP is the nx*nx by nx*nx block */
/*     tridiagonal matrix */

/*                  | T -I          | */
/*                  |-I  T -I       | */
/*             OP = |   -I  T       | */
/*                  |        ...  -I| */
/*                  |           -I T| */

/*     The subroutine TV is called to computed y<---T*x. */

int dssimp_av_(int32_t *nx, double *v, double *w)
{
    /* System generated locals */
    int32_t i__1;
    double d__1;

    /* Local variables */
    int32_t j;
    double h2;
    int32_t n2, lo;

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
	dssimp_tv_(nx, &v[1], &w[1]);
    daxpy_(nx, &c_b138, &v[*nx + 1], &c__1, &w[1], &c__1);

    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j) {
	lo = (j - 1) * *nx;
	dssimp_tv_(nx, &v[lo + 1], &w[lo + 1]);
	daxpy_(nx, &c_b138, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);
	daxpy_(nx, &c_b138, &v[lo + *nx + 1], &c__1, &w[lo + 1], &c__1);
/* L10: */
    }

    lo = (*nx - 1) * *nx;
	dssimp_tv_(nx, &v[lo + 1], &w[lo + 1]);
    daxpy_(nx, &c_b138, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);

/*     Scale the vector w by (1/h^2), where h is the mesh size */

    n2 = *nx * *nx;
    h2 = 1. / (double) ((*nx + 1) * (*nx + 1));
    d__1 = 1. / h2;
    dscal_(&n2, &d__1, &w[1], &c__1);
    return 0;
} /* av_ */

/* ------------------------------------------------------------------- */
int dssimp_tv_(int32_t *nx, double *x, double *y)
{
    /* System generated locals */
    int32_t i__1;

    /* Local variables */
    int32_t j;
    double dd, dl, du;

/*     Compute the matrix vector multiplication y<---T*x */
/*     where T is a nx by nx tridiagonal matrix with DD on the */
/*     diagonal, DL on the subdiagonal, and DU on the superdiagonal. */

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    dd = 4.;
    dl = -1.;
    du = -1.;

    y[1] = dd * x[1] + du * x[2];
    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j) {
	y[j] = dl * x[j - 1] + dd * x[j] + du * x[j + 1];
/* L10: */
    }
    y[*nx] = dl * x[*nx - 1] + dd * x[*nx];
    return 0;
} /* tv_ */

