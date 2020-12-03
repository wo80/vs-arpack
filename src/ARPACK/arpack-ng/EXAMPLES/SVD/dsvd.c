/* EXAMPLES\SVD\dsvd.f -- translated by f2c (version 20100827). */

#include "arpack.h"

int dsvd()
{
    /* System generated locals */
    int32_t i__1;
    double d__1;

    /* Builtin functions */

    int32_t s_wsle(cilist *), do_lio(int32_t *, int32_t *, char *, ftnlen), 
	    e_wsle(void);
    double sqrt(double);

    /* Local variables */
    int32_t j, m, n;
    double s[50]	/* was [25][2] */, u[5000]	/* was [500][10] */, 
	    v[6250]	/* was [250][25] */;
    double ax[500];
    int32_t ido, ncv, nev;
    double tol;
    char* bmat;
    int32_t info;
    bool rvec;
    int32_t ierr;
    double temp;
    int32_t mode1;
    double sigma;
    char* which;
    double resid[250];
    int32_t nconv;
    double workd[750];
    int32_t ipntr[11];
    double workl[825];
    int32_t iparam[11];
    bool select[25];
    int32_t ishfts, maxitr, lworkl;

    /* Fortran I/O blocks */

/*     This example program is intended to illustrate the */
/*     the use of ARPACK to compute the Singular Value Decomposition. */

/*     This code shows how to use ARPACK to find a few of the */
/*     largest singular values(sigma) and corresponding right singular */
/*     vectors (v) for the the matrix A by solving the symmetric problem: */

/*                        (A'*A)*v = sigma*v */

/*     where A is an m by n real matrix. */

/*     This code may be easily modified to estimate the 2-norm */
/*     condition number  largest(sigma)/smallest(sigma) by setting */
/*     which = 'BE' below.  This will ask for a few of the smallest */
/*     and a few of the largest singular values simultaneously. */
/*     The condition number could then be estimated by taking */
/*     the ratio of the largest and smallest singular values. */

/*     This formulation is appropriate when  m  .ge.  n. */
/*     Reverse the roles of A and A' in the case that  m .le. n. */

/*     The main points illustrated here are */

/*        1) How to declare sufficient memory to find NEV */
/*           largest singular values of A . */

/*        2) Illustration of the reverse communication interface */
/*           needed to utilize the top level ARPACK routine DSAUPD */
/*           that computes the quantities needed to construct */
/*           the desired singular values and vectors(if requested). */

/*        3) How to extract the desired singular values and vectors */
/*           using the ARPACK routine DSEUPD. */

/*        4) How to construct the left singular vectors U from the */
/*           right singular vectors V to obtain the decomposition */

/*                        A*V = U*S */

/*           where S = diag(sigma_1, sigma_2, ..., sigma_k). */

/*     The only thing that must be supplied in order to use this */
/*     routine on your problem is to change the array dimensions */
/*     appropriately, to specify WHICH singular values you want to */
/*     compute and to supply a the matrix-vector products */

/*                         w <-  Ax */
/*                         y <-  A'w */

/*     in place of the calls  to AV( ) and ATV( ) respectively below. */

/*     Further documentation is available in the header of DSAUPD */
/*     which may be found in the SRC directory. */

/*     This codes implements */

/* \Example-1 */
/*     ... Suppose we want to solve A'A*v = sigma*v in regular mode, */
/*         where A is derived from the simplest finite difference */
/*         discretization of the 2-dimensional kernel  K(s,t)dt  where */

/*                 K(s,t) =  s(t-1)   if 0 .le. s .le. t .le. 1, */
/*                           t(s-1)   if 0 .le. t .lt. s .le. 1. */

/*         See subroutines AV  and ATV for details. */
/*     ... OP = A'*A  and  B = I. */
/*     ... Assume "call av (n,x,y)" computes y = A*x */
/*     ... Assume "call atv (n,y,w)" computes w = A'*y */
/*     ... Assume exact shifts are used */
/*     ... */
/**
 * \BeginLib
 *
 * \Routines called:
 *     dsaupd  ARPACK reverse communication interface routine.
 *     dseupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     dnrm2   Level 1 BLAS that computes the norm of a vector.
 *     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
 *     dscal   Level 1 BLAS thst computes x <- x*alpha.
 *     dcopy   Level 1 BLAS thst computes y <- x.
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
 * FILE: svd.F   SID: 2.4   DATE OF SID: 10/17/00   RELEASE: 2
 *
 * \Remarks
 *     1. None
 *
 * \EndLib
 */
     /* ---------------------------------------------------- */
     /* Storage Declarations:                                */
     /*                                                      */
     /* It is assumed that A is M by N with M .ge. N.        */
     /*                                                      */
     /* The maximum dimensions for all arrays are            */
     /* set here to accommodate a problem size of            */
     /* M .le. MAXM  and  N .le. MAXN                        */
     /*                                                      */
     /* The NEV right singular vectors will be computed in   */
     /* the N by NCV array V.                                */
     /*                                                      */
     /* The NEV left singular vectors will be computed in    */
     /* the M by NEV array U.                                */
     /*                                                      */
     /* NEV is the number of singular values requested.      */
     /*     See specifications for ARPACK usage below.       */
     /*                                                      */
     /* NCV is the largest number of basis vectors that will */
     /*     be used in the Implicitly Restarted Arnoldi      */
     /*     Process.  Work per major iteration is            */
     /*     proportional to N*NCV*NCV.                       */
     /*                                                      */
     /* You must set:                                        */
     /*                                                      */
     /* MAXM:   Maximum number of rows of the A allowed.     */
     /* MAXN:   Maximum number of columns of the A allowed.  */
     /* MAXNEV: Maximum NEV allowed                          */
     /* MAXNCV: Maximum NCV allowed                          */
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

    m = 500;
    n = 100;

     /* ---------------------------------------------- */
     /* Specifications for ARPACK usage are set        */
     /* below:                                         */
     /*                                                */
     /*    1) NEV = 4 asks for 4 singular values to be */
     /*       computed.                                */
     /*                                                */
     /*    2) NCV = 20 sets the length of the Arnoldi  */
     /*       factorization                            */
     /*                                                */
     /*    3) This is a standard problem               */
     /*         (indicated by bmat  = 'I')             */
     /*                                                */
     /*    4) Ask for the NEV singular values of       */
     /*       largest magnitude                        */
     /*         (indicated by which = 'LM')            */
     /*       See documentation in DSAUPD for the      */
     /*       other options SM, BE.                    */
     /*                                                */
     /* Note: NEV and NCV must satisfy the following   */
     /*       conditions:                              */
     /*                 NEV <= MAXNEV,                 */
     /*             NEV + 1 <= NCV <= MAXNCV           */
     /* ---------------------------------------------- */

    nev = 4;
    ncv = 10;
    bmat = "I";
    which = "LM";

    if (n > 250) {
	printf(" ERROR with _SVD: N is greater than MAXN \n");
	return 0;
    } else if (m > 500) {
	printf(" ERROR with _SVD: M is greater than MAXM \n");
	return 0;
    } else if (nev > 10) {
	printf(" ERROR with _SVD: NEV is greater than MAXNEV \n");
	return 0;
    } else if (ncv > 25) {
	printf(" ERROR with _SVD: NCV is greater than MAXNCV \n");
	return 0;
    }

     /* --------------------------------------------------- */
     /* Specification of stopping rules and initial         */
     /* conditions before calling DSAUPD                    */
     /*                                                     */
     /*           abs(sigmaC - sigmaT) < TOL*abs(sigmaC)    */
     /*               computed   true                       */
     /*                                                     */
     /*      If TOL .le. 0,  then TOL <- macheps            */
     /*              (machine precision) is used.           */
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
     /* --------------------------------------------------- */

    lworkl = ncv * (ncv + 8);
    tol = 0.;
    info = 0;
    ido = 0;

     /* ------------------------------------------------- */
     /* Specification of Algorithm Mode:                  */
     /*                                                   */
     /* This program uses the exact shift strategy        */
     /* (indicated by setting IPARAM(1) = 1.)             */
     /* IPARAM(3) specifies the maximum number of Arnoldi */
     /* iterations allowed.  Mode 1 of DSAUPD is used     */
     /* (IPARAM(7) = 1). All these options can be changed */
     /* by the user. For details see the documentation in */
     /* DSAUPD.                                           */
     /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = n;
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

    dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__250, iparam, ipntr, workd, workl, &lworkl, &info);

    if (ido == -1 || ido == 1) {

           /* ------------------------------------- */
           /* Perform matrix vector multiplications */
           /*              w <--- A*x       (av())  */
           /*              y <--- A'*w      (atv()) */
           /* The user should supply his/her own    */
           /* matrix vector multiplication routines */
           /* here that takes workd(ipntr(1)) as    */
           /* the input, and returns the result in  */
           /* workd(ipntr(2)).                      */
           /* ------------------------------------- */

	dsvd_av_(&m, &n, &workd[ipntr[0] - 1], ax);
	dsvd_atv_(&m, &n, ax, &workd[ipntr[1] - 1]);

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

        /* ------------------------------------------ */
        /* No fatal errors occurred.                  */
        /* Post-Process using DSEUPD.                 */
        /*                                            */
        /* Computed singular values may be extracted. */
        /*                                            */
        /* Singular vectors may also be computed now  */
        /* if desired.  (indicated by rvec = .true.)  */
        /*                                            */
        /* The routine DSEUPD now called to do this   */
        /* post processing                            */
        /* ------------------------------------------ */

	rvec = true;

	dseupd_(&rvec, "All", select, s, v, &c__250, &sigma, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__250, iparam, ipntr, workd, workl, &lworkl, &ierr);

        /* --------------------------------------------- */
        /* Singular values are returned in the first     */
        /* column of the two dimensional array S         */
        /* and the corresponding right singular vectors  */
        /* are returned in the first NEV columns of the  */
        /* two dimensional array V as requested here.    */
        /* --------------------------------------------- */

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

		s[j - 1] = sqrt(s[j - 1]);

              /* --------------------------- */
              /* Compute the left singular   */
              /* vectors from the formula    */
              /*                             */
              /*     u = Av/sigma            */
              /*                             */
              /* u should have norm 1 so     */
              /* divide by norm(Av) instead. */
              /* --------------------------- */

		dsvd_av_(&m, &n, &v[j * 250 - 250], ax);
		dcopy_(&m, ax, &c__1, &u[j * 500 - 500], &c__1);
		temp = 1. / dnrm2_(&m, &u[j * 500 - 500], &c__1);
		dscal_(&m, &temp, &u[j * 500 - 500], &c__1);

              /* ------------------------- */
              /*                           */
              /* Compute the residual norm */
              /*                           */
              /*   ||  A*v - sigma*u ||    */
              /*                           */
              /* for the NCONV accurately  */
              /* computed singular values  */
              /* and vectors.  (iparam(5)  */
              /* indicates how many are    */
              /* accurate to the requested */
              /* tolerance).               */
              /* Store the result in 2nd   */
              /* column of array S.        */
              /* ------------------------- */

		d__1 = -s[j - 1];
		daxpy_(&m, &d__1, &u[j * 500 - 500], &c__1, ax, &c__1);
		s[j + 24] = dnrm2_(&m, ax, &c__1);

/* L20: */
	    }

           /* ----------------------------- */
           /* Display computed residuals    */
           /* ----------------------------- */

	    dmout_(&nconv, &c__2, s, &c__25, &c_n6, "Singular values and direct residuals");
	}

        /* ---------------------------------------- */
        /* Print additional convergence information */
        /* ---------------------------------------- */

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
	printf(" _SVD \n");
	printf(" ==== \n");
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

     /* ----------------------- */
     /* Done with program dsvd. */
     /* ----------------------- */

    return 0;
} /* MAIN__ */

/* ------------------------------------------------------------------ */
/*     matrix vector subroutines */

/*     The matrix A is derived from the simplest finite difference */
/*     discretization of the integral operator */

/*                     f(s) = integral(K(s,t)x(t)dt). */

/*     Thus, the matrix A is a discretization of the 2-dimensional kernel */
/*     K(s,t)dt, where */

/*                 K(s,t) =  s(t-1)   if 0 .le. s .le. t .le. 1, */
/*                           t(s-1)   if 0 .le. t .lt. s .le. 1. */

/*     Thus A is an m by n matrix with entries */

/*                 A(i,j) = k*(si)*(tj - 1)  if i .le. j, */
/*                          k*(tj)*(si - 1)  if i .gt. j */

/*     where si = i/(m+1)  and  tj = j/(n+1)  and k = 1/(n+1). */

/* ------------------------------------------------------------------- */

int dsvd_av_(int32_t *m, int32_t *n, double *x, double *w)
{
    /* System generated locals */
    int32_t i__1, i__2;

    /* Local variables */
    double h;
    int32_t i, j;
    double k, s, t;

/*     computes  w <- A*x */

    /* Parameter adjustments */
    --w;
    --x;

    /* Function Body */
    h = 1. / (double) (*m + 1);
    k = 1. / (double) (*n + 1);
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
	w[i] = 0.;
/* L5: */
    }
    t = 0.;

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	t += k;
	s = 0.;
	i__2 = j;
	for (i = 1; i <= i__2; ++i) {
	    s += h;
	    w[i] += k * s * (t - 1.) * x[j];
/* L10: */
	}
	i__2 = *m;
	for (i = j + 1; i <= i__2; ++i) {
	    s += h;
	    w[i] += k * t * (s - 1.) * x[j];
/* L20: */
	}
/* L30: */
    }

    return 0;
} /* av_ */

/* ------------------------------------------------------------------- */

int dsvd_atv_(int32_t *m, int32_t *n, double *w, double *y)
{
    /* System generated locals */
    int32_t i__1, i__2;

    /* Local variables */
    double h;
    int32_t i, j;
    double k, s, t;

/*     computes  y <- A'*w */

    /* Parameter adjustments */
    --w;
    --y;

    /* Function Body */
    h = 1. / (double) (*m + 1);
    k = 1. / (double) (*n + 1);
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	y[i] = 0.;
/* L5: */
    }
    t = 0.;

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	t += k;
	s = 0.;
	i__2 = j;
	for (i = 1; i <= i__2; ++i) {
	    s += h;
	    y[j] += k * s * (t - 1.) * w[i];
/* L10: */
	}
	i__2 = *m;
	for (i = j + 1; i <= i__2; ++i) {
	    s += h;
	    y[j] += k * t * (s - 1.) * w[i];
/* L20: */
	}
/* L30: */
    }

    return 0;
} /* atv_ */

