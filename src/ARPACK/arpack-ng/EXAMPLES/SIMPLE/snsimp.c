/* EXAMPLES\SIMPLE\snsimp.f -- translated by f2c (version 20100827). */

#include "arpack.h"

int snsimp()
{
    /* System generated locals */
    int32_t i__1;
    float r__1;

    /* Builtin functions */

    int32_t s_wsle(cilist *), do_lio(int32_t *, int32_t *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    float d[90]	/* was [30][3] */;
    int32_t j, n;
    float v[7680]	/* was [256][30] */;
    float ax[256];
    int32_t nx, ido, ncv, nev;
    float tol;
    char bmat[1];
    int32_t info;
    bool rvec;
    int32_t ierr, mode1;
    char which[2];
    float resid[256];
    int32_t nconv;
    float workd[768];
    bool first;
    int32_t ipntr[14];
    float workl[2880];
    int32_t iparam[11];
    float sigmai;
    bool select[30];
    float sigmar;
    int32_t ishfts, maxitr;
    int32_t lworkl;
    float workev[90];

    /* Fortran I/O blocks */
    static cilist io___7 = { 0, 6, 0, 0, 0 };
    static cilist io___8 = { 0, 6, 0, 0, 0 };
    static cilist io___9 = { 0, 6, 0, 0, 0 };
    static cilist io___23 = { 0, 6, 0, 0, 0 };
    static cilist io___24 = { 0, 6, 0, 0, 0 };
    static cilist io___25 = { 0, 6, 0, 0, 0 };
    static cilist io___26 = { 0, 6, 0, 0, 0 };
    static cilist io___34 = { 0, 6, 0, 0, 0 };
    static cilist io___35 = { 0, 6, 0, 0, 0 };
    static cilist io___36 = { 0, 6, 0, 0, 0 };
    static cilist io___37 = { 0, 6, 0, 0, 0 };
    static cilist io___42 = { 0, 6, 0, 0, 0 };
    static cilist io___43 = { 0, 6, 0, 0, 0 };
    static cilist io___44 = { 0, 6, 0, 0, 0 };
    static cilist io___45 = { 0, 6, 0, 0, 0 };
    static cilist io___46 = { 0, 6, 0, 0, 0 };
    static cilist io___47 = { 0, 6, 0, 0, 0 };
    static cilist io___48 = { 0, 6, 0, 0, 0 };
    static cilist io___49 = { 0, 6, 0, 0, 0 };
    static cilist io___50 = { 0, 6, 0, 0, 0 };
    static cilist io___51 = { 0, 6, 0, 0, 0 };
    static cilist io___52 = { 0, 6, 0, 0, 0 };
    static cilist io___53 = { 0, 6, 0, 0, 0 };
    static cilist io___54 = { 0, 6, 0, 0, 0 };
    static cilist io___55 = { 0, 6, 0, 0, 0 };
    static cilist io___56 = { 0, 6, 0, 0, 0 };
    static cilist io___57 = { 0, 6, 0, 0, 0 };
    static cilist io___58 = { 0, 6, 0, 0, 0 };
    static cilist io___59 = { 0, 6, 0, 0, 0 };
    static cilist io___60 = { 0, 6, 0, 0, 0 };

/*     This example program is intended to illustrate the */
/*     simplest case of using ARPACK in considerable detail. */
/*     This code may be used to understand basic usage of ARPACK */
/*     and as a template for creating an interface to ARPACK. */

/*     This code shows how to use ARPACK to find a few eigenvalues */
/*     (lambda) and corresponding eigenvectors (x) for the standard */
/*     eigenvalue problem: */

/*                        A*x = lambda*x */

/*     where A is a n by n real nonsymmetric matrix. */

/*     The main points illustrated here are */

/*        1) How to declare sufficient memory to find NEV */
/*           eigenvalues of largest magnitude.  Other options */
/*           are available. */

/*        2) Illustration of the reverse communication interface */
/*           needed to utilize the top level ARPACK routine SNAUPD */
/*           that computes the quantities needed to construct */
/*           the desired eigenvalues and eigenvectors(if requested). */

/*        3) How to extract the desired eigenvalues and eigenvectors */
/*           using the ARPACK routine SNEUPD. */

/*     The only thing that must be supplied in order to use this */
/*     routine on your problem is to change the array dimensions */
/*     appropriately, to specify WHICH eigenvalues you want to compute */
/*     and to supply a matrix-vector product */

/*                         w <-  Av */

/*     in place of the call to AV( )  below. */

/*     Once usage of this routine is understood, you may wish to explore */
/*     the other available options to improve convergence, to solve generalized */
/*     problems, etc.  Look at the file ex-nonsym.doc in DOCUMENTS directory. */
/*     This codes implements */

/* \Example-1 */
/*     ... Suppose we want to solve A*x = lambda*x in regular mode, */
/*         where A is obtained from the standard central difference */
/*         discretization of the convection-diffusion operator */
/*                 (Laplacian u) + rho*(du / dx) */
/*         on the unit square, with zero Dirichlet boundary condition. */

/*     ... OP = A  and  B = I. */
/*     ... Assume "call av (nx,x,y)" computes y = A*x */
/*     ... Use mode 1 of SNAUPD. */
/**
 * \BeginLib
 *
 * \Routines called:
 *     snaupd  ARPACK reverse communication interface routine.
 *     sneupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
 *     snrm2   Level 1 BLAS that computes the norm of a vector.
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
 * FILE: nsimp.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
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
     /* given by setting mnaupd = 1.                    */
     /* ----------------------------------------------- */

/* \SCCS Information: @(#) */
/* FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 */

     /* ------------------------------- */
     /* See debug.doc for documentation */
     /* ------------------------------- */
    debug_1.ndigit = -3;
    debug_1.logfil = 6;
    debug_1.mnaitr = 0;
    debug_1.mnapps = 0;
    debug_1.mnaupd = 1;
    debug_1.mnaup2 = 0;
    debug_1.mneigh = 0;
    debug_1.mneupd = 0;

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
     /*       factorization.                          */
     /*                                               */
     /*    3) This is a standard problem.             */
     /*         (indicated by bmat  = 'I')            */
     /*                                               */
     /*    4) Ask for the NEV eigenvalues of          */
     /*       largest magnitude.                      */
     /*         (indicated by which = 'LM')           */
     /*       See documentation in SNAUPD for the     */
     /*       other options SM, LR, SR, LI, SI.       */
     /*                                               */
     /* Note: NEV and NCV must satisfy the following  */
     /* conditions:                                   */
     /*              NEV <= MAXNEV                    */
     /*          NEV + 2 <= NCV <= MAXNCV             */
     /*                                               */
     /* --------------------------------------------- */

    nev = 4;
    ncv = 20;
    *bmat = 'I';
    strcpy(which, "LM");

    if (n > 256) {
	s_wsle(&io___7);
	do_lio(&c__9, &c__1, " ERROR with _NSIMP: N is greater than MAXN ", (ftnlen)43);
	e_wsle();
	goto L9000;
    } else if (nev > 12) {
	s_wsle(&io___8);
	do_lio(&c__9, &c__1, " ERROR with _NSIMP: NEV is greater than MAXNEV ", (ftnlen)47);
	e_wsle();
	goto L9000;
    } else if (ncv > 30) {
	s_wsle(&io___9);
	do_lio(&c__9, &c__1, " ERROR with _NSIMP: NCV is greater than MAXNCV ", (ftnlen)47);
	e_wsle();
	goto L9000;
    }

     /* --------------------------------------------------- */
     /*                                                     */
     /* Specification of stopping rules and initial         */
     /* conditions before calling SNAUPD                    */
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
     /*      from SNAUPD. (see usage below)                 */
     /*                                                     */
     /*      It MUST initially be set to 0 before the first */
     /*      call to SNAUPD.                                */
     /*                                                     */
     /* INFO on entry specifies starting vector information */
     /*      and on return indicates error codes            */
     /*                                                     */
     /*      Initially, setting INFO=0 indicates that a     */
     /*      random starting vector is requested to         */
     /*      start the ARNOLDI iteration.  Setting INFO to  */
     /*      a nonzero value on the initial call is used    */
     /*      if you want to specify your own starting       */
     /*      vector (This vector must be placed in RESID).  */
     /*                                                     */
     /* The work array WORKL is used in SNAUPD as           */
     /* workspace.  Its dimension LWORKL is set as          */
     /* illustrated below.                                  */
     /*                                                     */
     /* --------------------------------------------------- */

/* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 * 3 + ncv * 6;
    tol = 0.f;
    ido = 0;
    info = 0;

     /* ------------------------------------------------- */
     /* Specification of Algorithm Mode:                  */
     /*                                                   */
     /* This program uses the exact shift strategy        */
     /* (indicated by setting IPARAM(1) = 1).             */
     /* IPARAM(3) specifies the maximum number of Arnoldi */
     /* iterations allowed.  Mode 1 of SNAUPD is used     */
     /* (IPARAM(7) = 1). All these options can be changed */
     /* by the user. For details see the documentation in */
     /* SNAUPD.                                           */
     /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode1 = 1;

    iparam[0] = ishfts;

    iparam[2] = maxitr;

    iparam[6] = mode1;

     /* ----------------------------------------- */
     /* M A I N   L O O P (Reverse communication) */
     /* ----------------------------------------- */

L10:

        /* ------------------------------------------- */
        /* Repeatedly call the routine SNAUPD and take */
        /* actions indicated by parameter IDO until    */
        /* either convergence is indicated or maxitr   */
        /* has been exceeded.                          */
        /* ------------------------------------------- */

    snaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, 
	    iparam, ipntr, workd, workl, &lworkl, &info, (ftnlen)1, (ftnlen)2)
	    ;

    if (ido == -1 || ido == 1) {

           /* ----------------------------------------- */
           /* Perform matrix vector multiplication      */
           /*                y <--- Op*x                */
           /* The user should supply his/her own        */
           /* matrix vector multiplication routine here */
           /* that takes workd(ipntr(1)) as the input   */
           /* vector, and return the matrix vector      */
           /* product to workd(ipntr(2)).               */
           /* ----------------------------------------- */

	snsimp_av_(&nx, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

           /* --------------------------------------- */
           /* L O O P   B A C K to call SNAUPD again. */
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
        /* documentation in SNAUPD. */
        /* ------------------------ */

	s_wsle(&io___23);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___24);
	do_lio(&c__9, &c__1, " Error with _naupd, info = ", (ftnlen)27);
	do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(int32_t));
	e_wsle();
	s_wsle(&io___25);
	do_lio(&c__9, &c__1, " Check the documentation of _naupd", (ftnlen)34);
	e_wsle();
	s_wsle(&io___26);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();

    } else {

        /* ----------------------------------------- */
        /* No fatal errors occurred.                 */
        /* Post-Process using SNEUPD.                */
        /*                                           */
        /* Computed eigenvalues may be extracted.    */
        /*                                           */
        /* Eigenvectors may be also computed now if  */
        /* desired.  (indicated by rvec = .true.)    */
        /*                                           */
        /* The routine SNEUPD now called to do this  */
        /* post processing (Other modes may require  */
        /* more complicated post processing than     */
        /* mode1,)                                   */
        /*                                           */
        /* ----------------------------------------- */

	rvec = true;

	sneupd_(&rvec, "A", select, d, &d[30], v, &c__256, &sigmar, &
		sigmai, workev, bmat, &n, which, &nev, &tol, resid, &ncv, v, &
		c__256, iparam, ipntr, workd, workl, &lworkl, &ierr, (ftnlen)
		1, (ftnlen)1, (ftnlen)2);

        /* ---------------------------------------------- */
        /* The real parts of the eigenvalues are returned */
        /* in the first column of the two dimensional     */
        /* array D, and the IMAGINARY part are returned   */
        /* in the second column of D.  The corresponding  */
        /* eigenvectors are returned in the first         */
        /* NCONV (= IPARAM(5)) columns of the two         */
        /* dimensional array V if requested.  Otherwise,  */
        /* an orthogonal basis for the invariant subspace */
        /* corresponding to the eigenvalues in D is       */
        /* returned in V.                                 */
        /* ---------------------------------------------- */

	if (ierr != 0) {

           /* ---------------------------------- */
           /* Error condition:                   */
           /* Check the documentation of SNEUPD. */
           /* ---------------------------------- */

	    s_wsle(&io___34);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___35);
	    do_lio(&c__9, &c__1, " Error with _neupd, info = ", (ftnlen)27);
	    do_lio(&c__3, &c__1, (char *)&ierr, (ftnlen)sizeof(int32_t));
	    e_wsle();
	    s_wsle(&io___36);
	    do_lio(&c__9, &c__1, " Check the documentation of _neupd. ", (ftnlen)36);
	    e_wsle();
	    s_wsle(&io___37);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();

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
              /* eigenvectors.  (IPARAM(5) */
              /* indicates how many are    */
              /* accurate to the requested */
              /* tolerance)                */
              /* ------------------------- */

		if (d[j + 29] == 0.f) {

                 /* ------------------ */
                 /* Ritz value is real */
                 /* ------------------ */

			snsimp_av_(&nx, &v[(j << 8) - 256], ax);
		    r__1 = -d[j - 1];
		    saxpy_(&n, &r__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
		    d[j + 59] = snrm2_(&n, ax, &c__1);
		    d[j + 59] /= (r__1 = d[j - 1], dabs(r__1));

		} else if (first) {

                 /* ---------------------- */
                 /* Ritz value is complex. */
                 /* Residual of one Ritz   */
                 /* value of the conjugate */
                 /* pair is computed.      */
                 /* ---------------------- */

			snsimp_av_(&nx, &v[(j << 8) - 256], ax);
		    r__1 = -d[j - 1];
		    saxpy_(&n, &r__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
		    saxpy_(&n, &d[j + 29], &v[(j + 1 << 8) - 256], &c__1, 
			    ax, &c__1);
		    d[j + 59] = snrm2_(&n, ax, &c__1);
			snsimp_av_(&nx, &v[(j + 1 << 8) - 256], ax);
		    r__1 = -d[j + 29];
		    saxpy_(&n, &r__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
		    r__1 = -d[j - 1];
		    saxpy_(&n, &r__1, &v[(j + 1 << 8) - 256], &c__1, ax, &
			    c__1);
		    r__1 = snrm2_(&n, ax, &c__1);
		    d[j + 59] = slapy2_(&d[j + 59], &r__1);
		    d[j + 59] /= slapy2_(&d[j - 1], &d[j + 29]);
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

	    smout_(&c__6, &nconv, &c__3, d, &c__30, &c_n6, "Ritz values (Real, Imag) and residual residuals");
	}

        /* ----------------------------------------- */
        /* Print additional convergence information. */
        /* ----------------------------------------- */

	if (info == 1) {
	    s_wsle(&io___42);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___43);
	    do_lio(&c__9, &c__1, " Maximum number of iterations reached.", (ftnlen)38);
	    e_wsle();
	    s_wsle(&io___44);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	} else if (info == 3) {
	    s_wsle(&io___45);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___46);
	    do_lio(&c__9, &c__1, " No shifts could be applied during implicit", (ftnlen)43);
	    do_lio(&c__9, &c__1, " Arnoldi update, try increasing NCV.", (ftnlen)36);
	    e_wsle();
	    s_wsle(&io___47);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	}

	s_wsle(&io___48);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___49);
	do_lio(&c__9, &c__1, " _NSIMP ", (ftnlen)8);
	e_wsle();
	s_wsle(&io___50);
	do_lio(&c__9, &c__1, " ====== ", (ftnlen)8);
	e_wsle();
	s_wsle(&io___51);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___52);
	do_lio(&c__9, &c__1, " Size of the matrix is ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(int32_t));
	e_wsle();
	s_wsle(&io___53);
	do_lio(&c__9, &c__1, " The number of Ritz values requested is ", (ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&nev, (ftnlen)sizeof(int32_t));
	e_wsle();
	s_wsle(&io___54);
	do_lio(&c__9, &c__1, " The number of Arnoldi vectors generated", (ftnlen)40);
	do_lio(&c__9, &c__1, " (NCV) is ", (ftnlen)10);
	do_lio(&c__3, &c__1, (char *)&ncv, (ftnlen)sizeof(int32_t));
	e_wsle();
	s_wsle(&io___55);
	do_lio(&c__9, &c__1, " What portion of the spectrum: ", (ftnlen)31);
	do_lio(&c__9, &c__1, which, (ftnlen)2);
	e_wsle();
	s_wsle(&io___56);
	do_lio(&c__9, &c__1, " The number of converged Ritz values is ", (ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&nconv, (ftnlen)sizeof(int32_t));
	e_wsle();
	s_wsle(&io___57);
	do_lio(&c__9, &c__1, " The number of Implicit Arnoldi update", (ftnlen)38);
	do_lio(&c__9, &c__1, " iterations taken is ", (ftnlen)21);
	do_lio(&c__3, &c__1, (char *)&iparam[2], (ftnlen)sizeof(int32_t));
	e_wsle();
	s_wsle(&io___58);
	do_lio(&c__9, &c__1, " The number of OP*x is ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&iparam[8], (ftnlen)sizeof(int32_t));
	e_wsle();
	s_wsle(&io___59);
	do_lio(&c__9, &c__1, " The convergence criterion is ", (ftnlen)30);
	do_lio(&c__4, &c__1, (char *)&tol, (ftnlen)sizeof(float));
	e_wsle();
	s_wsle(&io___60);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();

    }

     /* ------------------------- */
     /* Done with program snsimp. */
     /* ------------------------- */

L9000:

    return 0;
} /* MAIN__ */

/* ========================================================================== */

/*     matrix vector subroutine */

/*     The matrix used is the 2 dimensional convection-diffusion */
/*     operator discretized using central difference. */

int snsimp_av_(int32_t *nx, float *v, float *w)
{
    /* System generated locals */
    int32_t i__1;
    float r__1;

    /* Local variables */
    int32_t j;
    float h2;
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
/*     has real eigenvalues.  When rho*h/2 > 1, it has complex */
/*     eigenvalues. */

/*     The subroutine TV is called to computed y<---T*x. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    h2 = 1.f / (float) ((*nx + 1) * (*nx + 1));

	snsimp_tv_(nx, &v[1], &w[1]);
    r__1 = -1.f / h2;
    saxpy_(nx, &r__1, &v[*nx + 1], &c__1, &w[1], &c__1);

    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j) {
	lo = (j - 1) * *nx;
	snsimp_tv_(nx, &v[lo + 1], &w[lo + 1]);
	r__1 = -1.f / h2;
	saxpy_(nx, &r__1, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);
	r__1 = -1.f / h2;
	saxpy_(nx, &r__1, &v[lo + *nx + 1], &c__1, &w[lo + 1], &c__1);
/* L10: */
    }

    lo = (*nx - 1) * *nx;
	snsimp_tv_(nx, &v[lo + 1], &w[lo + 1]);
    r__1 = -1.f / h2;
    saxpy_(nx, &r__1, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);

    return 0;
} /* av_ */

/* ========================================================================= */
int snsimp_tv_(int32_t *nx, float *x, float *y)
{
    /* System generated locals */
    int32_t i__1;

    /* Local variables */
    float h;
    int32_t j;
    float h2, dd, dl, du;

/*     Compute the matrix vector multiplication y<---T*x */
/*     where T is a nx by nx tridiagonal matrix with DD on the */
/*     diagonal, DL on the subdiagonal, and DU on the superdiagonal. */

/*     When rho*h/2 <= 1, the discrete convection-diffusion operator */
/*     has real eigenvalues.  When rho*h/2 > 1, it has complex */
/*     eigenvalues. */

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    h = 1.f / (float) (*nx + 1);
    h2 = h * h;
    dd = 4.f / h2;
    dl = -1.f / h2 - 50.f / h;
    du = -1.f / h2 + 50.f / h;

    y[1] = dd * x[1] + du * x[2];
    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j) {
	y[j] = dl * x[j - 1] + dd * x[j] + du * x[j + 1];
/* L10: */
    }
    y[*nx] = dl * x[*nx - 1] + dd * x[*nx];
    return 0;
} /* tv_ */

