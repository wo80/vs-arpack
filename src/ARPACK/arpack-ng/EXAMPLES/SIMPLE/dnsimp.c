/* EXAMPLES\SIMPLE\dnsimp.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    integer logfil, ndigit, mgetv0, msaupd, msaup2, msaitr, mseigt, msapps, 
	    msgets, mseupd, mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, 
	    mneupd, mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd;
} debug_;

#define debug_1 debug_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__256 = 256;
static integer c__3 = 3;
static integer c__6 = 6;
static integer c__30 = 30;
static integer c_n6 = -6;
static integer c__5 = 5;

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    doublereal d__[90]	/* was [30][3] */;
    integer j, n;
    doublereal v[7680]	/* was [256][30] */;
    extern /* Subroutine */ int av_(integer *, doublereal *, doublereal *);
    doublereal ax[256];
    integer nx, ido, ncv, nev;
    doublereal tol;
    char bmat[1];
    integer info;
    logical rvec;
    integer ierr, mode1;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    char which[2];
    doublereal resid[256];
    integer nconv;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    doublereal workd[768];
    logical first;
    extern /* Subroutine */ int dmout_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, char *, ftnlen);
    integer ipntr[14];
    doublereal workl[2880];
    extern doublereal dlapy2_(doublereal *, doublereal *);
    integer iparam[11];
    doublereal sigmai;
    logical select[30];
    extern /* Subroutine */ int dnaupd_(integer *, char *, integer *, char *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen);
    doublereal sigmar;
    extern /* Subroutine */ int dneupd_(logical *, char *, logical *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, char *, integer *, char *, integer *,
	     doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    integer ishfts, maxitr, lworkl;
    doublereal workev[90];

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
/*           needed to utilize the top level ARPACK routine DNAUPD */
/*           that computes the quantities needed to construct */
/*           the desired eigenvalues and eigenvectors(if requested). */

/*        3) How to extract the desired eigenvalues and eigenvectors */
/*           using the ARPACK routine DNEUPD. */

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
/*     ... Use mode 1 of DNAUPD. */

/* \BeginLib */

/* \Routines called: */
/*     dnaupd  ARPACK reverse communication interface routine. */
/*     dneupd  ARPACK routine that returns Ritz values and (optionally) */
/*             Ritz vectors. */
/*     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully. */
/*     daxpy   Level 1 BLAS that computes y <- alpha*x+y. */
/*     dnrm2   Level 1 BLAS that computes the norm of a vector. */
/*     av      Matrix vector multiplication routine that computes A*x. */
/*     tv      Matrix vector multiplication routine that computes T*x, */
/*             where T is a tridiagonal matrix.  It is used in routine */
/*             av. */

/* \Author */
/*     Richard Lehoucq */
/*     Danny Sorensen */
/*     Chao Yang */
/*     Dept. of Computational & */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \SCCS Information: @(#) */
/* FILE: nsimp.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2 */

/* \Remarks */
/*     1. None */

/* \EndLib */
/* --------------------------------------------------------------------------- */

/*     %------------------------------------------------------% */
/*     | Storage Declarations:                                | */
/*     |                                                      | */
/*     | The maximum dimensions for all arrays are            | */
/*     | set here to accommodate a problem size of            | */
/*     | N .le. MAXN                                          | */
/*     |                                                      | */
/*     | NEV is the number of eigenvalues requested.          | */
/*     |     See specifications for ARPACK usage below.       | */
/*     |                                                      | */
/*     | NCV is the largest number of basis vectors that will | */
/*     |     be used in the Implicitly Restarted Arnoldi      | */
/*     |     Process.  Work per major iteration is            | */
/*     |     proportional to N*NCV*NCV.                       | */
/*     |                                                      | */
/*     | You must set:                                        | */
/*     |                                                      | */
/*     | MAXN:   Maximum dimension of the A allowed.          | */
/*     | MAXNEV: Maximum NEV allowed.                         | */
/*     | MAXNCV: Maximum NCV allowed.                         | */
/*     %------------------------------------------------------% */


/*     %--------------% */
/*     | Local Arrays | */
/*     %--------------% */


/*     %---------------% */
/*     | Local Scalars | */
/*     %---------------% */


/*     %------------% */
/*     | Parameters | */
/*     %------------% */


/*     %-----------------------------% */
/*     | BLAS & LAPACK routines used | */
/*     %-----------------------------% */


/*     %--------------------% */
/*     | Intrinsic function | */
/*     %--------------------% */


/*     %-----------------------% */
/*     | Executable Statements | */
/*     %-----------------------% */

/*     %-------------------------------------------------% */
/*     | The following include statement and assignments | */
/*     | initiate trace output from the internal         | */
/*     | actions of ARPACK.  See debug.doc in the        | */
/*     | DOCUMENTS directory for usage.  Initially, the  | */
/*     | most useful information will be a breakdown of  | */
/*     | time spent in the various stages of computation | */
/*     | given by setting mnaupd = 1.                    | */
/*     %-------------------------------------------------% */


/* \SCCS Information: @(#) */
/* FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 */

/*     %---------------------------------% */
/*     | See debug.doc for documentation | */
/*     %---------------------------------% */
    debug_1.ndigit = -3;
    debug_1.logfil = 6;
    debug_1.mnaitr = 0;
    debug_1.mnapps = 0;
    debug_1.mnaupd = 1;
    debug_1.mnaup2 = 0;
    debug_1.mneigh = 0;
    debug_1.mneupd = 0;

/*     %-------------------------------------------------% */
/*     | The following sets dimensions for this problem. | */
/*     %-------------------------------------------------% */

    nx = 10;
    n = nx * nx;

/*     %-----------------------------------------------% */
/*     |                                               | */
/*     | Specifications for ARPACK usage are set       | */
/*     | below:                                        | */
/*     |                                               | */
/*     |    1) NEV = 4  asks for 4 eigenvalues to be   | */
/*     |       computed.                               | */
/*     |                                               | */
/*     |    2) NCV = 20 sets the length of the Arnoldi | */
/*     |       factorization.                          | */
/*     |                                               | */
/*     |    3) This is a standard problem.             | */
/*     |         (indicated by bmat  = 'I')            | */
/*     |                                               | */
/*     |    4) Ask for the NEV eigenvalues of          | */
/*     |       largest magnitude.                      | */
/*     |         (indicated by which = 'LM')           | */
/*     |       See documentation in DNAUPD for the     | */
/*     |       other options SM, LR, SR, LI, SI.       | */
/*     |                                               | */
/*     | Note: NEV and NCV must satisfy the following  | */
/*     | conditions:                                   | */
/*     |              NEV <= MAXNEV                    | */
/*     |          NEV + 2 <= NCV <= MAXNCV             | */
/*     |                                               | */
/*     %-----------------------------------------------% */

    nev = 4;
    ncv = 20;
    *(unsigned char *)bmat = 'I';
    s_copy(which, "LM", (ftnlen)2, (ftnlen)2);

    if (n > 256) {
	s_wsle(&io___7);
	do_lio(&c__9, &c__1, " ERROR with _NSIMP: N is greater than MAXN ", (
		ftnlen)43);
	e_wsle();
	goto L9000;
    } else if (nev > 12) {
	s_wsle(&io___8);
	do_lio(&c__9, &c__1, " ERROR with _NSIMP: NEV is greater than MAXNEV "
		, (ftnlen)47);
	e_wsle();
	goto L9000;
    } else if (ncv > 30) {
	s_wsle(&io___9);
	do_lio(&c__9, &c__1, " ERROR with _NSIMP: NCV is greater than MAXNCV "
		, (ftnlen)47);
	e_wsle();
	goto L9000;
    }

/*     %-----------------------------------------------------% */
/*     |                                                     | */
/*     | Specification of stopping rules and initial         | */
/*     | conditions before calling DNAUPD                    | */
/*     |                                                     | */
/*     | TOL  determines the stopping criterion.             | */
/*     |                                                     | */
/*     |      Expect                                         | */
/*     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) | */
/*     |               computed   true                       | */
/*     |                                                     | */
/*     |      If TOL .le. 0,  then TOL <- macheps            | */
/*     |           (machine precision) is used.              | */
/*     |                                                     | */
/*     | IDO  is the REVERSE COMMUNICATION parameter         | */
/*     |      used to specify actions to be taken on return  | */
/*     |      from DNAUPD. (see usage below)                 | */
/*     |                                                     | */
/*     |      It MUST initially be set to 0 before the first | */
/*     |      call to DNAUPD.                                | */
/*     |                                                     | */
/*     | INFO on entry specifies starting vector information | */
/*     |      and on return indicates error codes            | */
/*     |                                                     | */
/*     |      Initially, setting INFO=0 indicates that a     | */
/*     |      random starting vector is requested to         | */
/*     |      start the ARNOLDI iteration.  Setting INFO to  | */
/*     |      a nonzero value on the initial call is used    | */
/*     |      if you want to specify your own starting       | */
/*     |      vector (This vector must be placed in RESID).  | */
/*     |                                                     | */
/*     | The work array WORKL is used in DNAUPD as           | */
/*     | workspace.  Its dimension LWORKL is set as          | */
/*     | illustrated below.                                  | */
/*     |                                                     | */
/*     %-----------------------------------------------------% */

/* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 * 3 + ncv * 6;
    tol = 0.;
    ido = 0;
    info = 0;

/*     %---------------------------------------------------% */
/*     | Specification of Algorithm Mode:                  | */
/*     |                                                   | */
/*     | This program uses the exact shift strategy        | */
/*     | (indicated by setting IPARAM(1) = 1).             | */
/*     | IPARAM(3) specifies the maximum number of Arnoldi | */
/*     | iterations allowed.  Mode 1 of DNAUPD is used     | */
/*     | (IPARAM(7) = 1). All these options can be changed | */
/*     | by the user. For details see the documentation in | */
/*     | DNAUPD.                                           | */
/*     %---------------------------------------------------% */

    ishfts = 1;
    maxitr = 300;
    mode1 = 1;

    iparam[0] = ishfts;

    iparam[2] = maxitr;

    iparam[6] = mode1;

/*     %-------------------------------------------% */
/*     | M A I N   L O O P (Reverse communication) | */
/*     %-------------------------------------------% */

L10:

/*        %---------------------------------------------% */
/*        | Repeatedly call the routine DNAUPD and take | */
/*        | actions indicated by parameter IDO until    | */
/*        | either convergence is indicated or maxitr   | */
/*        | has been exceeded.                          | */
/*        %---------------------------------------------% */

    dnaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, 
	    iparam, ipntr, workd, workl, &lworkl, &info, (ftnlen)1, (ftnlen)2)
	    ;

    if (ido == -1 || ido == 1) {

/*           %-------------------------------------------% */
/*           | Perform matrix vector multiplication      | */
/*           |                y <--- Op*x                | */
/*           | The user should supply his/her own        | */
/*           | matrix vector multiplication routine here | */
/*           | that takes workd(ipntr(1)) as the input   | */
/*           | vector, and return the matrix vector      | */
/*           | product to workd(ipntr(2)).               | */
/*           %-------------------------------------------% */

	av_(&nx, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

/*           %-----------------------------------------% */
/*           | L O O P   B A C K to call DNAUPD again. | */
/*           %-----------------------------------------% */

	goto L10;

    }

/*     %----------------------------------------% */
/*     | Either we have convergence or there is | */
/*     | an error.                              | */
/*     %----------------------------------------% */

    if (info < 0) {

/*        %--------------------------% */
/*        | Error message, check the | */
/*        | documentation in DNAUPD. | */
/*        %--------------------------% */

	s_wsle(&io___23);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___24);
	do_lio(&c__9, &c__1, " Error with _naupd, info = ", (ftnlen)27);
	do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___25);
	do_lio(&c__9, &c__1, " Check the documentation of _naupd", (ftnlen)34)
		;
	e_wsle();
	s_wsle(&io___26);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();

    } else {

/*        %-------------------------------------------% */
/*        | No fatal errors occurred.                 | */
/*        | Post-Process using DNEUPD.                | */
/*        |                                           | */
/*        | Computed eigenvalues may be extracted.    | */
/*        |                                           | */
/*        | Eigenvectors may be also computed now if  | */
/*        | desired.  (indicated by rvec = .true.)    | */
/*        |                                           | */
/*        | The routine DNEUPD now called to do this  | */
/*        | post processing (Other modes may require  | */
/*        | more complicated post processing than     | */
/*        | mode1,)                                   | */
/*        |                                           | */
/*        %-------------------------------------------% */

	rvec = TRUE_;

	dneupd_(&rvec, "A", select, d__, &d__[30], v, &c__256, &sigmar, &
		sigmai, workev, bmat, &n, which, &nev, &tol, resid, &ncv, v, &
		c__256, iparam, ipntr, workd, workl, &lworkl, &ierr, (ftnlen)
		1, (ftnlen)1, (ftnlen)2);

/*        %------------------------------------------------% */
/*        | The real parts of the eigenvalues are returned | */
/*        | in the first column of the two dimensional     | */
/*        | array D, and the IMAGINARY part are returned   | */
/*        | in the second column of D.  The corresponding  | */
/*        | eigenvectors are returned in the first         | */
/*        | NCONV (= IPARAM(5)) columns of the two         | */
/*        | dimensional array V if requested.  Otherwise,  | */
/*        | an orthogonal basis for the invariant subspace | */
/*        | corresponding to the eigenvalues in D is       | */
/*        | returned in V.                                 | */
/*        %------------------------------------------------% */

	if (ierr != 0) {

/*           %------------------------------------% */
/*           | Error condition:                   | */
/*           | Check the documentation of DNEUPD. | */
/*           %------------------------------------% */

	    s_wsle(&io___34);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___35);
	    do_lio(&c__9, &c__1, " Error with _neupd, info = ", (ftnlen)27);
	    do_lio(&c__3, &c__1, (char *)&ierr, (ftnlen)sizeof(integer));
	    e_wsle();
	    s_wsle(&io___36);
	    do_lio(&c__9, &c__1, " Check the documentation of _neupd. ", (
		    ftnlen)36);
	    e_wsle();
	    s_wsle(&io___37);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();

	} else {

	    first = TRUE_;
	    nconv = iparam[4];
	    i__1 = nconv;
	    for (j = 1; j <= i__1; ++j) {

/*              %---------------------------% */
/*              | Compute the residual norm | */
/*              |                           | */
/*              |   ||  A*x - lambda*x ||   | */
/*              |                           | */
/*              | for the NCONV accurately  | */
/*              | computed eigenvalues and  | */
/*              | eigenvectors.  (IPARAM(5) | */
/*              | indicates how many are    | */
/*              | accurate to the requested | */
/*              | tolerance)                | */
/*              %---------------------------% */

		if (d__[j + 29] == 0.) {

/*                 %--------------------% */
/*                 | Ritz value is real | */
/*                 %--------------------% */

		    av_(&nx, &v[(j << 8) - 256], ax);
		    d__1 = -d__[j - 1];
		    daxpy_(&n, &d__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
		    d__[j + 59] = dnrm2_(&n, ax, &c__1);
		    d__[j + 59] /= (d__1 = d__[j - 1], abs(d__1));

		} else if (first) {

/*                 %------------------------% */
/*                 | Ritz value is complex. | */
/*                 | Residual of one Ritz   | */
/*                 | value of the conjugate | */
/*                 | pair is computed.      | */
/*                 %------------------------% */

		    av_(&nx, &v[(j << 8) - 256], ax);
		    d__1 = -d__[j - 1];
		    daxpy_(&n, &d__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
		    daxpy_(&n, &d__[j + 29], &v[(j + 1 << 8) - 256], &c__1, 
			    ax, &c__1);
		    d__[j + 59] = dnrm2_(&n, ax, &c__1);
		    av_(&nx, &v[(j + 1 << 8) - 256], ax);
		    d__1 = -d__[j + 29];
		    daxpy_(&n, &d__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
		    d__1 = -d__[j - 1];
		    daxpy_(&n, &d__1, &v[(j + 1 << 8) - 256], &c__1, ax, &
			    c__1);
		    d__1 = dnrm2_(&n, ax, &c__1);
		    d__[j + 59] = dlapy2_(&d__[j + 59], &d__1);
		    d__[j + 59] /= dlapy2_(&d__[j - 1], &d__[j + 29]);
		    d__[j + 60] = d__[j + 59];
		    first = FALSE_;
		} else {
		    first = TRUE_;
		}

/* L20: */
	    }

/*           %-----------------------------% */
/*           | Display computed residuals. | */
/*           %-----------------------------% */

	    dmout_(&c__6, &nconv, &c__3, d__, &c__30, &c_n6, "Ritz values (R"
		    "eal, Imag) and residual residuals", (ftnlen)47);
	}

/*        %-------------------------------------------% */
/*        | Print additional convergence information. | */
/*        %-------------------------------------------% */

	if (info == 1) {
	    s_wsle(&io___42);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___43);
	    do_lio(&c__9, &c__1, " Maximum number of iterations reached.", (
		    ftnlen)38);
	    e_wsle();
	    s_wsle(&io___44);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	} else if (info == 3) {
	    s_wsle(&io___45);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___46);
	    do_lio(&c__9, &c__1, " No shifts could be applied during implicit"
		    , (ftnlen)43);
	    do_lio(&c__9, &c__1, " Arnoldi update, try increasing NCV.", (
		    ftnlen)36);
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
	do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___53);
	do_lio(&c__9, &c__1, " The number of Ritz values requested is ", (
		ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&nev, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___54);
	do_lio(&c__9, &c__1, " The number of Arnoldi vectors generated", (
		ftnlen)40);
	do_lio(&c__9, &c__1, " (NCV) is ", (ftnlen)10);
	do_lio(&c__3, &c__1, (char *)&ncv, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___55);
	do_lio(&c__9, &c__1, " What portion of the spectrum: ", (ftnlen)31);
	do_lio(&c__9, &c__1, which, (ftnlen)2);
	e_wsle();
	s_wsle(&io___56);
	do_lio(&c__9, &c__1, " The number of converged Ritz values is ", (
		ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&nconv, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___57);
	do_lio(&c__9, &c__1, " The number of Implicit Arnoldi update", (
		ftnlen)38);
	do_lio(&c__9, &c__1, " iterations taken is ", (ftnlen)21);
	do_lio(&c__3, &c__1, (char *)&iparam[2], (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___58);
	do_lio(&c__9, &c__1, " The number of OP*x is ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&iparam[8], (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___59);
	do_lio(&c__9, &c__1, " The convergence criterion is ", (ftnlen)30);
	do_lio(&c__5, &c__1, (char *)&tol, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___60);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();

    }

/*     %---------------------------% */
/*     | Done with program dnsimp. | */
/*     %---------------------------% */

L9000:

    return 0;
} /* MAIN__ */


/* ========================================================================== */

/*     matrix vector subroutine */

/*     The matrix used is the 2 dimensional convection-diffusion */
/*     operator discretized using central difference. */

/* Subroutine */ int av_(integer *nx, doublereal *v, doublereal *w)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    integer j;
    doublereal h2;
    integer lo;
    extern /* Subroutine */ int tv_(integer *, doublereal *, doublereal *), 
	    daxpy_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *);


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
    h2 = 1. / (doublereal) ((*nx + 1) * (*nx + 1));

    tv_(nx, &v[1], &w[1]);
    d__1 = -1. / h2;
    daxpy_(nx, &d__1, &v[*nx + 1], &c__1, &w[1], &c__1);

    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j) {
	lo = (j - 1) * *nx;
	tv_(nx, &v[lo + 1], &w[lo + 1]);
	d__1 = -1. / h2;
	daxpy_(nx, &d__1, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);
	d__1 = -1. / h2;
	daxpy_(nx, &d__1, &v[lo + *nx + 1], &c__1, &w[lo + 1], &c__1);
/* L10: */
    }

    lo = (*nx - 1) * *nx;
    tv_(nx, &v[lo + 1], &w[lo + 1]);
    d__1 = -1. / h2;
    daxpy_(nx, &d__1, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);

    return 0;
} /* av_ */

/* ========================================================================= */
/* Subroutine */ int tv_(integer *nx, doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    doublereal h__;
    integer j;
    doublereal h2, dd, dl, du;




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
    h__ = 1. / (doublereal) (*nx + 1);
    h2 = h__ * h__;
    dd = 4. / h2;
    dl = -1. / h2 - 50. / h__;
    du = -1. / h2 + 50. / h__;

    y[1] = dd * x[1] + du * x[2];
    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j) {
	y[j] = dl * x[j - 1] + dd * x[j] + du * x[j + 1];
/* L10: */
    }
    y[*nx] = dl * x[*nx - 1] + dd * x[*nx];
    return 0;
} /* tv_ */

/* Main program alias */ int dnsimp_ () { MAIN__ (); return 0; }
