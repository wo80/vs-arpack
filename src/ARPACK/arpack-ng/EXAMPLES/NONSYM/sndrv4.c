/* EXAMPLES\NONSYM\sndrv4.f -- translated by f2c (version 20100827).
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
    real rho;
} convct_;

#define convct_1 convct_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__256 = 256;
static integer c__3 = 3;
static integer c__6 = 6;
static integer c__25 = 25;
static integer c_n6 = -6;
static integer c__4 = 4;

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    real d__[75]	/* was [25][3] */, h__;
    integer j, n;
    real s, v[6400]	/* was [256][25] */, s1, s2, s3, dd[256], dl[256];
    extern /* Subroutine */ int av_(integer *, real *, real *);
    real ax[256], du[256];
    extern /* Subroutine */ int mv_(integer *, real *, real *);
    real mx[256], du2[256];
    integer ido, ncv, nev;
    real tol;
    char bmat[1];
    integer mode, info;
    logical rvec;
    integer ierr, ipiv[256];
    extern doublereal snrm2_(integer *, real *, integer *);
    char which[2];
    real resid[256];
    integer nconv;
    real workd[768];
    logical first;
    integer ipntr[14];
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *);
    real workl[2025];
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *), smout_(integer *, integer *, integer *, real *
	    , integer *, integer *, char *, ftnlen);
    extern doublereal slapy2_(real *, real *);
    integer iparam[11];
    real sigmai;
    logical select[25];
    real sigmar;
    extern /* Subroutine */ int snaupd_(integer *, char *, integer *, char *, 
	    integer *, real *, real *, integer *, real *, integer *, integer *
	    , integer *, real *, real *, integer *, integer *, ftnlen, ftnlen)
	    , sneupd_(logical *, char *, logical *, real *, real *, real *, 
	    integer *, real *, real *, real *, char *, integer *, char *, 
	    integer *, real *, real *, integer *, real *, integer *, integer *
	    , integer *, real *, real *, integer *, integer *, ftnlen, ftnlen,
	     ftnlen);
    integer ishfts, maxitr;
    extern /* Subroutine */ int sgttrf_(integer *, real *, real *, real *, 
	    real *, integer *, integer *);
    integer lworkl;
    real workev[75];
    extern /* Subroutine */ int sgttrs_(char *, integer *, integer *, real *, 
	    real *, real *, real *, integer *, real *, integer *, integer *, 
	    ftnlen);

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 6, 0, 0, 0 };
    static cilist io___5 = { 0, 6, 0, 0, 0 };
    static cilist io___6 = { 0, 6, 0, 0, 0 };
    static cilist io___23 = { 0, 6, 0, 0, 0 };
    static cilist io___24 = { 0, 6, 0, 0, 0 };
    static cilist io___25 = { 0, 6, 0, 0, 0 };
    static cilist io___39 = { 0, 6, 0, 0, 0 };
    static cilist io___40 = { 0, 6, 0, 0, 0 };
    static cilist io___41 = { 0, 6, 0, 0, 0 };
    static cilist io___42 = { 0, 6, 0, 0, 0 };
    static cilist io___43 = { 0, 6, 0, 0, 0 };
    static cilist io___44 = { 0, 6, 0, 0, 0 };
    static cilist io___45 = { 0, 6, 0, 0, 0 };
    static cilist io___46 = { 0, 6, 0, 0, 0 };
    static cilist io___47 = { 0, 6, 0, 0, 0 };
    static cilist io___48 = { 0, 6, 0, 0, 0 };
    static cilist io___53 = { 0, 6, 0, 0, 0 };
    static cilist io___54 = { 0, 6, 0, 0, 0 };
    static cilist io___55 = { 0, 6, 0, 0, 0 };
    static cilist io___56 = { 0, 6, 0, 0, 0 };
    static cilist io___61 = { 0, 6, 0, 0, 0 };
    static cilist io___62 = { 0, 6, 0, 0, 0 };
    static cilist io___63 = { 0, 6, 0, 0, 0 };
    static cilist io___64 = { 0, 6, 0, 0, 0 };
    static cilist io___65 = { 0, 6, 0, 0, 0 };
    static cilist io___66 = { 0, 6, 0, 0, 0 };
    static cilist io___67 = { 0, 6, 0, 0, 0 };
    static cilist io___68 = { 0, 6, 0, 0, 0 };
    static cilist io___69 = { 0, 6, 0, 0, 0 };
    static cilist io___70 = { 0, 6, 0, 0, 0 };
    static cilist io___71 = { 0, 6, 0, 0, 0 };
    static cilist io___72 = { 0, 6, 0, 0, 0 };
    static cilist io___73 = { 0, 6, 0, 0, 0 };
    static cilist io___74 = { 0, 6, 0, 0, 0 };
    static cilist io___75 = { 0, 6, 0, 0, 0 };
    static cilist io___76 = { 0, 6, 0, 0, 0 };
    static cilist io___77 = { 0, 6, 0, 0, 0 };
    static cilist io___78 = { 0, 6, 0, 0, 0 };
    static cilist io___79 = { 0, 6, 0, 0, 0 };



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

/*     ... Use mode 3 of SNAUPD. */

/* \BeginLib */

/* \Routines called: */
/*     snaupd  ARPACK reverse communication interface routine. */
/*     sneupd  ARPACK routine that returns Ritz values and (optionally) */
/*             Ritz vectors. */
/*     sgttrf  LAPACK tridiagonal factorization routine. */
/*     sgttrs  LAPACK tridiagonal linear system solve routine. */
/*     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully. */
/*     saxpy   Level 1 BLAS that computes y <- alpha*x+y. */
/*     scopy   Level 1 BLAS that copies one vector to another. */
/*     sdot    Level 1 BLAS that computes the dot product of two vectors. */
/*     snrm2   Level 1 BLAS that computes the norm of a vector. */
/*     av      Matrix vector multiplication routine that computes A*x. */
/*     mv      Matrix vector multiplication routine that computes M*x. */

/* \Author */
/*     Richard Lehoucq */
/*     Danny Sorensen */
/*     Chao Yang */
/*     Dept. of Computational & */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \SCCS Information: @(#) */
/* FILE: ndrv4.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2 */

/* \Remarks */
/*     1. None */

/* \EndLib */
/* ----------------------------------------------------------------------- */

/*     %-----------------------------% */
/*     | Define leading dimensions   | */
/*     | for all arrays.             | */
/*     | MAXN:   Maximum dimension   | */
/*     |         of the A allowed.   | */
/*     | MAXNEV: Maximum NEV allowed | */
/*     | MAXNCV: Maximum NCV allowed | */
/*     %-----------------------------% */


/*     %--------------% */
/*     | Local Arrays | */
/*     %--------------% */


/*     %---------------% */
/*     | Local Scalars | */
/*     %---------------% */


/*     %-----------------------------% */
/*     | BLAS & LAPACK routines used | */
/*     %-----------------------------% */


/*     %--------------------% */
/*     | Intrinsic function | */
/*     %--------------------% */


/*     %------------% */
/*     | Parameters | */
/*     %------------% */


/*     %-----------------------% */
/*     | Executable statements | */
/*     %-----------------------% */

/*     %----------------------------------------------------% */
/*     | The number N is the dimension of the matrix.  A    | */
/*     | generalized eigenvalue problem is solved (BMAT =   | */
/*     | 'G').  NEV is the number of eigenvalues (closest   | */
/*     | to SIGMAR) to be approximated.  Since the          | */
/*     | shift-invert mode is used,  WHICH is set to 'LM'.  | */
/*     | The user can modify NEV, NCV, SIGMAR to solve      | */
/*     | problems of different sizes, and to get different  | */
/*     | parts of the spectrum.  However, The following     | */
/*     | conditions must be satisfied:                      | */
/*     |                     N <= MAXN,                     | */
/*     |                   NEV <= MAXNEV,                   | */
/*     |               NEV + 2 <= NCV <= MAXNCV             | */
/*     %----------------------------------------------------% */

    n = 100;
    nev = 4;
    ncv = 10;
    if (n > 256) {
	s_wsle(&io___4);
	do_lio(&c__9, &c__1, " ERROR with _NDRV4: N is greater than MAXN ", (
		ftnlen)43);
	e_wsle();
	goto L9000;
    } else if (nev > 10) {
	s_wsle(&io___5);
	do_lio(&c__9, &c__1, " ERROR with _NDRV4: NEV is greater than MAXNEV "
		, (ftnlen)47);
	e_wsle();
	goto L9000;
    } else if (ncv > 25) {
	s_wsle(&io___6);
	do_lio(&c__9, &c__1, " ERROR with _NDRV4: NCV is greater than MAXNCV "
		, (ftnlen)47);
	e_wsle();
	goto L9000;
    }
    *(unsigned char *)bmat = 'G';
    s_copy(which, "LM", (ftnlen)2, (ftnlen)2);
    sigmar = 1.f;
    sigmai = 0.f;

/*     %--------------------------------------------------% */
/*     | Construct C = A - SIGMA*M in real arithmetic,    | */
/*     | and factor C in real arithmetic (using LAPACK    | */
/*     | subroutine sgttrf). The matrix A is chosen to be | */
/*     | the tridiagonal matrix derived from the standard | */
/*     | central difference discretization of the 1-d     | */
/*     | convection-diffusion operator u" + rho*u' on the | */
/*     | interval [0, 1] with zero Dirichlet boundary     | */
/*     | condition.  The matrix M is the mass matrix      | */
/*     | formed by using piecewise linear elements on     | */
/*     | [0,1].                                           | */
/*     %--------------------------------------------------% */

    convct_1.rho = 10.f;
    h__ = 1.f / (real) (n + 1);
    s = convct_1.rho / 2.f;

    s1 = -1.f / h__ - s - sigmar * h__ / 6.f;
    s2 = 2.f / h__ - sigmar * 4.f * h__ / 6.f;
    s3 = -1.f / h__ + s - sigmar * h__ / 6.f;

    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j) {
	dl[j - 1] = s1;
	dd[j - 1] = s2;
	du[j - 1] = s3;
/* L10: */
    }
    dd[n - 1] = s2;

    sgttrf_(&n, dl, dd, du, du2, ipiv, &ierr);
    if (ierr != 0) {
	s_wsle(&io___23);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___24);
	do_lio(&c__9, &c__1, " ERROR with _gttrf in _NDRV4.", (ftnlen)29);
	e_wsle();
	s_wsle(&io___25);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	goto L9000;
    }

/*     %-----------------------------------------------------% */
/*     | The work array WORKL is used in SNAUPD as           | */
/*     | workspace.  Its dimension LWORKL is set as          | */
/*     | illustrated below.  The parameter TOL determines    | */
/*     | the stopping criterion. If TOL<=0, machine          | */
/*     | precision is used.  The variable IDO is used for    | */
/*     | reverse communication, and is initially set to 0.   | */
/*     | Setting INFO=0 indicates that a random vector is    | */
/*     | generated in SNAUPD to start the Arnoldi iteration. | */
/*     %-----------------------------------------------------% */

/* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 * 3 + ncv * 6;
    tol = 0.f;
    ido = 0;
    info = 0;

/*     %---------------------------------------------------% */
/*     | This program uses exact shifts with respect to    | */
/*     | the current Hessenberg matrix (IPARAM(1) = 1).    | */
/*     | IPARAM(3) specifies the maximum number of Arnoldi | */
/*     | iterations allowed.  Mode 3 of SNAUPD is used     | */
/*     | (IPARAM(7) = 3).  All these options can be        | */
/*     | changed by the user. For details, see the         | */
/*     | documentation in SNAUPD.                          | */
/*     %---------------------------------------------------% */

    ishfts = 1;
    maxitr = 300;
    mode = 3;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

/*     %------------------------------------------% */
/*     | M A I N   L O O P(Reverse communication) | */
/*     %------------------------------------------% */

L20:

/*        %---------------------------------------------% */
/*        | Repeatedly call the routine SNAUPD and take | */
/*        | actions indicated by parameter IDO until    | */
/*        | either convergence is indicated or maxitr   | */
/*        | has been exceeded.                          | */
/*        %---------------------------------------------% */

    snaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, 
	    iparam, ipntr, workd, workl, &lworkl, &info, (ftnlen)1, (ftnlen)2)
	    ;

    if (ido == -1) {

/*           %-------------------------------------------% */
/*           | Perform  y <--- OP*x = inv[A-SIGMA*M]*M*x | */
/*           | to force starting vector into the range   | */
/*           | of OP.   The user should supply his/her   | */
/*           | own matrix vector multiplication routine  | */
/*           | and a linear system solver.  The matrix   | */
/*           | vector multiplication routine should take | */
/*           | workd(ipntr(1)) as the input. The final   | */
/*           | result should be returned to              | */
/*           | workd(ipntr(2)).                          | */
/*           %-------------------------------------------% */

	mv_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
	sgttrs_("N", &n, &c__1, dl, dd, du, du2, ipiv, &workd[ipntr[1] - 1], &
		n, &ierr, (ftnlen)1);
	if (ierr != 0) {
	    s_wsle(&io___39);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___40);
	    do_lio(&c__9, &c__1, " ERROR with _gttrs in _NDRV4.", (ftnlen)29);
	    e_wsle();
	    s_wsle(&io___41);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    goto L9000;
	}

/*           %-----------------------------------------% */
/*           | L O O P   B A C K to call SNAUPD again. | */
/*           %-----------------------------------------% */

	goto L20;

    } else if (ido == 1) {

/*           %-----------------------------------------% */
/*           | Perform y <-- OP*x = inv[A-sigma*M]*M*x | */
/*           | M*x has been saved in workd(ipntr(3)).  | */
/*           | The user only need the linear system    | */
/*           | solver here that takes workd(ipntr(3))  | */
/*           | as input, and returns the result to     | */
/*           | workd(ipntr(2)).                        | */
/*           %-----------------------------------------% */

	scopy_(&n, &workd[ipntr[2] - 1], &c__1, &workd[ipntr[1] - 1], &c__1);
	sgttrs_("N", &n, &c__1, dl, dd, du, du2, ipiv, &workd[ipntr[1] - 1], &
		n, &ierr, (ftnlen)1);
	if (ierr != 0) {
	    s_wsle(&io___42);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___43);
	    do_lio(&c__9, &c__1, " ERROR with _gttrs in _NDRV4.", (ftnlen)29);
	    e_wsle();
	    s_wsle(&io___44);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    goto L9000;
	}

/*           %-----------------------------------------% */
/*           | L O O P   B A C K to call SNAUPD again. | */
/*           %-----------------------------------------% */

	goto L20;

    } else if (ido == 2) {

/*           %---------------------------------------------% */
/*           |          Perform  y <--- M*x                | */
/*           | Need matrix vector multiplication routine   | */
/*           | here that takes workd(ipntr(1)) as input    | */
/*           | and returns the result to workd(ipntr(2)).  | */
/*           %---------------------------------------------% */

	mv_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

/*           %-----------------------------------------% */
/*           | L O O P   B A C K to call SNAUPD again. | */
/*           %-----------------------------------------% */

	goto L20;

    }


/*     %-----------------------------------------% */
/*     | Either we have convergence, or there is | */
/*     | an error.                               | */
/*     %-----------------------------------------% */

    if (info < 0) {

/*        %--------------------------% */
/*        | Error message, check the | */
/*        | documentation in SNAUPD. | */
/*        %--------------------------% */

	s_wsle(&io___45);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___46);
	do_lio(&c__9, &c__1, " Error with _naupd, info = ", (ftnlen)27);
	do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___47);
	do_lio(&c__9, &c__1, " Check the documentation in _naupd.", (ftnlen)
		35);
	e_wsle();
	s_wsle(&io___48);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();

    } else {

/*        %-------------------------------------------% */
/*        | No fatal errors occurred.                 | */
/*        | Post-Process using SNEUPD.                | */
/*        |                                           | */
/*        | Computed eigenvalues may be extracted.    | */
/*        |                                           | */
/*        | Eigenvectors may also be computed now if  | */
/*        | desired.  (indicated by rvec = .true.)    | */
/*        %-------------------------------------------% */

	rvec = TRUE_;
	sneupd_(&rvec, "A", select, d__, &d__[25], v, &c__256, &sigmar, &
		sigmai, workev, bmat, &n, which, &nev, &tol, resid, &ncv, v, &
		c__256, iparam, ipntr, workd, workl, &lworkl, &ierr, (ftnlen)
		1, (ftnlen)1, (ftnlen)2);

/*        %-----------------------------------------------% */
/*        | The real part of the eigenvalue is returned   | */
/*        | in the first column of the two dimensional    | */
/*        | array D, and the IMAGINARY part is returned   | */
/*        | in the second column of D.  The corresponding | */
/*        | eigenvectors are returned in the first NEV    | */
/*        | columns of the two dimensional array V if     | */
/*        | requested.  Otherwise, an orthogonal basis    | */
/*        | for the invariant subspace corresponding to   | */
/*        | the eigenvalues in D is returned in V.        | */
/*        %-----------------------------------------------% */

	if (ierr != 0) {

/*            %------------------------------------% */
/*            | Error condition:                   | */
/*            | Check the documentation of SNEUPD. | */
/*            %------------------------------------% */

	    s_wsle(&io___53);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___54);
	    do_lio(&c__9, &c__1, " Error with _neupd, info = ", (ftnlen)27);
	    do_lio(&c__3, &c__1, (char *)&ierr, (ftnlen)sizeof(integer));
	    e_wsle();
	    s_wsle(&io___55);
	    do_lio(&c__9, &c__1, " Check the documentation of _neupd. ", (
		    ftnlen)36);
	    e_wsle();
	    s_wsle(&io___56);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();

	} else {

	    first = TRUE_;
	    nconv = iparam[4];
	    i__1 = nconv;
	    for (j = 1; j <= i__1; ++j) {

/*                %---------------------------% */
/*                | Compute the residual norm | */
/*                |                           | */
/*                |   ||  A*x - lambda*x ||   | */
/*                |                           | */
/*                | for the NCONV accurately  | */
/*                | computed eigenvalues and  | */
/*                | eigenvectors.  (iparam(5) | */
/*                | indicates how many are    | */
/*                | accurate to the requested | */
/*                | tolerance)                | */
/*                %---------------------------% */

		if (d__[j + 24] == 0.f) {

/*                    %--------------------% */
/*                    | Ritz value is real | */
/*                    %--------------------% */

		    av_(&n, &v[(j << 8) - 256], ax);
		    mv_(&n, &v[(j << 8) - 256], mx);
		    r__1 = -d__[j - 1];
		    saxpy_(&n, &r__1, mx, &c__1, ax, &c__1);
		    d__[j + 49] = snrm2_(&n, ax, &c__1);
		    d__[j + 49] /= (r__1 = d__[j - 1], dabs(r__1));

		} else if (first) {

/*                    %------------------------% */
/*                    | Ritz value is complex. | */
/*                    | Residual of one Ritz   | */
/*                    | value of the conjugate | */
/*                    | pair is computed.      | */
/*                    %------------------------% */

		    av_(&n, &v[(j << 8) - 256], ax);
		    mv_(&n, &v[(j << 8) - 256], mx);
		    r__1 = -d__[j - 1];
		    saxpy_(&n, &r__1, mx, &c__1, ax, &c__1);
		    mv_(&n, &v[(j + 1 << 8) - 256], mx);
		    saxpy_(&n, &d__[j + 24], mx, &c__1, ax, &c__1);
		    d__[j + 49] = snrm2_(&n, ax, &c__1);
		    av_(&n, &v[(j + 1 << 8) - 256], ax);
		    mv_(&n, &v[(j + 1 << 8) - 256], mx);
		    r__1 = -d__[j - 1];
		    saxpy_(&n, &r__1, mx, &c__1, ax, &c__1);
		    mv_(&n, &v[(j << 8) - 256], mx);
		    r__1 = -d__[j + 24];
		    saxpy_(&n, &r__1, mx, &c__1, ax, &c__1);
		    r__1 = snrm2_(&n, ax, &c__1);
		    d__[j + 49] = slapy2_(&d__[j + 49], &r__1);
		    d__[j + 49] /= slapy2_(&d__[j - 1], &d__[j + 24]);
		    d__[j + 50] = d__[j + 49];
		    first = FALSE_;
		} else {
		    first = TRUE_;
		}

/* L30: */
	    }

/*            %-----------------------------% */
/*            | Display computed residuals. | */
/*            %-----------------------------% */

	    smout_(&c__6, &nconv, &c__3, d__, &c__25, &c_n6, "Ritz values (R"
		    "eal,Imag) and relative residuals", (ftnlen)46);

	}

/*        %-------------------------------------------% */
/*        | Print additional convergence information. | */
/*        %-------------------------------------------% */

	if (info == 1) {
	    s_wsle(&io___61);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___62);
	    do_lio(&c__9, &c__1, " Maximum number of iterations reached.", (
		    ftnlen)38);
	    e_wsle();
	    s_wsle(&io___63);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	} else if (info == 3) {
	    s_wsle(&io___64);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___65);
	    do_lio(&c__9, &c__1, " No shifts could be applied during implicit"
		    , (ftnlen)43);
	    do_lio(&c__9, &c__1, " Arnoldi update, try increasing NCV.", (
		    ftnlen)36);
	    e_wsle();
	    s_wsle(&io___66);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	}

	s_wsle(&io___67);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___68);
	do_lio(&c__9, &c__1, " _NDRV4 ", (ftnlen)8);
	e_wsle();
	s_wsle(&io___69);
	do_lio(&c__9, &c__1, " ====== ", (ftnlen)8);
	e_wsle();
	s_wsle(&io___70);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___71);
	do_lio(&c__9, &c__1, " Size of the matrix is ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___72);
	do_lio(&c__9, &c__1, " The number of Ritz values requested is ", (
		ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&nev, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___73);
	do_lio(&c__9, &c__1, " The number of Arnoldi vectors generated", (
		ftnlen)40);
	do_lio(&c__9, &c__1, " (NCV) is ", (ftnlen)10);
	do_lio(&c__3, &c__1, (char *)&ncv, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___74);
	do_lio(&c__9, &c__1, " What portion of the spectrum: ", (ftnlen)31);
	do_lio(&c__9, &c__1, which, (ftnlen)2);
	e_wsle();
	s_wsle(&io___75);
	do_lio(&c__9, &c__1, " The number of converged Ritz values is ", (
		ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&nconv, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___76);
	do_lio(&c__9, &c__1, " The number of Implicit Arnoldi update", (
		ftnlen)38);
	do_lio(&c__9, &c__1, " iterations taken is ", (ftnlen)21);
	do_lio(&c__3, &c__1, (char *)&iparam[2], (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___77);
	do_lio(&c__9, &c__1, " The number of OP*x is ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&iparam[8], (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___78);
	do_lio(&c__9, &c__1, " The convergence criterion is ", (ftnlen)30);
	do_lio(&c__4, &c__1, (char *)&tol, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___79);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();

    }

/*     %---------------------------% */
/*     | Done with program sndrv4. | */
/*     %---------------------------% */

L9000:

    return 0;
} /* MAIN__ */


/* ========================================================================== */

/*     matrix vector multiplication subroutine */

/* Subroutine */ int mv_(integer *n, real *v, real *w)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    real h__;
    integer j;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);


/*     Compute the matrix vector multiplication y<---M*x */
/*     where M is mass matrix formed by using piecewise linear elements */
/*     on [0,1]. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    w[1] = (v[1] * 4.f + v[2] * 1.f) / 6.f;
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j) {
	w[j] = (v[j - 1] * 1.f + v[j] * 4.f + v[j + 1] * 1.f) / 6.f;
/* L10: */
    }
    w[*n] = (v[*n - 1] * 1.f + v[*n] * 4.f) / 6.f;

    h__ = 1.f / (real) (*n + 1);
    sscal_(n, &h__, &w[1], &c__1);
    return 0;
} /* mv_ */

/* ------------------------------------------------------------------ */
/* Subroutine */ int av_(integer *n, real *v, real *w)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    real h__;
    integer j;
    real s, dd, dl, du;


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
    h__ = 1.f / (real) (*n + 1);
    s = convct_1.rho / 2.f;
    dd = 2.f / h__;
    dl = -1.f / h__ - s;
    du = -1.f / h__ + s;

    w[1] = dd * v[1] + du * v[2];
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j) {
	w[j] = dl * v[j - 1] + dd * v[j] + du * v[j + 1];
/* L10: */
    }
    w[*n] = dl * v[*n - 1] + dd * v[*n];
    return 0;
} /* av_ */

/* Main program alias */ int sndrv4_ () { MAIN__ (); return 0; }