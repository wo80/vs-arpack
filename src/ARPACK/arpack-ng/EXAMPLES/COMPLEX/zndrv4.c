/* EXAMPLES\COMPLEX\zndrv4.f -- translated by f2c (version 20100827).
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
    doublecomplex rho;
} convct_;

#define convct_1 convct_

/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b3 = {2.,0.};
static doublecomplex c_b5 = {6.,0.};
static integer c__9 = 9;
static integer c__1 = 1;
static integer c__256 = 256;
static integer c__3 = 3;
static integer c__6 = 6;
static integer c__25 = 25;
static integer c_n6 = -6;
static integer c__5 = 5;

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
    double d_imag(doublecomplex *);

    /* Local variables */
    doublecomplex d__[25], h__;
    integer j, n;
    doublecomplex s, v[6400]	/* was [256][25] */, s1, s2, s3, dd[256], dl[
	    256];
    doublereal rd[75]	/* was [25][3] */;
    extern /* Subroutine */ int av_(integer *, doublecomplex *, doublecomplex 
	    *);
    doublecomplex ax[256], du[256];
    extern /* Subroutine */ int mv_(integer *, doublecomplex *, doublecomplex 
	    *);
    doublecomplex mx[256], du2[256];
    integer ido, ncv, nev;
    doublereal tol;
    char bmat[1];
    integer mode, info;
    logical rvec;
    integer ierr, ipiv[256];
    doublecomplex sigma;
    char which[2];
    doublecomplex resid[256];
    integer nconv;
    doublecomplex workd[768];
    extern /* Subroutine */ int dmout_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, char *, ftnlen);
    integer ipntr[14];
    doublecomplex workl[2000];
    doublereal rwork[256];
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern doublereal dlapy2_(doublereal *, doublereal *), dznrm2_(integer *, 
	    doublecomplex *, integer *);
    integer iparam[11];
    logical select[25];
    integer ishfts;
    extern /* Subroutine */ int znaupd_(integer *, char *, integer *, char *, 
	    integer *, doublereal *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen);
    integer maxitr;
    extern /* Subroutine */ int zneupd_(logical *, char *, logical *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, char *, integer *, char *, integer *, doublereal 
	    *, doublecomplex *, integer *, doublecomplex *, integer *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    integer lworkl;
    doublecomplex workev[50];
    extern /* Subroutine */ int zgttrf_(integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    integer *), zgttrs_(char *, integer *, integer *, doublecomplex *,
	     doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 6, 0, 0, 0 };
    static cilist io___5 = { 0, 6, 0, 0, 0 };
    static cilist io___6 = { 0, 6, 0, 0, 0 };
    static cilist io___22 = { 0, 6, 0, 0, 0 };
    static cilist io___23 = { 0, 6, 0, 0, 0 };
    static cilist io___24 = { 0, 6, 0, 0, 0 };
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
/*     in shift and invert mode for a generalized complex nonsymmetric */
/*     eigenvalue problem. */

/*     We implement example four of ex-complex.doc in DOCUMENTS directory */

/* \Example-4 */
/*     ... Suppose we want to solve A*x = lambda*B*x in shift-invert mode, */
/*         where A and B are derived from a finite element discretization */
/*         of a 1-dimensional convection-diffusion operator */
/*                         (d^2u/dx^2) + rho*(du/dx) */
/*         on the interval [0,1] with zero boundary condition using */
/*         piecewise linear elements. */

/*     ... where the shift sigma is a complex number. */

/*     ... OP = inv[A-SIGMA*M]*M  and  B = M. */

/*     ... Use mode 3 of ZNAUPD . */

/* \BeginLib */

/* \Routines called: */
/*     znaupd   ARPACK reverse communication interface routine. */
/*     zneupd   ARPACK routine that returns Ritz values and (optionally) */
/*             Ritz vectors. */
/*     zgttrf   LAPACK tridiagonal factorization routine. */
/*     zgttrs   LAPACK tridiagonal solve routine. */
/*     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully. */
/*     zaxpy    Level 1 BLAS that computes y <- alpha*x+y. */
/*     zcopy    Level 1 BLAS that copies one vector to another. */
/*     dznrm2   Level 1 BLAS that computes the norm of a complex vector. */
/*     av      Matrix vector multiplication routine that computes A*x. */
/*     mv      Matrix vector multiplication routine that computes M*x. */

/* \Author */
/*     Danny Sorensen */
/*     Richard Lehoucq */
/*     Chao Yang */
/*     Dept. of Computational & */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \SCCS Information: @(#) */
/* FILE: ndrv4.F   SID: 2.4   DATE OF SID: 10/18/00   RELEASE: 2 */

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
/*     | The user can modify NEV, NCV, SIGMA to solve       | */
/*     | problems of different sizes, and to get different  | */
/*     | parts of the spectrum.  However, The following     | */
/*     | conditions must be satisfied:                      | */
/*     |                     N <= MAXN,                     | */
/*     |                   NEV <= MAXNEV,                   | */
/*     |               NEV + 2 <= NCV <= MAXNCV             | */
/*     %----------------------------------------------------% */

    n = 100;
    nev = 4;
    ncv = 20;
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
    sigma.r = 1., sigma.i = 0.;

/*     %--------------------------------------------------% */
/*     | Construct C = A - SIGMA*M in COMPLEX arithmetic. | */
/*     | Factor C in COMPLEX arithmetic (using LAPACK     | */
/*     | subroutine zgttrf ). The matrix A is chosen to be | */
/*     | the tridiagonal matrix derived from the standard | */
/*     | central difference discretization of the 1-d     | */
/*     | convection-diffusion operator u``+ rho*u` on the | */
/*     | interval [0, 1] with zero Dirichlet boundary     | */
/*     | condition.  The matrix M is chosen to be the     | */
/*     | symmetric tridiagonal matrix with 4.0 on the     | */
/*     | diagonal and 1.0 on the off-diagonals.           | */
/*     %--------------------------------------------------% */

    convct_1.rho.r = 10., convct_1.rho.i = 0.;
    i__1 = n + 1;
    z__2.r = (doublereal) i__1, z__2.i = 0.;
    z_div(&z__1, &c_b1, &z__2);
    h__.r = z__1.r, h__.i = z__1.i;
    z_div(&z__1, &convct_1.rho, &c_b3);
    s.r = z__1.r, s.i = z__1.i;

    z__4.r = -1., z__4.i = -0.;
    z_div(&z__3, &z__4, &h__);
    z__2.r = z__3.r - s.r, z__2.i = z__3.i - s.i;
    z__6.r = sigma.r * h__.r - sigma.i * h__.i, z__6.i = sigma.r * h__.i + 
	    sigma.i * h__.r;
    z_div(&z__5, &z__6, &c_b5);
    z__1.r = z__2.r - z__5.r, z__1.i = z__2.i - z__5.i;
    s1.r = z__1.r, s1.i = z__1.i;
    z_div(&z__2, &c_b3, &h__);
    z__5.r = sigma.r * 4. - sigma.i * 0., z__5.i = sigma.i * 4. + sigma.r * 
	    0.;
    z__4.r = z__5.r * h__.r - z__5.i * h__.i, z__4.i = z__5.r * h__.i + 
	    z__5.i * h__.r;
    z_div(&z__3, &z__4, &c_b5);
    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
    s2.r = z__1.r, s2.i = z__1.i;
    z__4.r = -1., z__4.i = -0.;
    z_div(&z__3, &z__4, &h__);
    z__2.r = z__3.r + s.r, z__2.i = z__3.i + s.i;
    z__6.r = sigma.r * h__.r - sigma.i * h__.i, z__6.i = sigma.r * h__.i + 
	    sigma.i * h__.r;
    z_div(&z__5, &z__6, &c_b5);
    z__1.r = z__2.r - z__5.r, z__1.i = z__2.i - z__5.i;
    s3.r = z__1.r, s3.i = z__1.i;

    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j - 1;
	dl[i__2].r = s1.r, dl[i__2].i = s1.i;
	i__2 = j - 1;
	dd[i__2].r = s2.r, dd[i__2].i = s2.i;
	i__2 = j - 1;
	du[i__2].r = s3.r, du[i__2].i = s3.i;
/* L10: */
    }
    i__1 = n - 1;
    dd[i__1].r = s2.r, dd[i__1].i = s2.i;

    zgttrf_(&n, dl, dd, du, du2, ipiv, &ierr);
    if (ierr != 0) {
	s_wsle(&io___22);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___23);
	do_lio(&c__9, &c__1, " ERROR with _gttrf in _NDRV4.", (ftnlen)29);
	e_wsle();
	s_wsle(&io___24);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	goto L9000;
    }

/*     %-----------------------------------------------------% */
/*     | The work array WORKL is used in ZNAUPD  as           | */
/*     | workspace.  Its dimension LWORKL is set as          | */
/*     | illustrated below.  The parameter TOL determines    | */
/*     | the stopping criterion. If TOL<=0, machine          | */
/*     | precision is used.  The variable IDO is used for    | */
/*     | reverse communication, and is initially set to 0.   | */
/*     | Setting INFO=0 indicates that a random vector is    | */
/*     | generated in ZNAUPD  to start the Arnoldi iteration. | */
/*     %-----------------------------------------------------% */

/* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 * 3 + ncv * 5;
    tol = 0.f;
    ido = 0;
    info = 0;

/*     %---------------------------------------------------% */
/*     | This program uses exact shifts with respect to    | */
/*     | the current Hessenberg matrix (IPARAM(1) = 1).    | */
/*     | IPARAM(3) specifies the maximum number of Arnoldi | */
/*     | iterations allowed. Mode 3 of ZNAUPD  is used      | */
/*     | (IPARAM(7) = 3).  All these options can be        | */
/*     | changed by the user. For details see the          | */
/*     | documentation in ZNAUPD .                          | */
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
/*        | Repeatedly call the routine ZNAUPD  and take | */
/*        | actions indicated by parameter IDO until    | */
/*        | either convergence is indicated or maxitr   | */
/*        | has been exceeded.                          | */
/*        %---------------------------------------------% */

    znaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, 
	    iparam, ipntr, workd, workl, &lworkl, rwork, &info, (ftnlen)1, (
	    ftnlen)2);

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
	zgttrs_("N", &n, &c__1, dl, dd, du, du2, ipiv, &workd[ipntr[1] - 1], &
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
/*           | L O O P   B A C K to call ZNAUPD  again. | */
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

	zcopy_(&n, &workd[ipntr[2] - 1], &c__1, &workd[ipntr[1] - 1], &c__1);
	zgttrs_("N", &n, &c__1, dl, dd, du, du2, ipiv, &workd[ipntr[1] - 1], &
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
/*           | L O O P   B A C K to call ZNAUPD  again. | */
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
/*           | L O O P   B A C K to call ZNAUPD  again. | */
/*           %-----------------------------------------% */

	goto L20;

    }

/*     %-----------------------------------------% */
/*     | Either we have convergence, or there is | */
/*     | an error.                               | */
/*     %-----------------------------------------% */

    if (info < 0) {

/*        %----------------------------% */
/*        |  Error message, check the  | */
/*        |  documentation in ZNAUPD    | */
/*        %----------------------------% */

	s_wsle(&io___45);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___46);
	do_lio(&c__9, &c__1, " Error with _naupd, info = ", (ftnlen)27);
	do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___47);
	do_lio(&c__9, &c__1, " Check the documentation of _naupd.", (ftnlen)
		35);
	e_wsle();
	s_wsle(&io___48);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();

    } else {

/*        %-------------------------------------------% */
/*        | No fatal errors occurred.                 | */
/*        | Post-Process using ZNEUPD .                | */
/*        |                                           | */
/*        | Computed eigenvalues may be extracted.    | */
/*        |                                           | */
/*        | Eigenvectors may also be computed now if  | */
/*        | desired.  (indicated by rvec = .true.)    | */
/*        %-------------------------------------------% */

	rvec = TRUE_;

	zneupd_(&rvec, "A", select, d__, v, &c__256, &sigma, workev, bmat, &n,
		 which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, 
		workd, workl, &lworkl, rwork, &ierr, (ftnlen)1, (ftnlen)1, (
		ftnlen)2);

/*        %----------------------------------------------% */
/*        | Eigenvalues are returned in the one          | */
/*        | dimensional array D.  The corresponding      | */
/*        | eigenvectors are returned in the first NCONV | */
/*        | (=IPARAM(5)) columns of the two dimensional  | */
/*        | array V if requested.  Otherwise, an         | */
/*        | orthogonal basis for the invariant subspace  | */
/*        | corresponding to the eigenvalues in D is     | */
/*        | returned in V.                               | */
/*        %----------------------------------------------% */

	if (ierr != 0) {

/*           %------------------------------------% */
/*           | Error condition:                   | */
/*           | Check the documentation of ZNEUPD . | */
/*           %------------------------------------% */

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

	    nconv = iparam[4];
	    i__1 = nconv;
	    for (j = 1; j <= i__1; ++j) {

		av_(&n, &v[(j << 8) - 256], ax);
		mv_(&n, &v[(j << 8) - 256], mx);
		i__2 = j - 1;
		z__1.r = -d__[i__2].r, z__1.i = -d__[i__2].i;
		zaxpy_(&n, &z__1, mx, &c__1, ax, &c__1);
		i__2 = j - 1;
		rd[j - 1] = d__[i__2].r;
		rd[j + 24] = d_imag(&d__[j - 1]);
		rd[j + 49] = dznrm2_(&n, ax, &c__1);
		rd[j + 49] /= dlapy2_(&rd[j - 1], &rd[j + 24]);
/* L80: */
	    }

/*            %-----------------------------% */
/*            | Display computed residuals. | */
/*            %-----------------------------% */

	    dmout_(&c__6, &nconv, &c__3, rd, &c__25, &c_n6, "Ritz values (Re"
		    "al, Imag) and direct residuals", (ftnlen)45);

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
	do_lio(&c__9, &c__1, "_NDRV4 ", (ftnlen)7);
	e_wsle();
	s_wsle(&io___69);
	do_lio(&c__9, &c__1, "====== ", (ftnlen)7);
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
	do_lio(&c__5, &c__1, (char *)&tol, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___79);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();

    }

L9000:

    return 0;
} /* MAIN__ */


/* ========================================================================== */

/*     matrix vector multiplication subroutine */

/* Subroutine */ int mv_(integer *n, doublecomplex *v, doublecomplex *w)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    doublecomplex h__;
    integer j;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);


/*     Compute the matrix vector multiplication y<---M*x */
/*     where M is a n by n symmetric tridiagonal matrix with 4 on the */
/*     diagonal, 1 on the subdiagonal and superdiagonal. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    z__3.r = v[1].r * 4. - v[1].i * 0., z__3.i = v[1].i * 4. + v[1].r * 0.;
    z__4.r = v[2].r * 1. - v[2].i * 0., z__4.i = v[2].i * 1. + v[2].r * 0.;
    z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
    z_div(&z__1, &z__2, &c_b5);
    w[1].r = z__1.r, w[1].i = z__1.i;
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j - 1;
	z__4.r = v[i__3].r * 1. - v[i__3].i * 0., z__4.i = v[i__3].i * 1. + v[
		i__3].r * 0.;
	i__4 = j;
	z__5.r = v[i__4].r * 4. - v[i__4].i * 0., z__5.i = v[i__4].i * 4. + v[
		i__4].r * 0.;
	z__3.r = z__4.r + z__5.r, z__3.i = z__4.i + z__5.i;
	i__5 = j + 1;
	z__6.r = v[i__5].r * 1. - v[i__5].i * 0., z__6.i = v[i__5].i * 1. + v[
		i__5].r * 0.;
	z__2.r = z__3.r + z__6.r, z__2.i = z__3.i + z__6.i;
	z_div(&z__1, &z__2, &c_b5);
	w[i__2].r = z__1.r, w[i__2].i = z__1.i;
/* L40: */
    }
    i__1 = *n;
    i__2 = *n - 1;
    z__3.r = v[i__2].r * 1. - v[i__2].i * 0., z__3.i = v[i__2].i * 1. + v[
	    i__2].r * 0.;
    i__3 = *n;
    z__4.r = v[i__3].r * 4. - v[i__3].i * 0., z__4.i = v[i__3].i * 4. + v[
	    i__3].r * 0.;
    z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
    z_div(&z__1, &z__2, &c_b5);
    w[i__1].r = z__1.r, w[i__1].i = z__1.i;

    i__1 = *n + 1;
    z__2.r = (doublereal) i__1, z__2.i = 0.;
    z_div(&z__1, &c_b1, &z__2);
    h__.r = z__1.r, h__.i = z__1.i;
    zscal_(n, &h__, &w[1], &c__1);
    return 0;
} /* mv_ */

/* ------------------------------------------------------------------ */
/* Subroutine */ int av_(integer *n, doublecomplex *v, doublecomplex *w)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    doublecomplex h__;
    integer j;
    doublecomplex s, dd, dl, du;


    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    i__1 = *n + 1;
    z__2.r = (doublereal) i__1, z__2.i = 0.;
    z_div(&z__1, &c_b1, &z__2);
    h__.r = z__1.r, h__.i = z__1.i;
    z_div(&z__1, &convct_1.rho, &c_b3);
    s.r = z__1.r, s.i = z__1.i;
    z_div(&z__1, &c_b3, &h__);
    dd.r = z__1.r, dd.i = z__1.i;
    z__3.r = -1., z__3.i = -0.;
    z_div(&z__2, &z__3, &h__);
    z__1.r = z__2.r - s.r, z__1.i = z__2.i - s.i;
    dl.r = z__1.r, dl.i = z__1.i;
    z__3.r = -1., z__3.i = -0.;
    z_div(&z__2, &z__3, &h__);
    z__1.r = z__2.r + s.r, z__1.i = z__2.i + s.i;
    du.r = z__1.r, du.i = z__1.i;

    z__2.r = dd.r * v[1].r - dd.i * v[1].i, z__2.i = dd.r * v[1].i + dd.i * v[
	    1].r;
    z__3.r = du.r * v[2].r - du.i * v[2].i, z__3.i = du.r * v[2].i + du.i * v[
	    2].r;
    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
    w[1].r = z__1.r, w[1].i = z__1.i;
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j - 1;
	z__3.r = dl.r * v[i__3].r - dl.i * v[i__3].i, z__3.i = dl.r * v[i__3]
		.i + dl.i * v[i__3].r;
	i__4 = j;
	z__4.r = dd.r * v[i__4].r - dd.i * v[i__4].i, z__4.i = dd.r * v[i__4]
		.i + dd.i * v[i__4].r;
	z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
	i__5 = j + 1;
	z__5.r = du.r * v[i__5].r - du.i * v[i__5].i, z__5.i = du.r * v[i__5]
		.i + du.i * v[i__5].r;
	z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
	w[i__2].r = z__1.r, w[i__2].i = z__1.i;
/* L40: */
    }
    i__1 = *n;
    i__2 = *n - 1;
    z__2.r = dl.r * v[i__2].r - dl.i * v[i__2].i, z__2.i = dl.r * v[i__2].i + 
	    dl.i * v[i__2].r;
    i__3 = *n;
    z__3.r = dd.r * v[i__3].r - dd.i * v[i__3].i, z__3.i = dd.r * v[i__3].i + 
	    dd.i * v[i__3].r;
    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
    w[i__1].r = z__1.r, w[i__1].i = z__1.i;
    return 0;
} /* av_ */

/* Main program alias */ int zndrv4_ () { MAIN__ (); return 0; }