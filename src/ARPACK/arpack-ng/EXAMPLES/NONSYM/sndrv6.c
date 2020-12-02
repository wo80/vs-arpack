/* EXAMPLES\NONSYM\sndrv6.f -- translated by f2c (version 20100827).
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
    integer i__1, i__2, i__3;
    real r__1, r__2;
    complex q__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double r_imag(complex *);

    /* Local variables */
    real d__[75]	/* was [25][3] */;
    integer j, n;
    real v[6400]	/* was [256][25] */;
    complex c1, c2, c3;
    extern /* Subroutine */ int av_(integer *, real *, real *);
    real ax[256];
    extern /* Subroutine */ int mv_(integer *, real *, real *);
    real mx[256];
    complex cdd[256], cdl[256], cdu[256];
    integer ido, ncv, nev;
    real tol;
    complex cdu2[256];
    real deni;
    char bmat[1];
    integer mode;
    real denr;
    integer info;
    logical rvec;
    integer ierr, ipiv[256];
    real numi;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    real numr;
    extern doublereal snrm2_(integer *, real *, integer *);
    char which[2];
    real resid[256];
    complex ctemp[256];
    integer nconv;
    real workd[768];
    logical first;
    integer ipntr[14];
    real workl[2025];
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *), smout_(integer *, integer *, integer *, real *
	    , integer *, integer *, char *, ftnlen);
    extern doublereal slapy2_(real *, real *);
    integer iparam[11];
    real sigmai;
    logical select[25];
    real sigmar;
    extern /* Subroutine */ int cgttrf_(integer *, complex *, complex *, 
	    complex *, complex *, integer *, integer *), snaupd_(integer *, 
	    char *, integer *, char *, integer *, real *, real *, integer *, 
	    real *, integer *, integer *, integer *, real *, real *, integer *
	    , integer *, ftnlen, ftnlen), sneupd_(logical *, char *, logical *
	    , real *, real *, real *, integer *, real *, real *, real *, char 
	    *, integer *, char *, integer *, real *, real *, integer *, real *
	    , integer *, integer *, integer *, real *, real *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    integer ishfts, maxitr;
    extern /* Subroutine */ int cgttrs_(char *, integer *, integer *, complex 
	    *, complex *, complex *, complex *, integer *, complex *, integer 
	    *, integer *, ftnlen);
    integer lworkl;
    real workev[75];

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 6, 0, 0, 0 };
    static cilist io___5 = { 0, 6, 0, 0, 0 };
    static cilist io___6 = { 0, 6, 0, 0, 0 };
    static cilist io___21 = { 0, 6, 0, 0, 0 };
    static cilist io___22 = { 0, 6, 0, 0, 0 };
    static cilist io___23 = { 0, 6, 0, 0, 0 };
    static cilist io___38 = { 0, 6, 0, 0, 0 };
    static cilist io___39 = { 0, 6, 0, 0, 0 };
    static cilist io___40 = { 0, 6, 0, 0, 0 };
    static cilist io___41 = { 0, 6, 0, 0, 0 };
    static cilist io___42 = { 0, 6, 0, 0, 0 };
    static cilist io___43 = { 0, 6, 0, 0, 0 };
    static cilist io___44 = { 0, 6, 0, 0, 0 };
    static cilist io___45 = { 0, 6, 0, 0, 0 };
    static cilist io___46 = { 0, 6, 0, 0, 0 };
    static cilist io___47 = { 0, 6, 0, 0, 0 };
    static cilist io___52 = { 0, 6, 0, 0, 0 };
    static cilist io___53 = { 0, 6, 0, 0, 0 };
    static cilist io___54 = { 0, 6, 0, 0, 0 };
    static cilist io___55 = { 0, 6, 0, 0, 0 };
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
    static cilist io___80 = { 0, 6, 0, 0, 0 };
    static cilist io___81 = { 0, 6, 0, 0, 0 };
    static cilist io___82 = { 0, 6, 0, 0, 0 };



/*     Simple program to illustrate the idea of reverse communication */
/*     in shift-invert mode for a generalized nonsymmetric eigenvalue problem. */

/*     We implement example six of ex-nonsym.doc in DOCUMENTS directory */

/* \Example-6 */

/*     ... Suppose we want to solve A*x = lambda*B*x in shift-invert mode */
/*         The matrix A is the tridiagonal matrix with 2 on the diagonal, */
/*         -2 on the subdiagonal and 3 on the superdiagonal.  The matrix M */
/*         is the tridiagonal matrix with 4 on the diagonal and 1 on the */
/*         off-diagonals. */
/*     ... The shift sigma is a complex number (sigmar, sigmai). */
/*     ... OP = Imaginary_Part{inv[A-(SIGMAR,SIGMAI)*M]*M  and  B = M. */
/*     ... Use mode 4 of SNAUPD. */

/* \BeginLib */

/* \Routines called: */
/*     snaupd  ARPACK reverse communication interface routine. */
/*     sneupd  ARPACK routine that returns Ritz values and (optionally) */
/*             Ritz vectors. */
/*     cgttrf  LAPACK complex matrix factorization routine. */
/*     cgttrs  LAPACK complex linear system solve routine. */
/*     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully. */
/*     saxpy   Level 1 BLAS that computes y <- alpha*x+y. */
/*     sdot    Level 1 BLAS that computes the dot product of two vectors. */
/*     snrm2   Level 1 BLAS that computes the norm of a vector. */
/*     av      Matrix vector subroutine that computes A*x. */
/*     mv      Matrix vector subroutine that computes M*x. */

/* \Author */
/*     Richard Lehoucq */
/*     Danny Sorensen */
/*     Chao Yang */
/*     Dept. of Computational & */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \SCCS Information: @(#) */
/* FILE: ndrv6.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2 */

/* \Remarks */
/*     1. None */

/* \EndLib */
/* ------------------------------------------------------------------------- */

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


/*     %--------------------% */
/*     | Intrinsic Function | */
/*     %--------------------% */


/*     %-----------------------% */
/*     | Executable statements | */
/*     %-----------------------% */

/*     %----------------------------------------------------% */
/*     | The number N is the dimension of the matrix.  A    | */
/*     | generalized eigenvalue problem is solved (BMAT =   | */
/*     | 'G').  NEV is the number of eigenvalues (closest   | */
/*     | to the shift (SIGMAR,SIGMAI)) to be approximated.  | */
/*     | Since the shift-invert mode is used, WHICH is set  | */
/*     | to 'LM'.  The user can modify NEV, NCV, SIGMA to   | */
/*     | solve problems of different sizes, and to get      | */
/*     | different parts of the spectrum. However, The      | */
/*     | following conditions must be satisfied:            | */
/*     |                     N <= MAXN,                     | */
/*     |                   NEV <= MAXNEV,                   | */
/*     |               NEV + 2 <= NCV <= MAXNCV             | */
/*     %----------------------------------------------------% */

    n = 100;
    nev = 4;
    ncv = 20;
    if (n > 256) {
	s_wsle(&io___4);
	do_lio(&c__9, &c__1, " ERROR with _NDRV6: N is greater than MAXN ", (
		ftnlen)43);
	e_wsle();
	goto L9000;
    } else if (nev > 10) {
	s_wsle(&io___5);
	do_lio(&c__9, &c__1, " ERROR with _NDRV6: NEV is greater than MAXNEV "
		, (ftnlen)47);
	e_wsle();
	goto L9000;
    } else if (ncv > 25) {
	s_wsle(&io___6);
	do_lio(&c__9, &c__1, " ERROR with _NDRV6: NCV is greater than MAXNCV "
		, (ftnlen)47);
	e_wsle();
	goto L9000;
    }
    *(unsigned char *)bmat = 'G';
    s_copy(which, "LM", (ftnlen)2, (ftnlen)2);
    sigmar = .4f;
    sigmai = .6f;

/*     %----------------------------------------------------% */
/*     | Construct C = A - (SIGMAR,SIGMAI)*M in complex     | */
/*     | arithmetic, and factor C in complex arithmetic     | */
/*     | (using LAPACK subroutine cgttrf). The matrix A is  | */
/*     | chosen to be the tridiagonal matrix with -2 on the | */
/*     | subdiagonal, 2 on the diagonal and 3 on the        | */
/*     | superdiagonal. The matrix M is chosen to be the    | */
/*     | symmetric tridiagonal matrix with 4 on the         | */
/*     | diagonal and 1 on the off-diagonals.               | */
/*     %----------------------------------------------------% */

    r__1 = -2.f - sigmar;
    r__2 = -sigmai;
    q__1.r = r__1, q__1.i = r__2;
    c1.r = q__1.r, c1.i = q__1.i;
    r__1 = 2.f - sigmar * 4.f;
    r__2 = sigmai * -4.f;
    q__1.r = r__1, q__1.i = r__2;
    c2.r = q__1.r, c2.i = q__1.i;
    r__1 = 3.f - sigmar;
    r__2 = -sigmai;
    q__1.r = r__1, q__1.i = r__2;
    c3.r = q__1.r, c3.i = q__1.i;

    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j - 1;
	cdl[i__2].r = c1.r, cdl[i__2].i = c1.i;
	i__2 = j - 1;
	cdd[i__2].r = c2.r, cdd[i__2].i = c2.i;
	i__2 = j - 1;
	cdu[i__2].r = c3.r, cdu[i__2].i = c3.i;
/* L10: */
    }
    i__1 = n - 1;
    cdd[i__1].r = c2.r, cdd[i__1].i = c2.i;

    cgttrf_(&n, cdl, cdd, cdu, cdu2, ipiv, &ierr);
    if (ierr != 0) {
	s_wsle(&io___21);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___22);
	do_lio(&c__9, &c__1, " ERROR with _gttrf in _NDRV6.", (ftnlen)29);
	e_wsle();
	s_wsle(&io___23);
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
/*     | This program uses exact shift with respect to     | */
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

/*           %------------------------------------------------------------% */
/*           |                           Perform                          | */
/*           | y <--- OP*x = Imaginary_Part{inv[A-(SIGMAR,SIGMAI)*M]*M*x} | */
/*           | to force starting vector into the range of OP. The user    | */
/*           | should supply his/her own matrix vector multiplication     | */
/*           | routine and a complex linear system solver.  The matrix    | */
/*           | vector multiplication routine should take workd(ipntr(1))  | */
/*           | as the input. The final result (a real vector) should be   | */
/*           | returned to workd(ipntr(2)).                               | */
/*           %------------------------------------------------------------% */

	mv_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = j - 1;
	    i__3 = ipntr[1] + j - 2;
	    q__1.r = workd[i__3], q__1.i = 0.f;
	    ctemp[i__2].r = q__1.r, ctemp[i__2].i = q__1.i;
/* L30: */
	}

	cgttrs_("N", &n, &c__1, cdl, cdd, cdu, cdu2, ipiv, ctemp, &c__256, &
		ierr, (ftnlen)1);
	if (ierr != 0) {
	    s_wsle(&io___38);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___39);
	    do_lio(&c__9, &c__1, " ERROR with _gttrs in _NDRV6.", (ftnlen)29);
	    e_wsle();
	    s_wsle(&io___40);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    goto L9000;
	}
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
	    workd[ipntr[1] + j - 2] = r_imag(&ctemp[j - 1]);
/* L40: */
	}

/*           %-----------------------------------------% */
/*           | L O O P   B A C K to call SNAUPD again. | */
/*           %-----------------------------------------% */

	goto L20;

    } else if (ido == 1) {

/*           %------------------------------------------------------------% */
/*           |                          Perform                           | */
/*           | y <--- OP*x = Imaginary_Part{inv[A-(SIGMAR,SIGMAI)*M]*M*x} | */
/*           | M*x has been saved in workd(ipntr(3)). The user only need  | */
/*           | the complex linear system solver here that takes           | */
/*           | complex[workd(ipntr(3))] as input, and returns the result  | */
/*           | to workd(ipntr(2)).                                        | */
/*           %------------------------------------------------------------% */

	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = j - 1;
	    i__3 = ipntr[2] + j - 2;
	    q__1.r = workd[i__3], q__1.i = 0.f;
	    ctemp[i__2].r = q__1.r, ctemp[i__2].i = q__1.i;
/* L50: */
	}
	cgttrs_("N", &n, &c__1, cdl, cdd, cdu, cdu2, ipiv, ctemp, &c__256, &
		ierr, (ftnlen)1);
	if (ierr != 0) {
	    s_wsle(&io___41);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___42);
	    do_lio(&c__9, &c__1, " ERROR with _gttrs in _NDRV6.", (ftnlen)29);
	    e_wsle();
	    s_wsle(&io___43);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    goto L9000;
	}
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
	    workd[ipntr[1] + j - 2] = r_imag(&ctemp[j - 1]);
/* L60: */
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

	s_wsle(&io___44);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___45);
	do_lio(&c__9, &c__1, " Error with _naupd, info = ", (ftnlen)27);
	do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___46);
	do_lio(&c__9, &c__1, " Check the documentation of _naupd.", (ftnlen)
		35);
	e_wsle();
	s_wsle(&io___47);
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

/*           %------------------------------------% */
/*           | Error condition:                   | */
/*           | Check the documentation of SNEUPD. | */
/*           %------------------------------------% */

	    s_wsle(&io___52);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___53);
	    do_lio(&c__9, &c__1, " Error with _neupd, info = ", (ftnlen)27);
	    do_lio(&c__3, &c__1, (char *)&ierr, (ftnlen)sizeof(integer));
	    e_wsle();
	    s_wsle(&io___54);
	    do_lio(&c__9, &c__1, " Check the documentation of _neupd. ", (
		    ftnlen)36);
	    e_wsle();
	    s_wsle(&io___55);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();

	} else {

	    first = TRUE_;
	    nconv = iparam[4];
	    i__1 = nconv;
	    for (j = 1; j <= i__1; ++j) {

/*               %-------------------------------------% */
/*               | Use Rayleigh Quotient to recover    | */
/*               | eigenvalues of the original problem.| */
/*               %-------------------------------------% */

		if (d__[j + 24] == 0.f) {

/*                  %----------------------------% */
/*                  | Eigenvalue is real.        | */
/*                  | Compute d = x'(Ax)/x'(Mx). | */
/*                  %----------------------------% */

		    av_(&n, &v[(j << 8) - 256], ax);
		    numr = sdot_(&n, &v[(j << 8) - 256], &c__1, ax, &c__1);
		    mv_(&n, &v[(j << 8) - 256], ax);
		    denr = sdot_(&n, &v[(j << 8) - 256], &c__1, ax, &c__1);
		    d__[j - 1] = numr / denr;

		} else if (first) {

/*                  %------------------------% */
/*                  | Eigenvalue is complex. | */
/*                  | Compute the first one  | */
/*                  | of the conjugate pair. | */
/*                  %------------------------% */

/*                  %----------------% */
/*                  | Compute x'(Ax) | */
/*                  %----------------% */

		    av_(&n, &v[(j << 8) - 256], ax);
		    numr = sdot_(&n, &v[(j << 8) - 256], &c__1, ax, &c__1);
		    numi = sdot_(&n, &v[(j + 1 << 8) - 256], &c__1, ax, &c__1)
			    ;
		    av_(&n, &v[(j + 1 << 8) - 256], ax);
		    numr += sdot_(&n, &v[(j + 1 << 8) - 256], &c__1, ax, &
			    c__1);
		    numi = -numi + sdot_(&n, &v[(j << 8) - 256], &c__1, ax, &
			    c__1);

/*                  %----------------% */
/*                  | Compute x'(Mx) | */
/*                  %----------------% */

		    mv_(&n, &v[(j << 8) - 256], ax);
		    denr = sdot_(&n, &v[(j << 8) - 256], &c__1, ax, &c__1);
		    deni = sdot_(&n, &v[(j + 1 << 8) - 256], &c__1, ax, &c__1)
			    ;
		    mv_(&n, &v[(j + 1 << 8) - 256], ax);
		    denr += sdot_(&n, &v[(j + 1 << 8) - 256], &c__1, ax, &
			    c__1);
		    deni = -deni + sdot_(&n, &v[(j << 8) - 256], &c__1, ax, &
			    c__1);

/*                  %----------------% */
/*                  | d=x'(Ax)/x'(Mx)| */
/*                  %----------------% */

		    d__[j - 1] = (numr * denr + numi * deni) / slapy2_(&denr, 
			    &deni);
		    d__[j + 24] = (numi * denr - numr * deni) / slapy2_(&denr,
			     &deni);
		    first = FALSE_;

		} else {

/*                  %------------------------------% */
/*                  | Get the second eigenvalue of | */
/*                  | the conjugate pair by taking | */
/*                  | the conjugate of the last    | */
/*                  | eigenvalue computed.         | */
/*                  %------------------------------% */

		    d__[j - 1] = d__[j - 2];
		    d__[j + 24] = -d__[j + 23];
		    first = TRUE_;

		}

/* L70: */
	    }

/*           %---------------------------% */
/*           | Compute the residual norm | */
/*           |                           | */
/*           |   ||  A*x - lambda*x ||   | */
/*           |                           | */
/*           | for the NCONV accurately  | */
/*           | computed eigenvalues and  | */
/*           | eigenvectors.  (iparam(5) | */
/*           | indicates how many are    | */
/*           | accurate to the requested | */
/*           | tolerance)                | */
/*           %---------------------------% */

	    first = TRUE_;
	    nconv = iparam[4];
	    i__1 = nconv;
	    for (j = 1; j <= i__1; ++j) {

		if (d__[j + 24] == 0.f) {

/*                 %--------------------% */
/*                 | Ritz value is real | */
/*                 %--------------------% */

		    av_(&n, &v[(j << 8) - 256], ax);
		    mv_(&n, &v[(j << 8) - 256], mx);
		    r__1 = -d__[j - 1];
		    saxpy_(&n, &r__1, mx, &c__1, ax, &c__1);
		    d__[j + 49] = snrm2_(&n, ax, &c__1);
		    d__[j + 49] /= (r__1 = d__[j - 1], dabs(r__1));

		} else if (first) {

/*                 %------------------------% */
/*                 | Ritz value is complex  | */
/*                 | Residual of one Ritz   | */
/*                 | value of the conjugate | */
/*                 | pair is computed.      | */
/*                 %------------------------% */

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

/* L80: */
	    }

/*           %-----------------------------% */
/*           | Display computed residuals. | */
/*           %-----------------------------% */

	    smout_(&c__6, &nconv, &c__3, d__, &c__25, &c_n6, "Ritz values (R"
		    "eal,Imag) and relative residuals", (ftnlen)46);

	}

/*       %-------------------------------------------% */
/*       | Print additional convergence information. | */
/*       %-------------------------------------------% */

	if (info == 1) {
	    s_wsle(&io___64);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___65);
	    do_lio(&c__9, &c__1, " Maximum number of iterations reached.", (
		    ftnlen)38);
	    e_wsle();
	    s_wsle(&io___66);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	} else if (info == 3) {
	    s_wsle(&io___67);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___68);
	    do_lio(&c__9, &c__1, " No shifts could be applied during implicit"
		    , (ftnlen)43);
	    do_lio(&c__9, &c__1, " Arnoldi update, try increasing NCV.", (
		    ftnlen)36);
	    e_wsle();
	    s_wsle(&io___69);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	}

	s_wsle(&io___70);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___71);
	do_lio(&c__9, &c__1, " _NDRV6 ", (ftnlen)8);
	e_wsle();
	s_wsle(&io___72);
	do_lio(&c__9, &c__1, " ====== ", (ftnlen)8);
	e_wsle();
	s_wsle(&io___73);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___74);
	do_lio(&c__9, &c__1, " Size of the matrix is ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___75);
	do_lio(&c__9, &c__1, " The number of Ritz values requested is ", (
		ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&nev, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___76);
	do_lio(&c__9, &c__1, " The number of Arnoldi vectors generated", (
		ftnlen)40);
	do_lio(&c__9, &c__1, " (NCV) is ", (ftnlen)10);
	do_lio(&c__3, &c__1, (char *)&ncv, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___77);
	do_lio(&c__9, &c__1, " What portion of the spectrum: ", (ftnlen)31);
	do_lio(&c__9, &c__1, which, (ftnlen)2);
	e_wsle();
	s_wsle(&io___78);
	do_lio(&c__9, &c__1, " The number of converged Ritz values is ", (
		ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&nconv, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___79);
	do_lio(&c__9, &c__1, " The number of Implicit Arnoldi update", (
		ftnlen)38);
	do_lio(&c__9, &c__1, " iterations taken is ", (ftnlen)21);
	do_lio(&c__3, &c__1, (char *)&iparam[2], (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___80);
	do_lio(&c__9, &c__1, " The number of OP*x is ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&iparam[8], (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___81);
	do_lio(&c__9, &c__1, " The convergence criterion is ", (ftnlen)30);
	do_lio(&c__4, &c__1, (char *)&tol, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___82);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();

    }

/*     %---------------------------% */
/*     | Done with program sndrv6. | */
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
    integer j;


/*     Compute the matrix vector multiplication y<---M*x */
/*     where M is a n by n symmetric tridiagonal matrix with 4 on the */
/*     diagonal, 1 on the subdiagonal and superdiagonal. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    w[1] = v[1] * 4.f + v[2] * 1.f;
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j) {
	w[j] = v[j - 1] * 1.f + v[j] * 4.f + v[j + 1] * 1.f;
/* L10: */
    }
    w[*n] = v[*n - 1] * 1.f + v[*n] * 4.f;
    return 0;
} /* mv_ */

/* ------------------------------------------------------------------ */
/* Subroutine */ int av_(integer *n, real *v, real *w)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer j;


/*     Compute the matrix vector multiplication y<---A*x */
/*     where M is a n by n symmetric tridiagonal matrix with 2 on the */
/*     diagonal, -2 on the subdiagonal and 3 on the superdiagonal. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    w[1] = v[1] * 2.f + v[2] * 3.f;
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j) {
	w[j] = v[j - 1] * -2.f + v[j] * 2.f + v[j + 1] * 3.f;
/* L10: */
    }
    w[*n] = v[*n - 1] * -2.f + v[*n] * 2.f;
    return 0;
} /* av_ */

/* Main program alias */ int sndrv6_ () { MAIN__ (); return 0; }
