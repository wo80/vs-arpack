/* EXAMPLES\SYM\ssdrv2.f -- translated by f2c (version 20100827). */

#include "arpack.h"

int ssdrv2()
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    real d__[50]	/* was [25][2] */;
    integer j, n;
    real v[6400]	/* was [256][25] */, h2, ad[256];
    real ax[256], adl[256], adu[256];
    integer ido, ncv, nev;
    real tol, adu2[256];
    char bmat[1];
    integer mode, info;
    logical rvec;
    integer ierr, ipiv[256];
    real sigma;
    char which[2];
    real resid[256];
    integer nconv;
    real workd[768];
    integer ipntr[11];
    real workl[825];
    integer iparam[11];
    logical select[25];
    integer ishfts;
    integer maxitr;
    integer lworkl;

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 6, 0, 0, 0 };
    static cilist io___5 = { 0, 6, 0, 0, 0 };
    static cilist io___6 = { 0, 6, 0, 0, 0 };
    static cilist io___26 = { 0, 6, 0, 0, 0 };
    static cilist io___27 = { 0, 6, 0, 0, 0 };
    static cilist io___28 = { 0, 6, 0, 0, 0 };
    static cilist io___34 = { 0, 6, 0, 0, 0 };
    static cilist io___35 = { 0, 6, 0, 0, 0 };
    static cilist io___36 = { 0, 6, 0, 0, 0 };
    static cilist io___37 = { 0, 6, 0, 0, 0 };
    static cilist io___38 = { 0, 6, 0, 0, 0 };
    static cilist io___39 = { 0, 6, 0, 0, 0 };
    static cilist io___40 = { 0, 6, 0, 0, 0 };
    static cilist io___44 = { 0, 6, 0, 0, 0 };
    static cilist io___45 = { 0, 6, 0, 0, 0 };
    static cilist io___46 = { 0, 6, 0, 0, 0 };
    static cilist io___47 = { 0, 6, 0, 0, 0 };
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
    static cilist io___61 = { 0, 6, 0, 0, 0 };
    static cilist io___62 = { 0, 6, 0, 0, 0 };
    static cilist io___63 = { 0, 6, 0, 0, 0 };
    static cilist io___64 = { 0, 6, 0, 0, 0 };
    static cilist io___65 = { 0, 6, 0, 0, 0 };
    static cilist io___66 = { 0, 6, 0, 0, 0 };
    static cilist io___67 = { 0, 6, 0, 0, 0 };
    static cilist io___68 = { 0, 6, 0, 0, 0 };

/*     Program to illustrate the idea of reverse communication */
/*     in shift and invert mode for a standard symmetric eigenvalue */
/*     problem.  The following program uses the two LAPACK subroutines */
/*     sgttrf.f and sgttrs.f to factor and solve a tridiagonal system of */
/*     equations. */

/*     We implement example two of ex-sym.doc in DOCUMENTS directory */

/* \Example-2 */
/*     ... Suppose we want to solve A*x = lambda*x in shift-invert mode, */
/*         where A is derived from the central difference discretization */
/*         of the 1-dimensional Laplacian on [0,1]  with zero Dirichlet */
/*         boundary condition. */
/*     ... OP = (inv[A - sigma*I]) and  B = I. */
/*     ... Use mode 3 of SSAUPD. */

/* \BeginLib */

/* \Routines called: */
/*     ssaupd  ARPACK reverse communication interface routine. */
/*     sseupd  ARPACK routine that returns Ritz values and (optionally) */
/*             Ritz vectors. */
/*     sgttrf  LAPACK tridiagonal factorization routine. */
/*     sgttrs  LAPACK tridiagonal solve routine. */
/*     saxpy   saxpy   Level 1 BLAS that computes y <- alpha*x+y. */
/*     snrm2   Level 1 BLAS that computes the norm of a vector. */
/*     av      Matrix vector multiplication routine that computes A*x. */

/* \Author */
/*     Richard Lehoucq */
/*     Danny Sorensen */
/*     Chao Yang */
/*     Dept. of Computational & */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \SCCS Information: @(#) */
/* FILE: sdrv2.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2 */

/* \Remarks */
/*     1. None */

/* \EndLib */
/* ---------------------------------------------------------------------- */

/*     %-----------------------------% */
/*     | Define leading dimensions   | */
/*     | for all arrays.             | */
/*     | MAXN:   Maximum dimension   | */
/*     |         of the A allowed.   | */
/*     | MAXNEV: Maximum NEV allowed | */
/*     | MAXNCV: Maximum NCV allowed | */
/*     %-----------------------------% */

/*     %-----------------------% */
/*     | Executable Statements | */
/*     %-----------------------% */

/*     %----------------------------------------------------% */
/*     | The number N is the dimension of the matrix.  A    | */
/*     | standard eigenvalue problem is solved (BMAT = 'I'. | */
/*     | NEV is the number of eigenvalues (closest to       | */
/*     | SIGMA) to be approximated.  Since the shift-invert | */
/*     | mode is used, WHICH is set to 'LM'.  The user can  | */
/*     | modify NEV, NCV, SIGMA to solve problems of        | */
/*     | different sizes, and to get different parts of the | */
/*     | spectrum.  However, The following conditions must  | */
/*     | be satisfied:                                      | */
/*     |                   N <= MAXN,                       | */
/*     |                 NEV <= MAXNEV,                     | */
/*     |             NEV + 1 <= NCV <= MAXNCV               | */
/*     %----------------------------------------------------% */

    n = 100;
    nev = 4;
    ncv = 10;
    if (n > 256) {
	s_wsle(&io___4);
	do_lio(&c__9, &c__1, " ERROR with _SDRV2: N is greater than MAXN ", (
		ftnlen)43);
	e_wsle();
	goto L9000;
    } else if (nev > 10) {
	s_wsle(&io___5);
	do_lio(&c__9, &c__1, " ERROR with _SDRV2: NEV is greater than MAXNEV "
		, (ftnlen)47);
	e_wsle();
	goto L9000;
    } else if (ncv > 25) {
	s_wsle(&io___6);
	do_lio(&c__9, &c__1, " ERROR with _SDRV2: NCV is greater than MAXNCV "
		, (ftnlen)47);
	e_wsle();
	goto L9000;
    }

    *(unsigned char *)bmat = 'I';
    s_copy(which, "LM", (ftnlen)2, (ftnlen)2);
    sigma = 0.f;

/*     %--------------------------------------------------% */
/*     | The work array WORKL is used in SSAUPD as        | */
/*     | workspace.  Its dimension LWORKL is set as       | */
/*     | illustrated below.  The parameter TOL determines | */
/*     | the stopping criterion.  If TOL<=0, machine      | */
/*     | precision is used.  The variable IDO is used for | */
/*     | reverse communication and is initially set to 0. | */
/*     | Setting INFO=0 indicates that a random vector is | */
/*     | generated in SSAUPD to start the Arnoldi         | */
/*     | iteration.                                       | */
/*     %--------------------------------------------------% */

    lworkl = ncv * (ncv + 8);
    tol = 0.f;
    ido = 0;
    info = 0;

/*     %---------------------------------------------------% */
/*     | This program uses exact shifts with respect to    | */
/*     | the current Hessenberg matrix (IPARAM(1) = 1).    | */
/*     | IPARAM(3) specifies the maximum number of Arnoldi | */
/*     | iterations allowed.  Mode 3 of SSAUPD is used     | */
/*     | (IPARAM(7) = 3).  All these options may be        | */
/*     | changed by the user. For details, see the         | */
/*     | documentation in SSAUPD.                          | */
/*     %---------------------------------------------------% */

    ishfts = 1;
    maxitr = 300;
    mode = 3;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

/*     %-----------------------------------------------------% */
/*     | Call LAPACK routine to factor (A-SIGMA*I), where A  | */
/*     | is the 1-d Laplacian.                               | */
/*     %-----------------------------------------------------% */

    h2 = 1.f / (real) ((n + 1) * (n + 1));
    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	ad[j - 1] = 2.f / h2 - sigma;
	adl[j - 1] = -1.f / h2;
/* L20: */
    }
    scopy_(&n, adl, &c__1, adu, &c__1);
    sgttrf_(&n, adl, ad, adu, adu2, ipiv, &ierr);
    if (ierr != 0) {
	s_wsle(&io___26);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___27);
	do_lio(&c__9, &c__1, " Error with _gttrf in SDRV2.", (ftnlen)28);
	e_wsle();
	s_wsle(&io___28);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	goto L9000;
    }

/*     %-------------------------------------------% */
/*     | M A I N   L O O P (Reverse communication) | */
/*     %-------------------------------------------% */

L10:

/*        %---------------------------------------------% */
/*        | Repeatedly call the routine SSAUPD and take | */
/*        | actions indicated by parameter IDO until    | */
/*        | either convergence is indicated or maxitr   | */
/*        | has been exceeded.                          | */
/*        %---------------------------------------------% */

    ssaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, 
	    iparam, ipntr, workd, workl, &lworkl, &info, (ftnlen)1, (ftnlen)2)
	    ;

    if (ido == -1 || ido == 1) {

/*           %----------------------------------------% */
/*           | Perform y <-- OP*x = inv[A-sigma*I]*x. | */
/*           | The user only need the linear system   | */
/*           | solver here that takes workd(ipntr(1)) | */
/*           | as input, and returns the result to    | */
/*           | workd(ipntr(2)).                       | */
/*           %----------------------------------------% */

	scopy_(&n, &workd[ipntr[0] - 1], &c__1, &workd[ipntr[1] - 1], &c__1);

	sgttrs_("Notranspose", &n, &c__1, adl, ad, adu, adu2, ipiv, &workd[
		ipntr[1] - 1], &n, &ierr, (ftnlen)11);
	if (ierr != 0) {
	    s_wsle(&io___34);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___35);
	    do_lio(&c__9, &c__1, " Error with _gttrs in _SDRV2. ", (ftnlen)30)
		    ;
	    e_wsle();
	    s_wsle(&io___36);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    goto L9000;
	}

/*           %-----------------------------------------% */
/*           | L O O P   B A C K to call SSAUPD again. | */
/*           %-----------------------------------------% */

	goto L10;

    }

/*     %----------------------------------------% */
/*     | Either we have convergence or there is | */
/*     | an error.                              | */
/*     %----------------------------------------% */

    if (info < 0) {

/*        %----------------------------% */
/*        | Error message.  Check the  | */
/*        | documentation in SSAUPD    | */
/*        %----------------------------% */

	s_wsle(&io___37);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___38);
	do_lio(&c__9, &c__1, " Error with _saupd, info = ", (ftnlen)27);
	do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___39);
	do_lio(&c__9, &c__1, " Check documentation of _saupd ", (ftnlen)31);
	e_wsle();
	s_wsle(&io___40);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();

    } else {

/*        %-------------------------------------------% */
/*        | No fatal errors occurred.                 | */
/*        | Post-Process using SSEUPD.                | */
/*        |                                           | */
/*        | Computed eigenvalues may be extracted.    | */
/*        |                                           | */
/*        | Eigenvectors may also be computed now if  | */
/*        | desired.  (indicated by rvec = .true.)    | */
/*        %-------------------------------------------% */

	rvec = TRUE_;

	sseupd_(&rvec, "All", select, d__, v, &c__256, &sigma, bmat, &n, 
		which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, 
		workd, workl, &lworkl, &ierr, (ftnlen)3, (ftnlen)1, (ftnlen)2)
		;

/*        %----------------------------------------------% */
/*        | Eigenvalues are returned in the first column | */
/*        | of the two dimensional array D and the       | */
/*        | corresponding eigenvectors are returned in   | */
/*        | the first NEV columns of the two dimensional | */
/*        | array V if requested.  Otherwise, an         | */
/*        | orthogonal basis for the invariant subspace  | */
/*        | corresponding to the eigenvalues in D is     | */
/*        | returned in V.                               | */
/*        %----------------------------------------------% */
	if (ierr != 0) {

/*           %------------------------------------% */
/*           | Error condition:                   | */
/*           | Check the documentation of SSEUPD. | */
/*           %------------------------------------% */

	    s_wsle(&io___44);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___45);
	    do_lio(&c__9, &c__1, " Error with _seupd, info = ", (ftnlen)27);
	    do_lio(&c__3, &c__1, (char *)&ierr, (ftnlen)sizeof(integer));
	    e_wsle();
	    s_wsle(&io___46);
	    do_lio(&c__9, &c__1, " Check the documentation of _seupd ", (
		    ftnlen)35);
	    e_wsle();
	    s_wsle(&io___47);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();

	} else {

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
/*              | eigenvectors.  (iparam(5) | */
/*              | indicates how many are    | */
/*              | accurate to the requested | */
/*              | tolerance)                | */
/*              %---------------------------% */

		av_(&n, &v[(j << 8) - 256], ax);
		r__1 = -d__[j - 1];
		saxpy_(&n, &r__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
		d__[j + 24] = snrm2_(&n, ax, &c__1);
		d__[j + 24] /= (r__1 = d__[j - 1], dabs(r__1));

/* L30: */
	    }

/*           %-------------------------------% */
/*           | Display computed residuals    | */
/*           %-------------------------------% */

	    smout_(&c__6, &nconv, &c__2, d__, &c__25, &c_n6, "Ritz values an"
		    "d relative residuals", (ftnlen)34);
	}

/*        %------------------------------------------% */
/*        | Print additional convergence information | */
/*        %------------------------------------------% */

	if (info == 1) {
	    s_wsle(&io___50);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___51);
	    do_lio(&c__9, &c__1, " Maximum number of iterations reached.", (
		    ftnlen)38);
	    e_wsle();
	    s_wsle(&io___52);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	} else if (info == 3) {
	    s_wsle(&io___53);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___54);
	    do_lio(&c__9, &c__1, " No shifts could be applied during implicit"
		    , (ftnlen)43);
	    do_lio(&c__9, &c__1, " Arnoldi update, try increasing NCV.", (
		    ftnlen)36);
	    e_wsle();
	    s_wsle(&io___55);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	}

	s_wsle(&io___56);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___57);
	do_lio(&c__9, &c__1, " _SDRV2 ", (ftnlen)8);
	e_wsle();
	s_wsle(&io___58);
	do_lio(&c__9, &c__1, " ====== ", (ftnlen)8);
	e_wsle();
	s_wsle(&io___59);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___60);
	do_lio(&c__9, &c__1, " Size of the matrix is ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___61);
	do_lio(&c__9, &c__1, " The number of Ritz values requested is ", (
		ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&nev, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___62);
	do_lio(&c__9, &c__1, " The number of Arnoldi vectors generated", (
		ftnlen)40);
	do_lio(&c__9, &c__1, " (NCV) is ", (ftnlen)10);
	do_lio(&c__3, &c__1, (char *)&ncv, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___63);
	do_lio(&c__9, &c__1, " What portion of the spectrum: ", (ftnlen)31);
	do_lio(&c__9, &c__1, which, (ftnlen)2);
	e_wsle();
	s_wsle(&io___64);
	do_lio(&c__9, &c__1, " The number of converged Ritz values is ", (
		ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&nconv, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___65);
	do_lio(&c__9, &c__1, " The number of Implicit Arnoldi update", (
		ftnlen)38);
	do_lio(&c__9, &c__1, " iterations taken is ", (ftnlen)21);
	do_lio(&c__3, &c__1, (char *)&iparam[2], (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___66);
	do_lio(&c__9, &c__1, " The number of OP*x is ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&iparam[8], (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___67);
	do_lio(&c__9, &c__1, " The convergence criterion is ", (ftnlen)30);
	do_lio(&c__4, &c__1, (char *)&tol, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___68);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();

    }

/*     %---------------------------% */
/*     | Done with program ssdrv2. | */
/*     %---------------------------% */

L9000:

    return 0;
} /* MAIN__ */

/* ------------------------------------------------------------------------ */
/*     Matrix vector subroutine */
/*     where the matrix is the 1 dimensional discrete Laplacian on */
/*     the interval [0,1] with zero Dirichlet boundary condition. */

int av_(integer *n, real *v, real *w)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    integer j;
    real h2;

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    w[1] = v[1] * 2.f - v[2];
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j) {
	w[j] = -v[j - 1] + v[j] * 2.f - v[j + 1];
/* L100: */
    }
    j = *n;
    w[j] = -v[j - 1] + v[j] * 2.f;

/*     Scale the vector w by (1 / h^2). */

    h2 = 1.f / (real) ((*n + 1) * (*n + 1));
    r__1 = 1.f / h2;
    sscal_(n, &r__1, &w[1], &c__1);
    return 0;
} /* av_ */

