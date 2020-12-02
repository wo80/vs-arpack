/* EXAMPLES\NONSYM\sndrv3.f -- translated by f2c (version 20100827). */

#include "arpack.h"

int sndrv3()
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    real d__[75]	/* was [25][3] */, h__;
    integer j, n;
    real v[6400]	/* was [256][25] */, md[256], me[255];
    real ax[256];
    real mx[256];
    integer ido, ncv, nev;
    real tol;
    char bmat[1];
    integer mode, info;
    logical rvec;
    integer ierr;
    char which[2];
    real resid[256];
    integer nconv;
    real workd[768];
    logical first;
    integer ipntr[14];
    real workl[2025];
    integer iparam[11];
    real sigmai;
    logical select[25];
    real sigmar;
    integer ishfts, maxitr, lworkl;
    real workev[75];

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 6, 0, 0, 0 };
    static cilist io___5 = { 0, 6, 0, 0, 0 };
    static cilist io___6 = { 0, 6, 0, 0, 0 };
    static cilist io___14 = { 0, 6, 0, 0, 0 };
    static cilist io___15 = { 0, 6, 0, 0, 0 };
    static cilist io___16 = { 0, 6, 0, 0, 0 };
    static cilist io___30 = { 0, 6, 0, 0, 0 };
    static cilist io___31 = { 0, 6, 0, 0, 0 };
    static cilist io___32 = { 0, 6, 0, 0, 0 };
    static cilist io___33 = { 0, 6, 0, 0, 0 };
    static cilist io___34 = { 0, 6, 0, 0, 0 };
    static cilist io___35 = { 0, 6, 0, 0, 0 };
    static cilist io___36 = { 0, 6, 0, 0, 0 };
    static cilist io___43 = { 0, 6, 0, 0, 0 };
    static cilist io___44 = { 0, 6, 0, 0, 0 };
    static cilist io___45 = { 0, 6, 0, 0, 0 };
    static cilist io___46 = { 0, 6, 0, 0, 0 };
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
    static cilist io___69 = { 0, 6, 0, 0, 0 };

/*     Simple program to illustrate the idea of reverse communication */
/*     in inverse mode for a generalized nonsymmetric eigenvalue problem. */

/*     We implement example three of ex-nonsym.doc in DOCUMENTS directory */

/* \Example-3 */
/*     ... Suppose we want to solve A*x = lambda*B*x in inverse mode, */
/*         where A and B are derived from the finite element discretization */
/*         of the 1-dimensional convection-diffusion operator */
/*                           (d^2u / dx^2) + rho*(du/dx) */
/*         on the interval [0,1] with zero Dirichlet boundary condition */
/*         using linear elements. */

/*     ... So OP = inv[M]*A  and  B = M. */

/*     ... Use mode 2 of SNAUPD. */

/* \BeginLib */

/* \Routines called: */
/*     snaupd  ARPACK reverse communication interface routine. */
/*     sneupd  ARPACK routine that returns Ritz values and (optionally) */
/*             Ritz vectors. */
/*     spttrf  LAPACK symmetric positive definite tridiagonal factorization */
/*             routine. */
/*     spttrs  LAPACK symmetric positive definite tridiagonal solve routine. */
/*     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully. */
/*     saxpy   Level 1 BLAS that computes y <- alpha*x+y. */
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
/* FILE: ndrv3.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2 */

/* \Remarks */
/*     1. None */

/* \EndLib */
/* -------------------------------------------------------------------------- */

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
/*     | generalized eigenvalue problem is solved (BMAT =   | */
/*     | 'G').  NEV is the number of eigenvalues to be      | */
/*     | approximated.  The user can modify NEV, NCV, WHICH | */
/*     | to solve problems of different sizes, and to get   | */
/*     | different parts of the spectrum.  However, The     | */
/*     | following conditions must be satisfied:            | */
/*     |                    N <= MAXN,                      | */
/*     |                  NEV <= MAXNEV,                    | */
/*     |              NEV + 2 <= NCV <= MAXNCV              | */
/*     %----------------------------------------------------% */

    n = 100;
    nev = 4;
    ncv = 20;
    if (n > 256) {
	s_wsle(&io___4);
	do_lio(&c__9, &c__1, " ERROR with _NDRV3: N is greater than MAXN ", (
		ftnlen)43);
	e_wsle();
	goto L9000;
    } else if (nev > 10) {
	s_wsle(&io___5);
	do_lio(&c__9, &c__1, " ERROR with _NDRV3: NEV is greater than MAXNEV "
		, (ftnlen)47);
	e_wsle();
	goto L9000;
    } else if (ncv > 25) {
	s_wsle(&io___6);
	do_lio(&c__9, &c__1, " ERROR with _NDRV3: NCV is greater than MAXNCV "
		, (ftnlen)47);
	e_wsle();
	goto L9000;
    }
    *(unsigned char *)bmat = 'G';
    s_copy(which, "LM", (ftnlen)2, (ftnlen)2);

/*     %------------------------------------------------% */
/*     | M is the mass matrix formed by using piecewise | */
/*     | linear elements on [0,1].                      | */
/*     %------------------------------------------------% */

    h__ = 1.f / (real) (n + 1);
    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j) {
	md[j - 1] = h__ * 4.f;
	me[j - 1] = h__ * 1.f;
/* L20: */
    }
    md[n - 1] = h__ * 4.f;

    spttrf_(&n, md, me, &ierr);
    if (ierr != 0) {
	s_wsle(&io___14);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___15);
	do_lio(&c__9, &c__1, " ERROR with _pttrf. ", (ftnlen)20);
	e_wsle();
	s_wsle(&io___16);
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
/*     | iterations allowed.  Mode 2 of SNAUPD is used     | */
/*     | (IPARAM(7) = 2).  All these options can be        | */
/*     | changed by the user. For details, see the         | */
/*     | documentation in SNAUPD.                          | */
/*     %---------------------------------------------------% */

    ishfts = 1;
    maxitr = 300;
    mode = 2;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

/*     %-------------------------------------------% */
/*     | M A I N   L O O P (Reverse communication) | */
/*     %-------------------------------------------% */

L10:

/*        %---------------------------------------------% */
/*        | Repeatedly call the routine SNAUPD and take | */
/*        | actions indicated by parameter IDO until    | */
/*        | either convergence is indicated or maxitr   | */
/*        | has been exceeded.                          | */
/*        %---------------------------------------------% */

    snaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, 
	    iparam, ipntr, workd, workl, &lworkl, &info, (ftnlen)1, (ftnlen)2)
	    ;

    if (ido == -1 || ido == 1) {

/*           %----------------------------------------% */
/*           | Perform  y <--- OP*x = inv[M]*A*x      | */
/*           | The user should supply his/her own     | */
/*           | matrix vector routine and a linear     | */
/*           | system solver.  The matrix-vector      | */
/*           | subroutine should take workd(ipntr(1)) | */
/*           | as input, and the final result should  | */
/*           | be returned to workd(ipntr(2)).        | */
/*           %----------------------------------------% */

	sndrv3_av_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
	spttrs_(&n, &c__1, md, me, &workd[ipntr[1] - 1], &n, &ierr);
	if (ierr != 0) {
	    s_wsle(&io___30);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___31);
	    do_lio(&c__9, &c__1, " ERROR with _pttrs. ", (ftnlen)20);
	    e_wsle();
	    s_wsle(&io___32);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    goto L9000;
	}

/*           %-----------------------------------------% */
/*           | L O O P   B A C K to call SNAUPD again. | */
/*           %-----------------------------------------% */

	goto L10;

    } else if (ido == 2) {

/*           %-------------------------------------% */
/*           |        Perform  y <--- M*x          | */
/*           | The matrix vector multiplication    | */
/*           | routine should take workd(ipntr(1)) | */
/*           | as input and return the result to   | */
/*           | workd(ipntr(2)).                    | */
/*           %-------------------------------------% */

	sndrv3_mv_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

/*           %-----------------------------------------% */
/*           | L O O P   B A C K to call SNAUPD again. | */
/*           %-----------------------------------------% */

	goto L10;

    }

/*     %-----------------------------------------% */
/*     | Either we have convergence, or there is | */
/*     | an error.                               | */
/*     %-----------------------------------------% */

    if (info < 0) {

/*        %--------------------------% */
/*        | Error message. Check the | */
/*        | documentation in SNAUPD. | */
/*        %--------------------------% */

	s_wsle(&io___33);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___34);
	do_lio(&c__9, &c__1, " Error with _naupd, info = ", (ftnlen)27);
	do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___35);
	do_lio(&c__9, &c__1, " Check the documentation of _naupd.", (ftnlen)
		35);
	e_wsle();
	s_wsle(&io___36);
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

	    s_wsle(&io___43);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___44);
	    do_lio(&c__9, &c__1, " Error with _neupd, info = ", (ftnlen)27);
	    do_lio(&c__3, &c__1, (char *)&ierr, (ftnlen)sizeof(integer));
	    e_wsle();
	    s_wsle(&io___45);
	    do_lio(&c__9, &c__1, " Check the documentation of _neupd", (
		    ftnlen)34);
	    e_wsle();
	    s_wsle(&io___46);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();

	} else {

	    first = TRUE_;
	    nconv = iparam[4];
	    i__1 = iparam[4];
	    for (j = 1; j <= i__1; ++j) {

/*              %---------------------------% */
/*              | Compute the residual norm | */
/*              |                           | */
/*              |  ||  A*x - lambda*M*x ||  | */
/*              |                           | */
/*              | for the NCONV accurately  | */
/*              | computed eigenvalues and  | */
/*              | eigenvectors.  (iparam(5) | */
/*              | indicates how many are    | */
/*              | accurate to the requested | */
/*              | tolerance)                | */
/*              %---------------------------% */

		if (d__[j + 24] == 0.f) {

/*                 %--------------------% */
/*                 | Ritz value is real | */
/*                 %--------------------% */

			sndrv3_av_(&n, &v[(j << 8) - 256], ax);
			sndrv3_mv_(&n, &v[(j << 8) - 256], mx);
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

			sndrv3_av_(&n, &v[(j << 8) - 256], ax);
			sndrv3_mv_(&n, &v[(j << 8) - 256], mx);
		    r__1 = -d__[j - 1];
		    saxpy_(&n, &r__1, mx, &c__1, ax, &c__1);
			sndrv3_mv_(&n, &v[(j + 1 << 8) - 256], mx);
		    saxpy_(&n, &d__[j + 24], mx, &c__1, ax, &c__1);
/* Computing 2nd power */
		    r__1 = snrm2_(&n, ax, &c__1);
		    d__[j + 49] = r__1 * r__1;
			sndrv3_av_(&n, &v[(j + 1 << 8) - 256], ax);
			sndrv3_mv_(&n, &v[(j + 1 << 8) - 256], mx);
		    r__1 = -d__[j - 1];
		    saxpy_(&n, &r__1, mx, &c__1, ax, &c__1);
			sndrv3_mv_(&n, &v[(j << 8) - 256], mx);
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

/*           %-----------------------------% */
/*           | Display computed residuals. | */
/*           %-----------------------------% */

	    smout_(&c__6, &nconv, &c__3, d__, &c__25, &c_n6, "Ritz values (R"
		    "eal,Imag) and relative residuals", (ftnlen)46);

	}

/*        %------------------------------------------% */
/*        | Print additional convergence information | */
/*        %------------------------------------------% */

	if (info == 1) {
	    s_wsle(&io___51);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___52);
	    do_lio(&c__9, &c__1, " Maximum number of iterations reached.", (
		    ftnlen)38);
	    e_wsle();
	    s_wsle(&io___53);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	} else if (info == 3) {
	    s_wsle(&io___54);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___55);
	    do_lio(&c__9, &c__1, " No shifts could be applied during implicit"
		    , (ftnlen)43);
	    do_lio(&c__9, &c__1, " Arnoldi update, try increasing NCV.", (
		    ftnlen)36);
	    e_wsle();
	    s_wsle(&io___56);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	}

	s_wsle(&io___57);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___58);
	do_lio(&c__9, &c__1, " _NDRV3 ", (ftnlen)8);
	e_wsle();
	s_wsle(&io___59);
	do_lio(&c__9, &c__1, " ====== ", (ftnlen)8);
	e_wsle();
	s_wsle(&io___60);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___61);
	do_lio(&c__9, &c__1, " Size of the matrix is ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___62);
	do_lio(&c__9, &c__1, " The number of Ritz values requested is ", (
		ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&nev, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___63);
	do_lio(&c__9, &c__1, " The number of Arnoldi vectors generated", (
		ftnlen)40);
	do_lio(&c__9, &c__1, " (NCV) is ", (ftnlen)10);
	do_lio(&c__3, &c__1, (char *)&ncv, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___64);
	do_lio(&c__9, &c__1, " What portion of the spectrum: ", (ftnlen)31);
	do_lio(&c__9, &c__1, which, (ftnlen)2);
	e_wsle();
	s_wsle(&io___65);
	do_lio(&c__9, &c__1, " The number of converged Ritz values is ", (
		ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&nconv, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___66);
	do_lio(&c__9, &c__1, " The number of Implicit Arnoldi update", (
		ftnlen)38);
	do_lio(&c__9, &c__1, " iterations taken is ", (ftnlen)21);
	do_lio(&c__3, &c__1, (char *)&iparam[2], (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___67);
	do_lio(&c__9, &c__1, " The number of OP*x is ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&iparam[8], (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___68);
	do_lio(&c__9, &c__1, " The convergence criterion is ", (ftnlen)30);
	do_lio(&c__4, &c__1, (char *)&tol, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___69);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();

    }

/*     %---------------------------% */
/*     | Done with program sndrv3. | */
/*     %---------------------------% */

L9000:

    return 0;
} /* MAIN__ */

/* ========================================================================== */

/*     matrix vector multiplication subroutine */

int sndrv3_av_(integer *n, real *v, real *w)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    real h__;
    integer j;
    real s, dd, dl, du;

/*     Compute the matrix vector multiplication y<---A*x */
/*     where A is stiffness matrix obtained from the finite element */
/*     discretization of the 1-dimensional convection diffusion operator */
/*                           d^2u/dx^2 + rho*(du/dx) */
/*     on the interval [0,1] with zero Dirichlet boundary condition using */
/*     linear elements. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    h__ = 1.f / (real) (*n + 1);
    s = 5.f;
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

/* ------------------------------------------------------------------------ */
int sndrv3_mv_(integer *n, real *v, real *w)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    real h__;
    integer j;

/*     Compute the matrix vector multiplication y<---M*x */
/*     where M is the mass matrix formed by using piecewise linear */
/*     elements on [0,1]. */

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

    h__ = 1.f / (real) (*n + 1);
    sscal_(n, &h__, &w[1], &c__1);
    return 0;
} /* mv_ */

