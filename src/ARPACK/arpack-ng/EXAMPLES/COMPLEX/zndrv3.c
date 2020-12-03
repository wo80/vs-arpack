/* EXAMPLES\COMPLEX\zndrv3.f -- translated by f2c (version 20100827). */

#include "arpack.h"

int zndrv3()
{
    /* System generated locals */
    int32_t i__1, i__2;
    zomplex z__1, z__2;

    void z_div(zomplex *, zomplex *, zomplex *);
    double d_imag(zomplex *);

    /* Local variables */
    zomplex d__[25], h__;
    int32_t j, n;
    zomplex v[6400]	/* was [256][25] */, dd[256], dl[256];
    double rd[75]	/* was [25][3] */;
    zomplex ax[256], du[256];
    zomplex mx[256], du2[256];
    int32_t ido, ncv, nev;
    double tol;
    char bmat[1];
    int32_t mode, info;
    bool rvec;
    int32_t ierr, ipiv[256];
    zomplex sigma;
    char which[2];
    zomplex resid[256];
    int32_t nconv;
    zomplex workd[768];
    int32_t ipntr[14];
    zomplex workl[2000];
    double rwork[256];
    int32_t iparam[11];
    bool select[25];
    int32_t ishfts;
    int32_t maxitr;
    int32_t lworkl;
    zomplex workev[50];

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 6, 0, 0, 0 };
    static cilist io___5 = { 0, 6, 0, 0, 0 };
    static cilist io___6 = { 0, 6, 0, 0, 0 };
    static cilist io___18 = { 0, 6, 0, 0, 0 };
    static cilist io___19 = { 0, 6, 0, 0, 0 };
    static cilist io___20 = { 0, 6, 0, 0, 0 };
    static cilist io___35 = { 0, 6, 0, 0, 0 };
    static cilist io___36 = { 0, 6, 0, 0, 0 };
    static cilist io___37 = { 0, 6, 0, 0, 0 };
    static cilist io___38 = { 0, 6, 0, 0, 0 };
    static cilist io___39 = { 0, 6, 0, 0, 0 };
    static cilist io___40 = { 0, 6, 0, 0, 0 };
    static cilist io___41 = { 0, 6, 0, 0, 0 };
    static cilist io___46 = { 0, 6, 0, 0, 0 };
    static cilist io___47 = { 0, 6, 0, 0, 0 };
    static cilist io___48 = { 0, 6, 0, 0, 0 };
    static cilist io___49 = { 0, 6, 0, 0, 0 };
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
    static cilist io___70 = { 0, 6, 0, 0, 0 };
    static cilist io___71 = { 0, 6, 0, 0, 0 };
    static cilist io___72 = { 0, 6, 0, 0, 0 };

/*     Simple program to illustrate the idea of reverse communication */
/*     in inverse mode for a generalized complex nonsymmetric eigenvalue */
/*     problem. */

/*     We implement example three of ex-complex.doc in DOCUMENTS directory */

/* \Example-3 */
/*     ... Suppose we want to solve A*x = lambda*B*x in regular mode, */
/*         where A and B are derived from the finite element discretization */
/*         of the 1-dimensional convection-diffusion operator */
/*                   (d^2u/dx^2) + rho*(du/dx) */
/*         on the interval [0,1] with zero boundary condition using */
/*         piecewise linear elements. */

/*     ... OP = inv[M]*A  and  B = M. */

/*     ... Use mode 2 of ZNAUPD . */

/* \BeginLib */

/* \Routines called: */
/*     znaupd   ARPACK reverse communication interface routine. */
/*     zneupd   ARPACK routine that returns Ritz values and (optionally) */
/*             Ritz vectors. */
/*     zgttrf   LAPACK tridiagonal factorization routine. */
/*     zgttrs   LAPACK tridiagonal solve routine. */
/*     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully. */
/*     zaxpy    Level 1 BLAS that computes y <- alpha*x+y. */
/*     dznrm2   Level 1 BLAS that computes the norm of a vector. */
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
/* FILE: ndrv3.F   SID: 2.4   DATE OF SID: 10/18/00   RELEASE: 2 */

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
    sigma.r = 0., sigma.i = 0.;

/*     %-----------------------------------------------------% */
/*     | The matrix M is chosen to be the symmetric tri-     | */
/*     | diagonal matrix with 4 on the diagonal and 1 on the | */
/*     | off diagonals. It is factored by LAPACK subroutine  | */
/*     | zgttrf .                                             | */
/*     %-----------------------------------------------------% */

    i__1 = n + 1;
    z__2.r = (double) i__1, z__2.i = 0.;
    z_div(&z__1, &c_b137, &z__2);
    h__.r = z__1.r, h__.i = z__1.i;
    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j - 1;
	z__1.r = h__.r * 1. - h__.i * 0., z__1.i = h__.i * 1. + h__.r * 0.;
	dl[i__2].r = z__1.r, dl[i__2].i = z__1.i;
	i__2 = j - 1;
	z__1.r = h__.r * 4. - h__.i * 0., z__1.i = h__.r * 0. + h__.i * 4.;
	dd[i__2].r = z__1.r, dd[i__2].i = z__1.i;
	i__2 = j - 1;
	z__1.r = h__.r * 1. - h__.i * 0., z__1.i = h__.i * 1. + h__.r * 0.;
	du[i__2].r = z__1.r, du[i__2].i = z__1.i;
/* L20: */
    }
    i__1 = n - 1;
    z__1.r = h__.r * 4. - h__.i * 0., z__1.i = h__.r * 0. + h__.i * 4.;
    dd[i__1].r = z__1.r, dd[i__1].i = z__1.i;

    zgttrf_(&n, dl, dd, du, du2, ipiv, &ierr);
    if (ierr != 0) {
	s_wsle(&io___18);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___19);
	do_lio(&c__9, &c__1, " ERROR with _gttrf. ", (ftnlen)20);
	e_wsle();
	s_wsle(&io___20);
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
/*     | iterations allowed.  Mode 2 of ZNAUPD  is used     | */
/*     | (IPARAM(7) = 2).  All these options can be        | */
/*     | changed by the user. For details, see the         | */
/*     | documentation in ZNAUPD .                          | */
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
/*        | Repeatedly call the routine ZNAUPD  and take | */
/*        | actions indicated by parameter IDO until    | */
/*        | either convergence is indicated or maxitr   | */
/*        | has been exceeded.                          | */
/*        %---------------------------------------------% */

    znaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, 
	    iparam, ipntr, workd, workl, &lworkl, rwork, &info, (ftnlen)1, (
	    ftnlen)2);

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

	zndrv3_av_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
	zgttrs_("N", &n, &c__1, dl, dd, du, du2, ipiv, &workd[ipntr[1] - 1], &
		n, &ierr, (ftnlen)1);
	if (ierr != 0) {
	    s_wsle(&io___35);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___36);
	    do_lio(&c__9, &c__1, " ERROR with _gttrs. ", (ftnlen)20);
	    e_wsle();
	    s_wsle(&io___37);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    goto L9000;
	}

/*           %-----------------------------------------% */
/*           | L O O P   B A C K to call ZNAUPD  again. | */
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

	zndrv3_mv_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

/*           %-----------------------------------------% */
/*           | L O O P   B A C K to call ZNAUPD  again. | */
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
/*        | documentation in ZNAUPD . | */
/*        %--------------------------% */

	s_wsle(&io___38);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___39);
	do_lio(&c__9, &c__1, " Error with _naupd, info = ", (ftnlen)27);
	do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(int32_t));
	e_wsle();
	s_wsle(&io___40);
	do_lio(&c__9, &c__1, " Check the documentation of _naupd.", (ftnlen)
		35);
	e_wsle();
	s_wsle(&io___41);
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

	rvec = true;

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

	    s_wsle(&io___46);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___47);
	    do_lio(&c__9, &c__1, " Error with _neupd, info = ", (ftnlen)27);
	    do_lio(&c__3, &c__1, (char *)&ierr, (ftnlen)sizeof(int32_t));
	    e_wsle();
	    s_wsle(&io___48);
	    do_lio(&c__9, &c__1, " Check the documentation of _neupd", (
		    ftnlen)34);
	    e_wsle();
	    s_wsle(&io___49);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();

	} else {

	    nconv = iparam[4];
	    i__1 = nconv;
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

		zndrv3_av_(&n, &v[(j << 8) - 256], ax);
		zndrv3_mv_(&n, &v[(j << 8) - 256], mx);
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

/*           %-----------------------------% */
/*           | Display computed residuals. | */
/*           %-----------------------------% */

	    dmout_(&c__6, &nconv, &c__3, rd, &c__25, &c_n6, "Ritz values (Re"
		    "al, Imag) and relative residuals", (ftnlen)47);

	}

/*        %------------------------------------------% */
/*        | Print additional convergence information | */
/*        %------------------------------------------% */

	if (info == 1) {
	    s_wsle(&io___54);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___55);
	    do_lio(&c__9, &c__1, " Maximum number of iterations reached.", (
		    ftnlen)38);
	    e_wsle();
	    s_wsle(&io___56);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	} else if (info == 3) {
	    s_wsle(&io___57);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___58);
	    do_lio(&c__9, &c__1, " No shifts could be applied during implicit"
		    , (ftnlen)43);
	    do_lio(&c__9, &c__1, " Arnoldi update, try increasing NCV.", (
		    ftnlen)36);
	    e_wsle();
	    s_wsle(&io___59);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	}

	s_wsle(&io___60);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___61);
	do_lio(&c__9, &c__1, "_NDRV3 ", (ftnlen)7);
	e_wsle();
	s_wsle(&io___62);
	do_lio(&c__9, &c__1, "====== ", (ftnlen)7);
	e_wsle();
	s_wsle(&io___63);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___64);
	do_lio(&c__9, &c__1, " Size of the matrix is ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(int32_t));
	e_wsle();
	s_wsle(&io___65);
	do_lio(&c__9, &c__1, " The number of Ritz values requested is ", (
		ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&nev, (ftnlen)sizeof(int32_t));
	e_wsle();
	s_wsle(&io___66);
	do_lio(&c__9, &c__1, " The number of Arnoldi vectors generated ", (
		ftnlen)41);
	do_lio(&c__9, &c__1, " (NCV) is ", (ftnlen)10);
	do_lio(&c__3, &c__1, (char *)&ncv, (ftnlen)sizeof(int32_t));
	e_wsle();
	s_wsle(&io___67);
	do_lio(&c__9, &c__1, " What portion of the spectrum: ", (ftnlen)31);
	do_lio(&c__9, &c__1, which, (ftnlen)2);
	e_wsle();
	s_wsle(&io___68);
	do_lio(&c__9, &c__1, " The number of converged Ritz values is ", (
		ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&nconv, (ftnlen)sizeof(int32_t));
	e_wsle();
	s_wsle(&io___69);
	do_lio(&c__9, &c__1, " The number of Implicit Arnoldi update", (
		ftnlen)38);
	do_lio(&c__9, &c__1, " iterations taken is ", (ftnlen)21);
	do_lio(&c__3, &c__1, (char *)&iparam[2], (ftnlen)sizeof(int32_t));
	e_wsle();
	s_wsle(&io___70);
	do_lio(&c__9, &c__1, " The number of OP*x is ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&iparam[8], (ftnlen)sizeof(int32_t));
	e_wsle();
	s_wsle(&io___71);
	do_lio(&c__9, &c__1, " The convergence criterion is ", (ftnlen)30);
	do_lio(&c__5, &c__1, (char *)&tol, (ftnlen)sizeof(double));
	e_wsle();
	s_wsle(&io___72);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();

    }

L9000:

    return 0;
} /* MAIN__ */

/* ========================================================================== */

/*     matrix vector multiplication subroutine */

int zndrv3_av_(int32_t *n, zomplex *v, zomplex *w)
{
    /* System generated locals */
    int32_t i__1, i__2, i__3, i__4, i__5;
    zomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    void z_div(zomplex *, zomplex *, zomplex *);

    /* Local variables */
    zomplex h__;
    int32_t j;
    zomplex s, dd, dl, du;

/*     Compute the matrix vector multiplication y<---A*x */
/*     where A is the stiffness matrix formed by using piecewise linear */
/*     elements on [0,1]. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    i__1 = *n + 1;
    z__2.r = (double) i__1, z__2.i = 0.;
    z_div(&z__1, &c_b137, &z__2);
    h__.r = z__1.r, h__.i = z__1.i;
    z_div(&z__1, &c_b164_dx, &c_b3_dx);
    s.r = z__1.r, s.i = z__1.i;
    z_div(&z__1, &c_b3_dx, &h__);
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
/* L10: */
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

/* ------------------------------------------------------------------------ */
int zndrv3_mv_(int32_t *n, zomplex *v, zomplex *w)
{
    /* System generated locals */
    int32_t i__1, i__2, i__3, i__4, i__5;
    zomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    void z_div(zomplex *, zomplex *, zomplex *);

    /* Local variables */
    zomplex h__;
    int32_t j;

/*     Compute the matrix vector multiplication y<---M*x */
/*     where M is the mass matrix formed by using piecewise linear elements */
/*     on [0,1]. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    z__2.r = v[1].r * 4. - v[1].i * 0., z__2.i = v[1].i * 4. + v[1].r * 0.;
    z__3.r = v[2].r * 1. - v[2].i * 0., z__3.i = v[2].i * 1. + v[2].r * 0.;
    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
    w[1].r = z__1.r, w[1].i = z__1.i;
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j - 1;
	z__3.r = v[i__3].r * 1. - v[i__3].i * 0., z__3.i = v[i__3].i * 1. + v[
		i__3].r * 0.;
	i__4 = j;
	z__4.r = v[i__4].r * 4. - v[i__4].i * 0., z__4.i = v[i__4].i * 4. + v[
		i__4].r * 0.;
	z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
	i__5 = j + 1;
	z__5.r = v[i__5].r * 1. - v[i__5].i * 0., z__5.i = v[i__5].i * 1. + v[
		i__5].r * 0.;
	z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
	w[i__2].r = z__1.r, w[i__2].i = z__1.i;
/* L10: */
    }
    i__1 = *n;
    i__2 = *n - 1;
    z__2.r = v[i__2].r * 1. - v[i__2].i * 0., z__2.i = v[i__2].i * 1. + v[
	    i__2].r * 0.;
    i__3 = *n;
    z__3.r = v[i__3].r * 4. - v[i__3].i * 0., z__3.i = v[i__3].i * 4. + v[
	    i__3].r * 0.;
    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
    w[i__1].r = z__1.r, w[i__1].i = z__1.i;

    i__1 = *n + 1;
    z__2.r = (double) i__1, z__2.i = 0.;
    z_div(&z__1, &c_b137, &z__2);
    h__.r = z__1.r, h__.i = z__1.i;
    zscal_(n, &h__, &w[1], &c__1);
    return 0;
} /* mv_ */

