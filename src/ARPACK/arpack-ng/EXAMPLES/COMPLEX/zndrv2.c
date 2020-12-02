/* EXAMPLES\COMPLEX\zndrv2.f -- translated by f2c (version 20100827). */

#include "arpack.h"

struct {
    doublecomplex rho;
} convct_;

#define convct_1 convct_

int main()
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1, z__2, z__3, z__4;

    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
    double d_imag(doublecomplex *);

    /* Local variables */
    doublecomplex d__[25], h__;
    integer j, n;
    doublecomplex s, v[6400]	/* was [256][25] */, h2, s1, s2, s3, dd[256], 
	    dl[256];
    doublereal rd[75]	/* was [25][3] */;
    doublecomplex ax[256], du[256], du2[256];
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
    integer ipntr[14];
    doublecomplex workl[2000];
    doublereal rwork[256];
    integer iparam[11];
    logical select[25];
    integer ishfts;
    integer maxitr;
    integer lworkl;
    doublecomplex workev[50];

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 6, 0, 0, 0 };
    static cilist io___5 = { 0, 6, 0, 0, 0 };
    static cilist io___6 = { 0, 6, 0, 0, 0 };
    static cilist io___23 = { 0, 6, 0, 0, 0 };
    static cilist io___24 = { 0, 6, 0, 0, 0 };
    static cilist io___25 = { 0, 6, 0, 0, 0 };
    static cilist io___40 = { 0, 6, 0, 0, 0 };
    static cilist io___41 = { 0, 6, 0, 0, 0 };
    static cilist io___42 = { 0, 6, 0, 0, 0 };
    static cilist io___43 = { 0, 6, 0, 0, 0 };
    static cilist io___44 = { 0, 6, 0, 0, 0 };
    static cilist io___45 = { 0, 6, 0, 0, 0 };
    static cilist io___46 = { 0, 6, 0, 0, 0 };
    static cilist io___51 = { 0, 6, 0, 0, 0 };
    static cilist io___52 = { 0, 6, 0, 0, 0 };
    static cilist io___53 = { 0, 6, 0, 0, 0 };
    static cilist io___54 = { 0, 6, 0, 0, 0 };
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
    static cilist io___73 = { 0, 6, 0, 0, 0 };
    static cilist io___74 = { 0, 6, 0, 0, 0 };
    static cilist io___75 = { 0, 6, 0, 0, 0 };
    static cilist io___76 = { 0, 6, 0, 0, 0 };

/*     Simple program to illustrate the idea of reverse communication */
/*     in shift-invert mode for a standard complex nonsymmetric eigenvalue */
/*     problem. */

/*     We implement example two of ex-complex.doc in DOCUMENTS directory */

/* \Example-2 */
/*     ... Suppose we want to solve A*x = lambda*x in shift-invert mode, */
/*         where A is derived from the central difference discretization */
/*         of the 1-dimensional convection-diffusion operator */
/*                   (d^2u/dx^2) + rho*(du/dx) */
/*         on the interval [0,1] with zero Dirichlet boundary condition. */
/*     ... The shift sigma is a complex number. */

/*     ... OP = inv[A-sigma*I] and  B = I. */

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
/*     dznrm2   Level 1 BLAS that computes the norm of a vector. */
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
/* FILE: ndrv2.F   SID: 2.6   DATE OF SID: 10/18/00   RELEASE: 2 */

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

/*     %--------------------------------------------------% */
/*     | The number N is the dimension of the matrix.  A  | */
/*     | standard eigenvalue problem is solved (BMAT =    | */
/*     | 'I').  NEV is the number of eigenvalues (closest | */
/*     | to the shift SIGMA) to be approximated.  Since   | */
/*     | the shift-invert mode is used, WHICH is set to   | */
/*     | 'LM'.  The user can modify NEV, NCV, SIGMA to    | */
/*     | solve problems of different sizes, and to get    | */
/*     | different parts of the spectrum.  However, The   | */
/*     | following conditions must be satisfied:          | */
/*     |                 N <= MAXN,                       | */
/*     |               NEV <= MAXNEV,                     | */
/*     |           NEV + 2 <= NCV <= MAXNCV               | */
/*     %--------------------------------------------------% */

    n = 100;
    nev = 4;
    ncv = 20;
    if (n > 256) {
	s_wsle(&io___4);
	do_lio(&c__9, &c__1, " ERROR with _NDRV2: N is greater than MAXN ", (
		ftnlen)43);
	e_wsle();
	goto L9000;
    } else if (nev > 10) {
	s_wsle(&io___5);
	do_lio(&c__9, &c__1, " ERROR with _NDRV2: NEV is greater than MAXNEV "
		, (ftnlen)47);
	e_wsle();
	goto L9000;
    } else if (ncv > 25) {
	s_wsle(&io___6);
	do_lio(&c__9, &c__1, " ERROR with _NDRV2: NCV is greater than MAXNCV "
		, (ftnlen)47);
	e_wsle();
	goto L9000;
    }
    *(unsigned char *)bmat = 'I';
    s_copy(which, "LM", (ftnlen)2, (ftnlen)2);
    sigma.r = 0., sigma.i = 0.;

/*     %----------------------------------------------------% */
/*     | Construct C = A - SIGMA*I, factor C in complex     | */
/*     | arithmetic (using LAPACK subroutine zgttrf ). The   | */
/*     | matrix A is chosen to be the tridiagonal matrix    | */
/*     | derived from standard central difference of the    | */
/*     | 1-d convection diffusion operator - u``+ rho*u` on | */
/*     | the interval [0, 1] with zero Dirichlet boundary   | */
/*     | condition.                                         | */
/*     %----------------------------------------------------% */

    convct_1.rho.r = 10., convct_1.rho.i = 0.;
    i__1 = n + 1;
    z__2.r = (doublereal) i__1, z__2.i = 0.;
    z_div(&z__1, &c_b137, &z__2);
    h__.r = z__1.r, h__.i = z__1.i;
    z__1.r = h__.r * h__.r - h__.i * h__.i, z__1.i = h__.r * h__.i + h__.i * 
	    h__.r;
    h2.r = z__1.r, h2.i = z__1.i;
    z_div(&z__1, &convct_1.rho, &c_b3_dx);
    s.r = z__1.r, s.i = z__1.i;

    z__3.r = -1., z__3.i = -0.;
    z_div(&z__2, &z__3, &h2);
    z_div(&z__4, &s, &h__);
    z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - z__4.i;
    s1.r = z__1.r, s1.i = z__1.i;
    z_div(&z__2, &c_b3_dx, &h2);
    z__1.r = z__2.r - sigma.r, z__1.i = z__2.i - sigma.i;
    s2.r = z__1.r, s2.i = z__1.i;
    z__3.r = -1., z__3.i = -0.;
    z_div(&z__2, &z__3, &h2);
    z_div(&z__4, &s, &h__);
    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
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
	s_wsle(&io___23);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___24);
	do_lio(&c__9, &c__1, " ERROR with _gttrf in _NDRV2.", (ftnlen)29);
	e_wsle();
	s_wsle(&io___25);
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

/*     %-------------------------------------------% */
/*     | M A I N   L O O P (Reverse communication) | */
/*     %-------------------------------------------% */

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

    if (ido == -1 || ido == 1) {

/*           %-------------------------------------------% */
/*           | Perform  y <--- OP*x = inv[A-SIGMA*I]*x   | */
/*           | The user should supply his/her own linear | */
/*           | system solver here that takes             | */
/*           | workd(ipntr(1)) as the input, and returns | */
/*           | the result to workd(ipntr(2)).            | */
/*           %-------------------------------------------% */

	zcopy_(&n, &workd[ipntr[0] - 1], &c__1, &workd[ipntr[1] - 1], &c__1);

	zgttrs_("N", &n, &c__1, dl, dd, du, du2, ipiv, &workd[ipntr[1] - 1], &
		n, &ierr, (ftnlen)1);
	if (ierr != 0) {
	    s_wsle(&io___40);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___41);
	    do_lio(&c__9, &c__1, " ERROR with _gttrs in _NDRV2.", (ftnlen)29);
	    e_wsle();
	    s_wsle(&io___42);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    goto L9000;
	}

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

/*        %--------------------------% */
/*        | Error message, check the | */
/*        | documentation in ZNAUPD   | */
/*        %--------------------------% */

	s_wsle(&io___43);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___44);
	do_lio(&c__9, &c__1, " Error with _naupd, info = ", (ftnlen)27);
	do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___45);
	do_lio(&c__9, &c__1, " Check the documentation in _naupd.", (ftnlen)
		35);
	e_wsle();
	s_wsle(&io___46);
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

	    s_wsle(&io___51);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___52);
	    do_lio(&c__9, &c__1, " Error with _neupd, info = ", (ftnlen)27);
	    do_lio(&c__3, &c__1, (char *)&ierr, (ftnlen)sizeof(integer));
	    e_wsle();
	    s_wsle(&io___53);
	    do_lio(&c__9, &c__1, " Check the documentation of _neupd. ", (
		    ftnlen)36);
	    e_wsle();
	    s_wsle(&io___54);
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
		i__2 = j - 1;
		z__1.r = -d__[i__2].r, z__1.i = -d__[i__2].i;
		zaxpy_(&n, &z__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
		i__2 = j - 1;
		rd[j - 1] = d__[i__2].r;
		rd[j + 24] = d_imag(&d__[j - 1]);
		rd[j + 49] = dznrm2_(&n, ax, &c__1);
		rd[j + 49] /= dlapy2_(&rd[j - 1], &rd[j + 24]);
/* L60: */
	    }

/*           %-----------------------------% */
/*           | Display computed residuals. | */
/*           %-----------------------------% */

	    dmout_(&c__6, &nconv, &c__3, rd, &c__25, &c_n6, "Ritz values (Re"
		    "al, Imag) and relative residuals", (ftnlen)47);

	}

/*        %-------------------------------------------% */
/*        | Print additional convergence information. | */
/*        %-------------------------------------------% */

	if (info == 1) {
	    s_wsle(&io___58);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___59);
	    do_lio(&c__9, &c__1, " Maximum number of iterations reached.", (
		    ftnlen)38);
	    e_wsle();
	    s_wsle(&io___60);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	} else if (info == 3) {
	    s_wsle(&io___61);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___62);
	    do_lio(&c__9, &c__1, " No shifts could be applied during implicit"
		    , (ftnlen)43);
	    do_lio(&c__9, &c__1, " Arnoldi update, try increasing NCV.", (
		    ftnlen)36);
	    e_wsle();
	    s_wsle(&io___63);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	}

	s_wsle(&io___64);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___65);
	do_lio(&c__9, &c__1, "_NDRV2 ", (ftnlen)7);
	e_wsle();
	s_wsle(&io___66);
	do_lio(&c__9, &c__1, "====== ", (ftnlen)7);
	e_wsle();
	s_wsle(&io___67);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___68);
	do_lio(&c__9, &c__1, " Size of the matrix is ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___69);
	do_lio(&c__9, &c__1, " The number of Ritz values requested is ", (
		ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&nev, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___70);
	do_lio(&c__9, &c__1, " The number of Arnoldi vectors generated", (
		ftnlen)40);
	do_lio(&c__9, &c__1, " (NCV) is ", (ftnlen)10);
	do_lio(&c__3, &c__1, (char *)&ncv, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___71);
	do_lio(&c__9, &c__1, " What portion of the spectrum: ", (ftnlen)31);
	do_lio(&c__9, &c__1, which, (ftnlen)2);
	e_wsle();
	s_wsle(&io___72);
	do_lio(&c__9, &c__1, " The number of converged Ritz values is ", (
		ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&nconv, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___73);
	do_lio(&c__9, &c__1, " The number of Implicit Arnoldi update", (
		ftnlen)38);
	do_lio(&c__9, &c__1, " iterations taken is ", (ftnlen)21);
	do_lio(&c__3, &c__1, (char *)&iparam[2], (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___74);
	do_lio(&c__9, &c__1, " The number of OP*x is ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&iparam[8], (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___75);
	do_lio(&c__9, &c__1, " The convergence criterion is ", (ftnlen)30);
	do_lio(&c__5, &c__1, (char *)&tol, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___76);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();

    }

/*     %---------------------------% */
/*     | Done with program zndrv2 . | */
/*     %---------------------------% */

L9000:

    return 0;
} /* MAIN__ */

/* ------------------------------------------------------------------- */

/*     matrix vector multiplication subroutine */

int av_(integer *n, doublecomplex *v, doublecomplex *w)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    doublecomplex h__;
    integer j;
    doublecomplex s, h2, dd, dl, du;

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    i__1 = *n + 1;
    z__2.r = (doublereal) i__1, z__2.i = 0.;
    z_div(&z__1, &c_b137, &z__2);
    h__.r = z__1.r, h__.i = z__1.i;
    z__1.r = h__.r * h__.r - h__.i * h__.i, z__1.i = h__.r * h__.i + h__.i * 
	    h__.r;
    h2.r = z__1.r, h2.i = z__1.i;
    z_div(&z__1, &convct_1.rho, &c_b3_dx);
    s.r = z__1.r, s.i = z__1.i;
    z_div(&z__1, &c_b3_dx, &h2);
    dd.r = z__1.r, dd.i = z__1.i;
    z__3.r = -1., z__3.i = -0.;
    z_div(&z__2, &z__3, &h2);
    z_div(&z__4, &s, &h__);
    z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - z__4.i;
    dl.r = z__1.r, dl.i = z__1.i;
    z__3.r = -1., z__3.i = -0.;
    z_div(&z__2, &z__3, &h2);
    z_div(&z__4, &s, &h__);
    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
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

