/* EXAMPLES\COMPLEX\zndrv1.f -- translated by f2c (version 20100827). */

#include "arpack.h"

int zndrv1()
{
    /* System generated locals */
    int32_t i__1, i__2;
    zomplex z__1;

    double d_imag(zomplex *);

    /* Local variables */
    zomplex d__[30];
    int32_t j, n;
    zomplex v[7680]	/* was [256][30] */;
    double rd[90]	/* was [30][3] */;
    zomplex ax[256];
    int32_t nx, ido, ncv, nev;
    double tol;
    char bmat[1];
    int32_t mode, info;
    bool rvec;
    int32_t ierr;
    zomplex sigma;
    char which[2];
    zomplex resid[256];
    int32_t nconv;
    zomplex workd[768];
    int32_t ipntr[14];
    zomplex workl[2850];
    double rwork[30];
    int32_t iparam[11];
    bool select[30];
    int32_t ishfts;
    int32_t maxitr;
    int32_t lworkl;
    zomplex workev[90];

    /* Fortran I/O blocks */
    static cilist io___5 = { 0, 6, 0, 0, 0 };
    static cilist io___6 = { 0, 6, 0, 0, 0 };
    static cilist io___7 = { 0, 6, 0, 0, 0 };
    static cilist io___24 = { 0, 6, 0, 0, 0 };
    static cilist io___25 = { 0, 6, 0, 0, 0 };
    static cilist io___26 = { 0, 6, 0, 0, 0 };
    static cilist io___27 = { 0, 6, 0, 0, 0 };
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

/*     Example program to illustrate the idea of reverse communication */
/*     for a standard complex nonsymmetric eigenvalue problem. */

/*     We implement example one of ex-complex.doc in DOCUMENTS directory */

/* \Example-1 */
/*     ... Suppose we want to solve A*x = lambda*x in regular mode, */
/*         where A is obtained from the standard central difference */
/*         discretization of the convection-diffusion operator */
/*                 (Laplacian u) + rho*(du / dx) */
/*         on the unit squre [0,1]x[0,1] with zero Dirichlet boundary */
/*         condition. */

/*     ... OP = A  and  B = I. */

/*     ... Assume "call av (nx,x,y)" computes y = A*x */

/*     ... Use mode 1 of ZNAUPD . */
/**
 * \BeginLib
 *
 * \Routines called
 *     znaupd   ARPACK reverse communication interface routine.
 *     zneupd   ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     dznrm2   Level 1 BLAS that computes the norm of a complex vector.
 *     zaxpy    Level 1 BLAS that computes y <- alpha*x+y.
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
 * FILE: ndrv1.F   SID: 2.4   DATE OF SID: 10/17/00   RELEASE: 2
 *
 * \Remarks
 *     1. None
 *
 * \EndLib
 */
     /* --------------------------- */
     /* Define maximum dimensions   */
     /* for all arrays.             */
     /* MAXN:   Maximum dimension   */
     /*         of the A allowed.   */
     /* MAXNEV: Maximum NEV allowed */
     /* MAXNCV: Maximum NCV allowed */
     /* --------------------------- */

     /* --------------------- */
     /* Executable Statements */
     /* --------------------- */

     /* ------------------------------------------------ */
     /* The number NX is the number of interior points   */
     /* in the discretization of the 2-dimensional       */
     /* convection-diffusion operator on the unit        */
     /* square with zero Dirichlet boundary condition.   */
     /* The number N(=NX*NX) is the dimension of the     */
     /* matrix.  A standard eigenvalue problem is        */
     /* solved (BMAT = 'I').  NEV is the number of       */
     /* eigenvalues to be approximated.  The user can    */
     /* modify NX, NEV, NCV, WHICH to solve problems of  */
     /* different sizes, and to get different parts of   */
     /* the spectrum.  However, The following            */
     /* conditions must be satisfied:                    */
     /*                   N <= MAXN                      */
     /*                 NEV <= MAXNEV                    */
     /*           NEV + 2 <= NCV <= MAXNCV               */
     /* ------------------------------------------------ */

    nx = 10;
    n = nx * nx;
    nev = 4;
    ncv = 20;
    if (n > 256) {
	s_wsle(&io___5);
	do_lio(&c__9, &c__1, " ERROR with _NDRV1: N is greater than MAXN ", (ftnlen)43);
	e_wsle();
	goto L9000;
    } else if (nev > 12) {
	s_wsle(&io___6);
	do_lio(&c__9, &c__1, " ERROR with _NDRV1: NEV is greater than MAXNEV ", (ftnlen)47);
	e_wsle();
	goto L9000;
    } else if (ncv > 30) {
	s_wsle(&io___7);
	do_lio(&c__9, &c__1, " ERROR with _NDRV1: NCV is greater than MAXNCV ", (ftnlen)47);
	e_wsle();
	goto L9000;
    }
    *(unsigned char *)bmat = 'I';
    strcpy(which, "LM");

     /* ------------------------------------------------- */
     /* The work array WORKL is used in ZNAUPD  as         */
     /* workspace.  Its dimension LWORKL is set as        */
     /* illustrated below.  The parameter TOL determines  */
     /* the stopping criterion. If TOL<=0, machine        */
     /* precision is used.  The variable IDO is used for  */
     /* reverse communication, and is initially set to 0. */
     /* Setting INFO=0 indicates that a random vector is  */
     /* generated to start the ARNOLDI iteration.         */
     /* ------------------------------------------------- */

/* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 * 3 + ncv * 5;
    tol = 0.f;
    ido = 0;
    info = 0;

     /* ------------------------------------------------- */
     /* This program uses exact shift with respect to     */
     /* the current Hessenberg matrix (IPARAM(1) = 1).    */
     /* IPARAM(3) specifies the maximum number of Arnoldi */
     /* iterations allowed.  Mode 1 of ZNAUPD  is used     */
     /* (IPARAM(7) = 1). All these options can be changed */
     /* by the user. For details see the documentation in */
     /* ZNAUPD .                                           */
     /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode = 1;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

     /* ----------------------------------------- */
     /* M A I N   L O O P (Reverse communication) */
     /* ----------------------------------------- */

L10:

        /* ------------------------------------------- */
        /* Repeatedly call the routine ZNAUPD  and take */
        /* actions indicated by parameter IDO until    */
        /* either convergence is indicated or maxitr   */
        /* has been exceeded.                          */
        /* ------------------------------------------- */

    znaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, 
	    iparam, ipntr, workd, workl, &lworkl, rwork, &info, (ftnlen)1, (
	    ftnlen)2);

    if (ido == -1 || ido == 1) {

           /* ----------------------------------------- */
           /* Perform matrix vector multiplication      */
           /*                y <--- OP*x                */
           /* The user should supply his/her own        */
           /* matrix vector multiplication routine here */
           /* that takes workd(ipntr(1)) as the input   */
           /* vector, and return the matrix vector      */
           /* product to workd(ipntr(2)).               */
           /* ----------------------------------------- */

	zndrv1_av_(&nx, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

           /* --------------------------------------- */
           /* L O O P   B A C K to call ZNAUPD  again. */
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
        /* documentation in ZNAUPD   */
        /* ------------------------ */

	s_wsle(&io___24);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___25);
	do_lio(&c__9, &c__1, " Error with _naupd, info = ", (ftnlen)27);
	do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(int32_t));
	e_wsle();
	s_wsle(&io___26);
	do_lio(&c__9, &c__1, " Check the documentation of _naupd", (ftnlen)34);
	e_wsle();
	s_wsle(&io___27);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();

    } else {

        /* ----------------------------------------- */
        /* No fatal errors occurred.                 */
        /* Post-Process using ZNEUPD .                */
        /*                                           */
        /* Computed eigenvalues may be extracted.    */
        /*                                           */
        /* Eigenvectors may also be computed now if  */
        /* desired.  (indicated by rvec = .true.)    */
        /* ----------------------------------------- */

	rvec = true;

	zneupd_(&rvec, "A", select, d__, v, &c__256, &sigma, workev, bmat, &n,
		 which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, 
		workd, workl, &lworkl, rwork, &ierr, (ftnlen)1, (ftnlen)1, (
		ftnlen)2);

        /* -------------------------------------------- */
        /* Eigenvalues are returned in the one          */
        /* dimensional array D.  The corresponding      */
        /* eigenvectors are returned in the first NCONV */
        /* (=IPARAM(5)) columns of the two dimensional  */
        /* array V if requested.  Otherwise, an         */
        /* orthogonal basis for the invariant subspace  */
        /* corresponding to the eigenvalues in D is     */
        /* returned in V.                               */
        /* -------------------------------------------- */

	if (ierr != 0) {

           /* ---------------------------------- */
           /* Error condition:                   */
           /* Check the documentation of ZNEUPD . */
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
               /* eigenvectors.  (iparam(5) */
               /* indicates how many are    */
               /* accurate to the requested */
               /* tolerance)                */
               /* ------------------------- */

		zndrv1_av_(&nx, &v[(j << 8) - 256], ax);
		i__2 = j - 1;
		z__1.r = -d__[i__2].r, z__1.i = -d__[i__2].i;
		zaxpy_(&n, &z__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
		i__2 = j - 1;
		rd[j - 1] = d__[i__2].r;
		rd[j + 29] = d_imag(&d__[j - 1]);
		rd[j + 59] = dznrm2_(&n, ax, &c__1);
		rd[j + 59] /= dlapy2_(&rd[j - 1], &rd[j + 29]);
/* L20: */
	    }

            /* --------------------------- */
            /* Display computed residuals. */
            /* --------------------------- */

	    dmout_(&c__6, &nconv, &c__3, rd, &c__30, &c_n6, "Ritz values (Real, Imag) and relative residuals");
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
	do_lio(&c__9, &c__1, "_NDRV1", (ftnlen)6);
	e_wsle();
	s_wsle(&io___50);
	do_lio(&c__9, &c__1, "====== ", (ftnlen)7);
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
	do_lio(&c__5, &c__1, (char *)&tol, (ftnlen)sizeof(double));
	e_wsle();
	s_wsle(&io___60);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();

    }

     /* ------------------------- */
     /* Done with program zndrv1 . */
     /* ------------------------- */

L9000:

    return 0;
} /* MAIN__ */

/* ========================================================================== */

/*     matrix vector subroutine */

/*     The matrix used is the convection-diffusion operator */
/*     discretized using centered difference. */

int zndrv1_av_(int32_t *nx, zomplex *v, zomplex *w)
{
    /* System generated locals */
    int32_t i__1;
    zomplex z__1, z__2;

    /* Builtin functions */
    void z_div(zomplex *, zomplex *, zomplex *);

    /* Local variables */
    int32_t j;
    zomplex h2;
    int32_t lo;

/*     Computes w <--- OP*v, where OP is the nx*nx by nx*nx block */
/*     tridiagonal matrix */

/*                  | T -I          | */
/*                  |-I  T -I       | */
/*             OP = |   -I  T       | */
/*                  |        ...  -I| */
/*                  |           -I T| */

/*     derived from the standard central difference  discretization */
/*     of the convection-diffusion operator (Laplacian u) + rho*(du/dx) */
/*     with zero boundary condition. */

/*     The subroutine TV is called to computed y<---T*x. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    i__1 = (*nx + 1) * (*nx + 1);
    z__2.r = (double) i__1, z__2.i = 0.;
    z_div(&z__1, &c_b137, &z__2);
    h2.r = z__1.r, h2.i = z__1.i;

	zndrv1_tv_(nx, &v[1], &w[1]);
    z__2.r = -1., z__2.i = -0.;
    z_div(&z__1, &z__2, &h2);
    zaxpy_(nx, &z__1, &v[*nx + 1], &c__1, &w[1], &c__1);

    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j) {
	lo = (j - 1) * *nx;
	zndrv1_tv_(nx, &v[lo + 1], &w[lo + 1]);
	z__2.r = -1., z__2.i = -0.;
	z_div(&z__1, &z__2, &h2);
	zaxpy_(nx, &z__1, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);
	z__2.r = -1., z__2.i = -0.;
	z_div(&z__1, &z__2, &h2);
	zaxpy_(nx, &z__1, &v[lo + *nx + 1], &c__1, &w[lo + 1], &c__1);
/* L10: */
    }

    lo = (*nx - 1) * *nx;
	zndrv1_tv_(nx, &v[lo + 1], &w[lo + 1]);
    z__2.r = -1., z__2.i = -0.;
    z_div(&z__1, &z__2, &h2);
    zaxpy_(nx, &z__1, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);

    return 0;
} /* av_ */

/* ========================================================================= */
int zndrv1_tv_(int32_t *nx, zomplex *x, zomplex *y)
{
    /* System generated locals */
    int32_t i__1, i__2, i__3, i__4, i__5;
    zomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    void z_div(zomplex *, zomplex *, zomplex *);

    /* Local variables */
    zomplex h__;
    int32_t j;
    zomplex h2, dd, dl, du;

/*     Compute the matrix vector multiplication y<---T*x */
/*     where T is a nx by nx tridiagonal matrix with DD on the */
/*     diagonal, DL on the subdiagonal, and DU on the superdiagonal */

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    i__1 = *nx + 1;
    z__2.r = (double) i__1, z__2.i = 0.;
    z_div(&z__1, &c_b137, &z__2);
    h__.r = z__1.r, h__.i = z__1.i;
    z__1.r = h__.r * h__.r - h__.i * h__.i, z__1.i = h__.r * h__.i + h__.i * 
	    h__.r;
    h2.r = z__1.r, h2.i = z__1.i;
    z_div(&z__1, &c_b151_dx, &h2);
    dd.r = z__1.r, dd.i = z__1.i;
    z__3.r = -1., z__3.i = -0.;
    z_div(&z__2, &z__3, &h2);
    z__5.r = 50., z__5.i = 0.;
    z_div(&z__4, &z__5, &h__);
    z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - z__4.i;
    dl.r = z__1.r, dl.i = z__1.i;
    z__3.r = -1., z__3.i = -0.;
    z_div(&z__2, &z__3, &h2);
    z__5.r = 50., z__5.i = 0.;
    z_div(&z__4, &z__5, &h__);
    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
    du.r = z__1.r, du.i = z__1.i;

    z__2.r = dd.r * x[1].r - dd.i * x[1].i, z__2.i = dd.r * x[1].i + dd.i * x[
	    1].r;
    z__3.r = du.r * x[2].r - du.i * x[2].i, z__3.i = du.r * x[2].i + du.i * x[
	    2].r;
    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
    y[1].r = z__1.r, y[1].i = z__1.i;
    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j - 1;
	z__3.r = dl.r * x[i__3].r - dl.i * x[i__3].i, z__3.i = dl.r * x[i__3]
		.i + dl.i * x[i__3].r;
	i__4 = j;
	z__4.r = dd.r * x[i__4].r - dd.i * x[i__4].i, z__4.i = dd.r * x[i__4]
		.i + dd.i * x[i__4].r;
	z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
	i__5 = j + 1;
	z__5.r = du.r * x[i__5].r - du.i * x[i__5].i, z__5.i = du.r * x[i__5]
		.i + du.i * x[i__5].r;
	z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
	y[i__2].r = z__1.r, y[i__2].i = z__1.i;
/* L10: */
    }
    i__1 = *nx;
    i__2 = *nx - 1;
    z__2.r = dl.r * x[i__2].r - dl.i * x[i__2].i, z__2.i = dl.r * x[i__2].i + 
	    dl.i * x[i__2].r;
    i__3 = *nx;
    z__3.r = dd.r * x[i__3].r - dd.i * x[i__3].i, z__3.i = dd.r * x[i__3].i + 
	    dd.i * x[i__3].r;
    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
    y[i__1].r = z__1.r, y[i__1].i = z__1.i;
    return 0;
} /* tv_ */

