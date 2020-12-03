/* EXAMPLES\NONSYM\dndrv6.f -- translated by f2c (version 20100827). */

#include "arpack.h"

int dndrv6()
{
    /* System generated locals */
    int32_t i__1, i__2, i__3;
    double d__1, d__2;
    zomplex z__1;

    double d_imag(zomplex *);

    /* Local variables */
    double d[75]	/* was [25][3] */;
    int32_t j, n;
    double v[6400]	/* was [256][25] */;
    zomplex c1, c2, c3;
    double ax[256];
    double mx[256];
    zomplex cdd[256], cdl[256], cdu[256];
    int32_t ido, ncv, nev;
    double tol;
    zomplex cdu2[256];
    double deni;
    char* bmat;
    int32_t mode;
    double denr;
    int32_t info;
    bool rvec;
    int32_t ierr, ipiv[256];
    double numi, numr;
    char* which;
    double resid[256];
    zomplex ctemp[256];
    int32_t nconv;
    double workd[768];
    bool first;
    int32_t ipntr[14];
    double workl[2025];
    int32_t iparam[11];
    double sigmai;
    bool select[25];
    double sigmar;
    int32_t ishfts, maxitr, lworkl;
    double workev[75];

    /* Fortran I/O blocks */

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
/*     ... Use mode 4 of DNAUPD. */
/**
 * \BeginLib
 *
 * \Routines called:
 *     dnaupd  ARPACK reverse communication interface routine.
 *     dneupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     zgttrf  LAPACK complex matrix factorization routine.
 *     zgttrs  LAPACK complex linear system solve routine.
 *     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
 *     ddot    Level 1 BLAS that computes the dot product of two vectors.
 *     dnrm2   Level 1 BLAS that computes the norm of a vector.
 *     av      Matrix vector subroutine that computes A*x.
 *     mv      Matrix vector subroutine that computes M*x.
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
 * FILE: ndrv6.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
 *
 * \Remarks
 *     1. None
 *
 * \EndLib
 */
     /* --------------------------- */
     /* Define leading dimensions   */
     /* for all arrays.             */
     /* MAXN:   Maximum dimension   */
     /*         of the A allowed.   */
     /* MAXNEV: Maximum NEV allowed */
     /* MAXNCV: Maximum NCV allowed */
     /* --------------------------- */

     /* --------------------- */
     /* Executable Statements */
     /* --------------------- */

     /* -------------------------------------------------- */
     /* The number N is the dimension of the matrix.  A    */
     /* generalized eigenvalue problem is solved (BMAT =   */
     /* 'G').  NEV is the number of eigenvalues (closest   */
     /* to the shift (SIGMAR,SIGMAI)) to be approximated.  */
     /* Since the shift-invert mode is used, WHICH is set  */
     /* to 'LM'.  The user can modify NEV, NCV, SIGMA to   */
     /* solve problems of different sizes, and to get      */
     /* different parts of the spectrum. However, The      */
     /* following conditions must be satisfied:            */
     /*                     N <= MAXN,                     */
     /*                   NEV <= MAXNEV,                   */
     /*               NEV + 2 <= NCV <= MAXNCV             */
     /* -------------------------------------------------- */

    n = 100;
    nev = 4;
    ncv = 20;
    if (n > 256) {
	printf(" ERROR with _NDRV6: N is greater than MAXN \n");
	return 0;
    } else if (nev > 10) {
	printf(" ERROR with _NDRV6: NEV is greater than MAXNEV \n");
	return 0;
    } else if (ncv > 25) {
	printf(" ERROR with _NDRV6: NCV is greater than MAXNCV \n");
	return 0;
    }
    bmat = "G";
    strcpy(which, "LM");
    sigmar = .4;
    sigmai = .6;

     /* -------------------------------------------------- */
     /* Construct C = A - (SIGMAR,SIGMAI)*M in complex     */
     /* arithmetic, and factor C in complex arithmetic     */
     /* (using LAPACK subroutine zgttrf). The matrix A is  */
     /* chosen to be the tridiagonal matrix with -2 on the */
     /* subdiagonal, 2 on the diagonal and 3 on the        */
     /* superdiagonal. The matrix M is chosen to be the    */
     /* symmetric tridiagonal matrix with 4 on the         */
     /* diagonal and 1 on the off-diagonals.               */
     /* -------------------------------------------------- */

    d__1 = -2. - sigmar;
    d__2 = -sigmai;
    z__1.r = d__1, z__1.i = d__2;
    c1.r = z__1.r, c1.i = z__1.i;
    d__1 = 2. - sigmar * 4.;
    d__2 = sigmai * -4.;
    z__1.r = d__1, z__1.i = d__2;
    c2.r = z__1.r, c2.i = z__1.i;
    d__1 = 3. - sigmar;
    d__2 = -sigmai;
    z__1.r = d__1, z__1.i = d__2;
    c3.r = z__1.r, c3.i = z__1.i;

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

    zgttrf_(&n, cdl, cdd, cdu, cdu2, ipiv, &ierr);
    if (ierr != 0) {
	printf(" \n");
	printf(" ERROR with _gttrf in _NDRV6.\n");
	printf(" \n");
	return 0;
    }

     /* --------------------------------------------------- */
     /* The work array WORKL is used in DNAUPD as           */
     /* workspace.  Its dimension LWORKL is set as          */
     /* illustrated below.  The parameter TOL determines    */
     /* the stopping criterion. If TOL<=0, machine          */
     /* precision is used.  The variable IDO is used for    */
     /* reverse communication, and is initially set to 0.   */
     /* Setting INFO=0 indicates that a random vector is    */
     /* generated in DNAUPD to start the Arnoldi iteration. */
     /* --------------------------------------------------- */

/* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 * 3 + ncv * 6;
    tol = 0.;
    ido = 0;
    info = 0;

     /* ------------------------------------------------- */
     /* This program uses exact shift with respect to     */
     /* the current Hessenberg matrix (IPARAM(1) = 1).    */
     /* IPARAM(3) specifies the maximum number of Arnoldi */
     /* iterations allowed.  Mode 3 of DNAUPD is used     */
     /* (IPARAM(7) = 3).  All these options can be        */
     /* changed by the user. For details, see the         */
     /* documentation in DNAUPD.                          */
     /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode = 3;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

     /* ---------------------------------------- */
     /* M A I N   L O O P(Reverse communication) */
     /* ---------------------------------------- */

L20:

        /* ------------------------------------------- */
        /* Repeatedly call the routine DNAUPD and take */
        /* actions indicated by parameter IDO until    */
        /* either convergence is indicated or maxitr   */
        /* has been exceeded.                          */
        /* ------------------------------------------- */

    dnaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, 
	    iparam, ipntr, workd, workl, &lworkl, &info, (ftnlen)1, (ftnlen)2)
	    ;

    if (ido == -1) {

           /* ---------------------------------------------------------- */
           /*                           Perform                          */
           /* y <--- OP*x = Imaginary_Part{inv[A-(SIGMAR,SIGMAI)*M]*M*x} */
           /* to force starting vector into the range of OP. The user    */
           /* should supply his/her own matrix vector multiplication     */
           /* routine and a complex linear system solver.  The matrix    */
           /* vector multiplication routine should take workd(ipntr(1))  */
           /* as the input. The final result (a real vector) should be   */
           /* returned to workd(ipntr(2)).                               */
           /* ---------------------------------------------------------- */

	dndrv6_mv_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = j - 1;
	    i__3 = ipntr[1] + j - 2;
	    z__1.r = workd[i__3], z__1.i = 0.;
	    ctemp[i__2].r = z__1.r, ctemp[i__2].i = z__1.i;
/* L30: */
	}

	zgttrs_("N", &n, &c__1, cdl, cdd, cdu, cdu2, ipiv, ctemp, &c__256, &
		ierr, (ftnlen)1);
	if (ierr != 0) {
	    printf(" \n");
	    printf(" ERROR with _gttrs in _NDRV6.\n");
	    printf(" \n");
	    return ierr;
	}
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
	    workd[ipntr[1] + j - 2] = d_imag(&ctemp[j - 1]);
/* L40: */
	}

           /* --------------------------------------- */
           /* L O O P   B A C K to call DNAUPD again. */
           /* --------------------------------------- */

	goto L20;

    } else if (ido == 1) {

           /* ---------------------------------------------------------- */
           /*                          Perform                           */
           /* y <--- OP*x = Imaginary_Part{inv[A-(SIGMAR,SIGMAI)*M]*M*x} */
           /* M*x has been saved in workd(ipntr(3)). The user only need  */
           /* the complex linear system solver here that takes           */
           /* complex[workd(ipntr(3))] as input, and returns the result  */
           /* to workd(ipntr(2)).                                        */
           /* ---------------------------------------------------------- */

	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = j - 1;
	    i__3 = ipntr[2] + j - 2;
	    z__1.r = workd[i__3], z__1.i = 0.;
	    ctemp[i__2].r = z__1.r, ctemp[i__2].i = z__1.i;
/* L50: */
	}
	zgttrs_("N", &n, &c__1, cdl, cdd, cdu, cdu2, ipiv, ctemp, &c__256, &
		ierr, (ftnlen)1);
	if (ierr != 0) {
	    printf(" \n");
	    printf(" ERROR with _gttrs in _NDRV6.\n");
	    printf(" \n");
	    return ierr;
	}
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
	    workd[ipntr[1] + j - 2] = d_imag(&ctemp[j - 1]);
/* L60: */
	}

           /* --------------------------------------- */
           /* L O O P   B A C K to call DNAUPD again. */
           /* --------------------------------------- */

	goto L20;

    } else if (ido == 2) {

           /* ------------------------------------------- */
           /*          Perform  y <--- M*x                */
           /* Need matrix vector multiplication routine   */
           /* here that takes workd(ipntr(1)) as input    */
           /* and returns the result to workd(ipntr(2)).  */
           /* ------------------------------------------- */

	dndrv6_mv_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

           /* --------------------------------------- */
           /* L O O P   B A C K to call DNAUPD again. */
           /* --------------------------------------- */

	goto L20;

    }

     /* --------------------------------------- */
     /* Either we have convergence, or there is */
     /* an error.                               */
     /* --------------------------------------- */

    if (info < 0) {

        /* ------------------------ */
        /* Error message, check the */
        /* documentation in DNAUPD. */
        /* ------------------------ */

	printf(" \n");
	printf(" Error with _naupd info = %d\n", info);
	printf(" Check the documentation of _naupd.\n");
	printf(" \n");

    } else {

        /* ----------------------------------------- */
        /* No fatal errors occurred.                 */
        /* Post-Process using DNEUPD.                */
        /*                                           */
        /* Computed eigenvalues may be extracted.    */
        /*                                           */
        /* Eigenvectors may also be computed now if  */
        /* desired.  (indicated by rvec = .true.)    */
        /* ----------------------------------------- */

	rvec = true;
	dneupd_(&rvec, "A", select, d, &d[25], v, &c__256, &sigmar, &
		sigmai, workev, bmat, &n, which, &nev, &tol, resid, &ncv, v, &
		c__256, iparam, ipntr, workd, workl, &lworkl, &ierr, (ftnlen)
		1, (ftnlen)1, (ftnlen)2);

        /* --------------------------------------------- */
        /* The real part of the eigenvalue is returned   */
        /* in the first column of the two dimensional    */
        /* array D, and the IMAGINARY part is returned   */
        /* in the second column of D.  The corresponding */
        /* eigenvectors are returned in the first NEV    */
        /* columns of the two dimensional array V if     */
        /* requested.  Otherwise, an orthogonal basis    */
        /* for the invariant subspace corresponding to   */
        /* the eigenvalues in D is returned in V.        */
        /* --------------------------------------------- */

	if (ierr != 0) {

           /* ---------------------------------- */
           /* Error condition:                   */
           /* Check the documentation of DNEUPD. */
           /* ---------------------------------- */

	    printf(" \n");
	    printf(" Error with _neupd info = %d\n", ierr);
	    printf(" Check the documentation of _neupd. \n");
	    printf(" \n");

	} else {

	    first = true;
	    nconv = iparam[4];
	    i__1 = nconv;
	    for (j = 1; j <= i__1; ++j) {

               /* ----------------------------------- */
               /* Use Rayleigh Quotient to recover    */
               /* eigenvalues of the original problem.*/
               /* ----------------------------------- */

		if (d[j + 24] == 0.) {

                  /* -------------------------- */
                  /* Eigenvalue is real.        */
                  /* Compute d = x'(Ax)/x'(Mx). */
                  /* -------------------------- */

			dndrv6_av_(&n, &v[(j << 8) - 256], ax);
		    numr = ddot_(&n, &v[(j << 8) - 256], &c__1, ax, &c__1);
			dndrv6_mv_(&n, &v[(j << 8) - 256], ax);
		    denr = ddot_(&n, &v[(j << 8) - 256], &c__1, ax, &c__1);
		    d[j - 1] = numr / denr;

		} else if (first) {

                  /* ---------------------- */
                  /* Eigenvalue is complex. */
                  /* Compute the first one  */
                  /* of the conjugate pair. */
                  /* ---------------------- */

                  /* -------------- */
                  /* Compute x'(Ax) */
                  /* -------------- */

			dndrv6_av_(&n, &v[(j << 8) - 256], ax);
		    numr = ddot_(&n, &v[(j << 8) - 256], &c__1, ax, &c__1);
		    numi = ddot_(&n, &v[(j + 1 << 8) - 256], &c__1, ax, &c__1)
			    ;
			dndrv6_av_(&n, &v[(j + 1 << 8) - 256], ax);
		    numr += ddot_(&n, &v[(j + 1 << 8) - 256], &c__1, ax, &
			    c__1);
		    numi = -numi + ddot_(&n, &v[(j << 8) - 256], &c__1, ax, &
			    c__1);

                  /* -------------- */
                  /* Compute x'(Mx) */
                  /* -------------- */

			dndrv6_mv_(&n, &v[(j << 8) - 256], ax);
		    denr = ddot_(&n, &v[(j << 8) - 256], &c__1, ax, &c__1);
		    deni = ddot_(&n, &v[(j + 1 << 8) - 256], &c__1, ax, &c__1)
			    ;
			dndrv6_mv_(&n, &v[(j + 1 << 8) - 256], ax);
		    denr += ddot_(&n, &v[(j + 1 << 8) - 256], &c__1, ax, &
			    c__1);
		    deni = -deni + ddot_(&n, &v[(j << 8) - 256], &c__1, ax, &
			    c__1);

                  /* -------------- */
                  /* d=x'(Ax)/x'(Mx)*/
                  /* -------------- */

		    d[j - 1] = (numr * denr + numi * deni) / dlapy2_(&denr, 
			    &deni);
		    d[j + 24] = (numi * denr - numr * deni) / dlapy2_(&denr,
			     &deni);
		    first = false;

		} else {

                  /* ---------------------------- */
                  /* Get the second eigenvalue of */
                  /* the conjugate pair by taking */
                  /* the conjugate of the last    */
                  /* eigenvalue computed.         */
                  /* ---------------------------- */

		    d[j - 1] = d[j - 2];
		    d[j + 24] = -d[j + 23];
		    first = true;

		}

/* L70: */
	    }

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

	    first = true;
	    nconv = iparam[4];
	    i__1 = nconv;
	    for (j = 1; j <= i__1; ++j) {

		if (d[j + 24] == 0.) {

                 /* ------------------ */
                 /* Ritz value is real */
                 /* ------------------ */

			dndrv6_av_(&n, &v[(j << 8) - 256], ax);
			dndrv6_mv_(&n, &v[(j << 8) - 256], mx);
		    d__1 = -d[j - 1];
		    daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
		    d[j + 49] = dnrm2_(&n, ax, &c__1);
		    d[j + 49] /= (d__1 = d[j - 1], abs(d__1));

		} else if (first) {

                 /* ---------------------- */
                 /* Ritz value is complex  */
                 /* Residual of one Ritz   */
                 /* value of the conjugate */
                 /* pair is computed.      */
                 /* ---------------------- */

			dndrv6_av_(&n, &v[(j << 8) - 256], ax);
			dndrv6_mv_(&n, &v[(j << 8) - 256], mx);
		    d__1 = -d[j - 1];
		    daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
			dndrv6_mv_(&n, &v[(j + 1 << 8) - 256], mx);
		    daxpy_(&n, &d[j + 24], mx, &c__1, ax, &c__1);
		    d[j + 49] = dnrm2_(&n, ax, &c__1);
			dndrv6_av_(&n, &v[(j + 1 << 8) - 256], ax);
			dndrv6_mv_(&n, &v[(j + 1 << 8) - 256], mx);
		    d__1 = -d[j - 1];
		    daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
			dndrv6_mv_(&n, &v[(j << 8) - 256], mx);
		    d__1 = -d[j + 24];
		    daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
		    d__1 = dnrm2_(&n, ax, &c__1);
		    d[j + 49] = dlapy2_(&d[j + 49], &d__1);
		    d[j + 49] /= dlapy2_(&d[j - 1], &d[j + 24]);
		    d[j + 50] = d[j + 49];
		    first = false;
		} else {
		    first = true;
		}

/* L80: */
	    }

           /* --------------------------- */
           /* Display computed residuals. */
           /* --------------------------- */

	    dmout_(&nconv, &c__3, d, &c__25, &c_n6, "Ritz values (Real,Imag) and relative residuals");

	}

       /* ----------------------------------------- */
       /* Print additional convergence information. */
       /* ----------------------------------------- */

	if (info == 1) {
	    printf(" \n");
	    printf(" Maximum number of iterations reached.\n");
	    printf(" \n");
	} else if (info == 3) {
	    printf(" \n");
	    printf(" No shifts could be applied during implicit\n");
	    printf(" Arnoldi update try increasing NCV.\n");
	    printf(" \n");
	}

	printf(" \n");
	printf(" _NDRV6 \n");
	printf(" ====== \n");
	printf(" \n");
	printf(" Size of the matrix is %d\n", n);
	printf(" The number of Ritz values requested is %d\n", nev);
	printf(" The number of Arnoldi vectors generated (NCV) is %d\n", ncv);
	printf(" What portion of the spectrum: %s\n", which);
	printf(" The number of converged Ritz values is %d\n", nconv);
	printf(" The number of Implicit Arnoldi update iterations taken is %d\n", iparam[2]);
	printf(" The number of OP*x is %d\n", iparam[8]);
	printf(" The convergence criterion is %f\n", tol);
	printf(" \n");

    }

     /* ------------------------- */
     /* Done with program dndrv6. */
     /* ------------------------- */

    return 0;
} /* MAIN__ */

/* ========================================================================== */

/*     matrix vector multiplication subroutine */

int dndrv6_mv_(int32_t *n, double *v, double *w)
{
    /* System generated locals */
    int32_t i__1;

    /* Local variables */
    int32_t j;

/*     Compute the matrix vector multiplication y<---M*x */
/*     where M is a n by n symmetric tridiagonal matrix with 4 on the */
/*     diagonal, 1 on the subdiagonal and superdiagonal. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    w[1] = v[1] * 4. + v[2] * 1.;
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j) {
	w[j] = v[j - 1] * 1. + v[j] * 4. + v[j + 1] * 1.;
/* L10: */
    }
    w[*n] = v[*n - 1] * 1. + v[*n] * 4.;
    return 0;
} /* mv_ */

/* ------------------------------------------------------------------ */
int dndrv6_av_(int32_t *n, double *v, double *w)
{
    /* System generated locals */
    int32_t i__1;

    /* Local variables */
    int32_t j;

/*     Compute the matrix vector multiplication y<---A*x */
/*     where M is a n by n symmetric tridiagonal matrix with 2 on the */
/*     diagonal, -2 on the subdiagonal and 3 on the superdiagonal. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    w[1] = v[1] * 2. + v[2] * 3.;
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j) {
	w[j] = v[j - 1] * -2. + v[j] * 2. + v[j + 1] * 3.;
/* L10: */
    }
    w[*n] = v[*n - 1] * -2. + v[*n] * 2.;
    return 0;
} /* av_ */

