/* D:\Projekte\ARPACK\arpack-ng\SRC\dgetv0.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \BeginDoc
 *
 * \Name: dgetv0
 *
 * \Description:
 *  Generate a random initial residual vector for the Arnoldi process.
 *  Force the residual vector to be in the range of the operator OP.
 *
 * \Usage:
 *  call dgetv0
 *     ( IDO, BMAT, ITRY, INITV, N, J, V, LDV, RESID, RNORM,
 *       IPNTR, WORKD, IERR )
 *
 * \Arguments
 *  IDO     Integer.  (INPUT/OUTPUT)
 *          Reverse communication flag.  IDO must be zero on the first
 *          call to dgetv0.
 *          -------------------------------------------------------------
 *          IDO =  0: first call to the reverse communication interface
 *          IDO = -1: compute  Y = OP * X  where
 *                    IPNTR(1) is the pointer into WORKD for X,
 *                    IPNTR(2) is the pointer into WORKD for Y.
 *                    This is for the initialization phase to force the
 *                    starting vector into the range of OP.
 *          IDO =  2: compute  Y = B * X  where
 *                    IPNTR(1) is the pointer into WORKD for X,
 *                    IPNTR(2) is the pointer into WORKD for Y.
 *          IDO = 99: done
 *          -------------------------------------------------------------
 *
 *  BMAT    Character*1.  (INPUT)
 *          BMAT specifies the type of the matrix B in the (generalized)
 *          eigenvalue problem A*x = lambda*B*x.
 *          B = 'I' -> standard eigenvalue problem A*x = lambda*x
 *          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
 *
 *  ITRY    Integer.  (INPUT)
 *          ITRY counts the number of times that dgetv0 is called.
 *          It should be set to 1 on the initial call to dgetv0.
 *
 *  INITV   Logical variable.  (INPUT)
 *          .TRUE.  => the initial residual vector is given in RESID.
 *          .FALSE. => generate a random initial residual vector.
 *
 *  N       Integer.  (INPUT)
 *          Dimension of the problem.
 *
 *  J       Integer.  (INPUT)
 *          Index of the residual vector to be generated, with respect to
 *          the Arnoldi process.  J > 1 in case of a "restart".
 *
 *  V       Double precision N by J array.  (INPUT)
 *          The first J-1 columns of V contain the current Arnoldi basis
 *          if this is a "restart".
 *
 *  LDV     Integer.  (INPUT)
 *          Leading dimension of V exactly as declared in the calling
 *          program.
 *
 *  RESID   Double precision array of length N.  (INPUT/OUTPUT)
 *          Initial residual vector to be generated.  If RESID is
 *          provided, force RESID into the range of the operator OP.
 *
 *  RNORM   Double precision scalar.  (OUTPUT)
 *          B-norm of the generated residual.
 *
 *  IPNTR   Integer array of length 3.  (OUTPUT)
 *
 *  WORKD   Double precision work array of length 2*N.  (REVERSE COMMUNICATION).
 *          On exit, WORK(1:N) = B*RESID to be used in SSAITR.
 *
 *  IERR    Integer.  (OUTPUT)
 *          =  0: Normal exit.
 *          = -1: Cannot generate a nontrivial restarted residual vector
 *                in the range of the operator OP.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Local variables:
 *     xxxxxx  float
 *
 * \References:
 *  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
 *     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
 *     pp 357-385.
 *  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
 *     Restarted Arnoldi Iteration", Rice University Technical Report
 *     TR95-13, Department of Computational and Applied Mathematics.
 *
 * \Routines called:
 *     arscnd  ARPACK utility routine for timing.
 *     dvout   ARPACK utility routine for vector output.
 *     dlarnv  LAPACK routine for generating a random vector.
 *     dgemv   Level 2 BLAS routine for matrix vector multiplication.
 *     dcopy   Level 1 BLAS that copies one vector to another.
 *     ddot    Level 1 BLAS that computes the scalar product of two vectors.
 *     dnrm2   Level 1 BLAS that computes the norm of a vector.
 *
 * \Author
 *     Danny Sorensen               Phuong Vu
 *     Richard Lehoucq              CRPC / Rice University
 *     Dept. of Computational &     Houston, Texas
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * \SCCS Information: @(#)
 * FILE: getv0.F   SID: 2.7   DATE OF SID: 04/07/99   RELEASE: 2
 *
 * \EndLib
 */

int dgetv0_(int32_t *ido, char *bmat, int32_t *itry, bool 
	*initv, int32_t *n, int32_t *j, double *v, int32_t *ldv, 
	double *resid, double *rnorm, int32_t *ipntr, double *
	workd, int32_t *ierr)
{
    /* Initialized data */

    static bool inits = true;

    /* System generated locals */
    int32_t v_dim1, v_offset, i__1;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static float t0, t1, t2, t3;
    int32_t jj;
    static int32_t iter;
    static bool orth;
    static int32_t iseed[4];
    int32_t idist;
    static bool first;
    static double rnorm0;
    static int32_t msglvl;

/*     %-----------------% */
/*     | Data Statements | */
/*     %-----------------% */

    /* Parameter adjustments */
    --workd;
    --resid;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --ipntr;

    /* Function Body */

/*     %-----------------------% */
/*     | Executable Statements | */
/*     %-----------------------% */

/*     %-----------------------------------% */
/*     | Initialize the seed of the LAPACK | */
/*     | random number generator           | */
/*     %-----------------------------------% */

    if (inits) {
	iseed[0] = 1;
	iseed[1] = 3;
	iseed[2] = 5;
	iseed[3] = 7;
	inits = false;
    }

    if (*ido == 0) {

/*        %-------------------------------% */
/*        | Initialize timing statistics  | */
/*        | & message level for debugging | */
/*        %-------------------------------% */

	arscnd_(&t0);
	msglvl = debug_1.mgetv0;

	*ierr = 0;
	iter = 0;
	first = false;
	orth = false;

/*        %-----------------------------------------------------% */
/*        | Possibly generate a random starting vector in RESID | */
/*        | Use a LAPACK random number generator used by the    | */
/*        | matrix generation routines.                         | */
/*        |    idist = 1: uniform (0,1)  distribution;          | */
/*        |    idist = 2: uniform (-1,1) distribution;          | */
/*        |    idist = 3: normal  (0,1)  distribution;          | */
/*        %-----------------------------------------------------% */

	if (! (*initv)) {
	    idist = 2;
	    dlarnv_(&idist, iseed, n, &resid[1]);
	}

/*        %----------------------------------------------------------% */
/*        | Force the starting vector into the range of OP to handle | */
/*        | the generalized problem when B is possibly (singular).   | */
/*        %----------------------------------------------------------% */

	arscnd_(&t2);
	if (*itry == 1) {
	    ++timing_1.nopx;
	    ipntr[1] = 1;
	    ipntr[2] = *n + 1;
	    dcopy_(n, &resid[1], &c__1, &workd[1], &c__1);
	    *ido = -1;
	    goto L9000;
	} else if (*itry > 1 && *(unsigned char *)bmat == 'G') {
	    dcopy_(n, &resid[1], &c__1, &workd[*n + 1], &c__1);
	}
    }

/*     %-----------------------------------------% */
/*     | Back from computing OP*(initial-vector) | */
/*     %-----------------------------------------% */

    if (first) {
	goto L20;
    }

/*     %-----------------------------------------------% */
/*     | Back from computing OP*(orthogonalized-vector) | */
/*     %-----------------------------------------------% */

    if (orth) {
	goto L40;
    }

    if (*(unsigned char *)bmat == 'G') {
	arscnd_(&t3);
	timing_1.tmvopx += t3 - t2;
    }

/*     %------------------------------------------------------% */
/*     | Starting vector is now in the range of OP; r = OP*r; | */
/*     | Compute B-norm of starting vector.                   | */
/*     %------------------------------------------------------% */

    arscnd_(&t2);
    first = true;
    if (*itry == 1) {
	dcopy_(n, &workd[*n + 1], &c__1, &resid[1], &c__1);
    }
    if (*(unsigned char *)bmat == 'G') {
	++timing_1.nbx;
	ipntr[1] = *n + 1;
	ipntr[2] = 1;
	*ido = 2;
	goto L9000;
    } else if (*(unsigned char *)bmat == 'I') {
	dcopy_(n, &resid[1], &c__1, &workd[1], &c__1);
    }

L20:

    if (*(unsigned char *)bmat == 'G') {
	arscnd_(&t3);
	timing_1.tmvbx += t3 - t2;
    }

    first = false;
    if (*(unsigned char *)bmat == 'G') {
	rnorm0 = ddot_(n, &resid[1], &c__1, &workd[1], &c__1);
	rnorm0 = sqrt((abs(rnorm0)));
    } else if (*(unsigned char *)bmat == 'I') {
	rnorm0 = dnrm2_(n, &resid[1], &c__1);
    }
    *rnorm = rnorm0;

/*     %---------------------------------------------% */
/*     | Exit if this is the very first Arnoldi step | */
/*     %---------------------------------------------% */

    if (*j == 1) {
	goto L50;
    }

/*     %---------------------------------------------------------------- */
/*     | Otherwise need to B-orthogonalize the starting vector against | */
/*     | the current Arnoldi basis using Gram-Schmidt with iter. ref.  | */
/*     | This is the case where an invariant subspace is encountered   | */
/*     | in the middle of the Arnoldi factorization.                   | */
/*     |                                                               | */
/*     |       s = V^{T}*B*r;   r = r - V*s;                           | */
/*     |                                                               | */
/*     | Stopping criteria used for iter. ref. is discussed in         | */
/*     | Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   | */
/*     %---------------------------------------------------------------% */

    orth = true;
L30:

    i__1 = *j - 1;
    dgemv_("T", n, &i__1, &d_one, &v[v_offset], ldv, &workd[1], &c__1, &d_zero,
	     &workd[*n + 1], &c__1);
    i__1 = *j - 1;
    dgemv_("N", n, &i__1, &d_m1, &v[v_offset], ldv, &workd[*n + 1], &c__1, &
	    d_one, &resid[1], &c__1);

/*     %----------------------------------------------------------% */
/*     | Compute the B-norm of the orthogonalized starting vector | */
/*     %----------------------------------------------------------% */

    arscnd_(&t2);
    if (*(unsigned char *)bmat == 'G') {
	++timing_1.nbx;
	dcopy_(n, &resid[1], &c__1, &workd[*n + 1], &c__1);
	ipntr[1] = *n + 1;
	ipntr[2] = 1;
	*ido = 2;
	goto L9000;
    } else if (*(unsigned char *)bmat == 'I') {
	dcopy_(n, &resid[1], &c__1, &workd[1], &c__1);
    }

L40:

    if (*(unsigned char *)bmat == 'G') {
	arscnd_(&t3);
	timing_1.tmvbx += t3 - t2;
    }

    if (*(unsigned char *)bmat == 'G') {
	*rnorm = ddot_(n, &resid[1], &c__1, &workd[1], &c__1);
	*rnorm = sqrt((abs(*rnorm)));
    } else if (*(unsigned char *)bmat == 'I') {
	*rnorm = dnrm2_(n, &resid[1], &c__1);
    }

/*     %--------------------------------------% */
/*     | Check for further orthogonalization. | */
/*     %--------------------------------------% */

    if (msglvl > 2) {
	dvout_(&debug_1.logfil, &c__1, &rnorm0, &debug_1.ndigit, "_getv0: re"
		"-orthonalization ; rnorm0 is", (ftnlen)38);
	dvout_(&debug_1.logfil, &c__1, rnorm, &debug_1.ndigit, "_getv0: re-o"
		"rthonalization ; rnorm is", (ftnlen)37);
    }

    if (*rnorm > rnorm0 * .717f) {
	goto L50;
    }

    ++iter;
    if (iter <= 5) {

/*        %-----------------------------------% */
/*        | Perform iterative refinement step | */
/*        %-----------------------------------% */

	rnorm0 = *rnorm;
	goto L30;
    } else {

/*        %------------------------------------% */
/*        | Iterative refinement step "failed" | */
/*        %------------------------------------% */

	i__1 = *n;
	for (jj = 1; jj <= i__1; ++jj) {
	    resid[jj] = 0.;
/* L45: */
	}
	*rnorm = 0.;
	*ierr = -1;
    }

L50:

    if (msglvl > 0) {
	dvout_(&debug_1.logfil, &c__1, rnorm, &debug_1.ndigit, "_getv0: B-no"
		"rm of initial / restarted starting vector", (ftnlen)53);
    }
    if (msglvl > 3) {
	dvout_(&debug_1.logfil, n, &resid[1], &debug_1.ndigit, "_getv0: init"
		"ial / restarted starting vector", (ftnlen)43);
    }
    *ido = 99;

    arscnd_(&t1);
    timing_1.tgetv0 += t1 - t0;

L9000:
    return 0;

/*     %---------------% */
/*     | End of dgetv0 | */
/*     %---------------% */

} /* dgetv0_ */

