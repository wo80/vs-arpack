/* arpack-ng\SRC\sneigh.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \BeginDoc
 *
 * \Name: sneigh
 *
 * \Description:
 *  Compute the eigenvalues of the current upper Hessenberg matrix
 *  and the corresponding Ritz estimates given the current residual norm.
 *
 * \Usage:
 *  call sneigh
 *     ( RNORM, N, H, LDH, RITZR, RITZI, BOUNDS, Q, LDQ, WORKL, IERR )
 *
 * \Arguments
 *  RNORM   Real scalar.  (INPUT)
 *          Residual norm corresponding to the current upper Hessenberg
 *          matrix H.
 *
 *  N       Integer.  (INPUT)
 *          Size of the matrix H.
 *
 *  H       Real N by N array.  (INPUT)
 *          H contains the current upper Hessenberg matrix.
 *
 *  LDH     Integer.  (INPUT)
 *          Leading dimension of H exactly as declared in the calling
 *          program.
 *
 *  RITZR,  Real arrays of length N.  (OUTPUT)
 *  RITZI   On output, RITZR(1:N) (resp. RITZI(1:N)) contains the float
 *          (respectively imaginary) parts of the eigenvalues of H.
 *
 *  BOUNDS  Real array of length N.  (OUTPUT)
 *          On output, BOUNDS contains the Ritz estimates associated with
 *          the eigenvalues RITZR and RITZI.  This is equal to RNORM
 *          times the last components of the eigenvectors corresponding
 *          to the eigenvalues in RITZR and RITZI.
 *
 *  Q       Real N by N array.  (WORKSPACE)
 *          Workspace needed to store the eigenvectors of H.
 *
 *  LDQ     Integer.  (INPUT)
 *          Leading dimension of Q exactly as declared in the calling
 *          program.
 *
 *  WORKL   Real work array of length N**2 + 3*N.  (WORKSPACE)
 *          Private (replicated) array on each PE or array allocated on
 *          the front end.  This is needed to keep the full Schur form
 *          of H and also in the calculation of the eigenvectors of H.
 *
 *  IERR    Integer.  (OUTPUT)
 *          Error exit flag from slahqr or strevc.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Local variables:
 *     xxxxxx  float
 *
 * \Routines called:
 *     slahqr  LAPACK routine to compute the float Schur form of an
 *             upper Hessenberg matrix and last row of the Schur vectors.
 *     arscnd  ARPACK utility routine for timing.
 *     smout   ARPACK utility routine that prints matrices
 *     svout   ARPACK utility routine that prints vectors.
 *     slacpy  LAPACK matrix copy routine.
 *     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     strevc  LAPACK routine to compute the eigenvectors of a matrix
 *             in upper quasi-triangular form
 *     sgemv   Level 2 BLAS routine for matrix vector multiplication.
 *     scopy   Level 1 BLAS that copies one vector to another .
 *     snrm2   Level 1 BLAS that computes the norm of a vector.
 *     sscal   Level 1 BLAS that scales a vector.
 *
 *
 * \Author
 *     Danny Sorensen               Phuong Vu
 *     Richard Lehoucq              CRPC / Rice University
 *     Dept. of Computational &     Houston, Texas
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * \Revision history:
 *     xx/xx/92: Version ' 2.1'
 *
 * \SCCS Information: @(#)
 * FILE: neigh.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
 *
 * \Remarks
 *     None
 *
 * \EndLib
 */

int sneigh_(float *rnorm, int32_t *n, float *h, int32_t *ldh,
	 float *ritzr, float *ritzi, float *bounds, float *q, int32_t *ldq, float *
	workl, int32_t *ierr)
{
    /* System generated locals */
    int32_t h_dim1, h_offset, q_dim1, q_offset, i__1;
    float r__1, r__2;

    /* Local variables */
    int32_t i, j;
    static float t0, t1;
    float vl[1], temp;
    int32_t iconj;
    bool select[1];
    int32_t msglvl;

     /* --------------------- */
     /* Executable Statements */
     /* --------------------- */

     /* ----------------------------- */
     /* Initialize timing statistics  */
     /* & message level for debugging */
     /* ----------------------------- */

    /* Parameter adjustments */
    --workl;
    --bounds;
    --ritzi;
    --ritzr;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h -= h_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;

    /* Function Body */
#ifndef NO_TIMER
    arscnd_(&t0);
#endif

    msglvl = debug_1.mneigh;

#ifndef NO_TRACE
    if (msglvl > 2) {
	smout_(n, n, &h[h_offset], ldh, &debug_1.ndigit, "_neigh: Entering upper Hessenberg matrix H ");
    }
#endif

     /* --------------------------------------------------------- */
     /* 1. Compute the eigenvalues, the last components of the    */
     /*    corresponding Schur vectors and the full Schur form T  */
     /*    of the current upper Hessenberg matrix H.              */
     /* slahqr returns the full Schur form of H in WORKL(1:N**2)  */
     /* and the last components of the Schur vectors in BOUNDS.   */
     /* --------------------------------------------------------- */

    slacpy_("A", n, n, &h[h_offset], ldh, &workl[1], n);
    i__1 = *n - 1;
    for (j = 1; j <= i__1; ++j) {
	bounds[j] = 0.f;
/* L5: */
    }
    bounds[*n] = 1.f;
    slahqr_(&c_true, &c_true, n, &c__1, n, &workl[1], n, &ritzr[1], &ritzi[1],
	     &c__1, &c__1, &bounds[1], &c__1, ierr);
    if (*ierr != 0) {
	goto L9000;
    }

#ifndef NO_TRACE
    if (msglvl > 1) {
	svout_(n, &bounds[1], &debug_1.ndigit, "_neigh: last row of the Schur matrix for H");
    }
#endif

     /* --------------------------------------------------------- */
     /* 2. Compute the eigenvectors of the full Schur form T and  */
     /*    apply the last components of the Schur vectors to get  */
     /*    the last components of the corresponding eigenvectors. */
     /* Remember that if the i-th and (i+1)-st eigenvalues are    */
     /* complex conjugate pairs, then the real & imaginary part   */
     /* of the eigenvector components are split across adjacent   */
     /* columns of Q.                                             */
     /* --------------------------------------------------------- */

    strevc_("R", "A", select, n, &workl[1], n, vl, n, &q[q_offset], ldq, n, n,
	     &workl[*n * *n + 1], ierr);

    if (*ierr != 0) {
	goto L9000;
    }

     /* ---------------------------------------------- */
     /* Scale the returning eigenvectors so that their */
     /* euclidean norms are all one. LAPACK subroutine */
     /* strevc returns each eigenvector normalized so  */
     /* that the element of largest magnitude has      */
     /* magnitude 1; here the magnitude of a complex   */
     /* number (x,y) is taken to be |x| + |y|.         */
     /* ---------------------------------------------- */

    iconj = 0;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	if ((r__1 = ritzi[i], dabs(r__1)) <= 0.f) {

           /* -------------------- */
           /* Real eigenvalue case */
           /* -------------------- */

	    temp = snrm2_(n, &q[i * q_dim1 + 1], &c__1);
	    r__1 = 1.f / temp;
	    sscal_(n, &r__1, &q[i * q_dim1 + 1], &c__1);
	} else {

           /* ----------------------------------------- */
           /* Complex conjugate pair case. Note that    */
           /* since the real and imaginary part of      */
           /* the eigenvector are stored in consecutive */
           /* columns, we further normalize by the      */
           /* square root of two.                       */
           /* ----------------------------------------- */

	    if (iconj == 0) {
		r__1 = snrm2_(n, &q[i * q_dim1 + 1], &c__1);
		r__2 = snrm2_(n, &q[(i + 1) * q_dim1 + 1], &c__1);
		temp = slapy2_(&r__1, &r__2);
		r__1 = 1.f / temp;
		sscal_(n, &r__1, &q[i * q_dim1 + 1], &c__1);
		r__1 = 1.f / temp;
		sscal_(n, &r__1, &q[(i + 1) * q_dim1 + 1], &c__1);
		iconj = 1;
	    } else {
		iconj = 0;
	    }
	}
/* L10: */
    }

    sgemv_("T", n, n, &s_one, &q[q_offset], ldq, &bounds[1], &c__1, &s_zero, &
	    workl[1], &c__1);

#ifndef NO_TRACE
    if (msglvl > 1) {
	svout_(n, &workl[1], &debug_1.ndigit, "_neigh: Last row of the eigenvector matrix for H");
    }
#endif

     /* -------------------------- */
     /* Compute the Ritz estimates */
     /* -------------------------- */

    iconj = 0;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	if ((r__1 = ritzi[i], dabs(r__1)) <= 0.f) {

           /* -------------------- */
           /* Real eigenvalue case */
           /* -------------------- */

	    bounds[i] = *rnorm * (r__1 = workl[i], dabs(r__1));
	} else {

           /* ----------------------------------------- */
           /* Complex conjugate pair case. Note that    */
           /* since the real and imaginary part of      */
           /* the eigenvector are stored in consecutive */
           /* columns, we need to take the magnitude    */
           /* of the last components of the two vectors */
           /* ----------------------------------------- */

	    if (iconj == 0) {
		bounds[i] = *rnorm * slapy2_(&workl[i], &workl[i + 1]);
		bounds[i + 1] = bounds[i];
		iconj = 1;
	    } else {
		iconj = 0;
	    }
	}
/* L20: */
    }

#ifndef NO_TRACE
    if (msglvl > 2) {
	svout_(n, &ritzr[1], &debug_1.ndigit, "_neigh: Real part of the eigenvalues of H");
	svout_(n, &ritzi[1], &debug_1.ndigit, "_neigh: Imaginary part of the eigenvalues of H");
	svout_(n, &bounds[1], &debug_1.ndigit, "_neigh: Ritz estimates for the eigenvalues of H");
    }
#endif

#ifndef NO_TIMER
    arscnd_(&t1);
    timing_1.tneigh += t1 - t0;
#endif

L9000:
    return 0;

     /* ------------- */
     /* End of sneigh */
     /* ------------- */

} /* sneigh_ */

