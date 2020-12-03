/* D:\Projekte\ARPACK\arpack-ng\SRC\zneigh.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \BeginDoc
 *
 * \Name: zneigh
 *
 * \Description:
 *  Compute the eigenvalues of the current upper Hessenberg matrix
 *  and the corresponding Ritz estimates given the current residual norm.
 *
 * \Usage:
 *  call zneigh
 *     ( RNORM, N, H, LDH, RITZ, BOUNDS, Q, LDQ, WORKL, RWORK, IERR )
 *
 * \Arguments
 *  RNORM   Double precision scalar.  (INPUT)
 *          Residual norm corresponding to the current upper Hessenberg
 *          matrix H.
 *
 *  N       Integer.  (INPUT)
 *          Size of the matrix H.
 *
 *  H       Complex*16 N by N array.  (INPUT)
 *          H contains the current upper Hessenberg matrix.
 *
 *  LDH     Integer.  (INPUT)
 *          Leading dimension of H exactly as declared in the calling
 *          program.
 *
 *  RITZ    Complex*16 array of length N.  (OUTPUT)
 *          On output, RITZ(1:N) contains the eigenvalues of H.
 *
 *  BOUNDS  Complex*16 array of length N.  (OUTPUT)
 *          On output, BOUNDS contains the Ritz estimates associated with
 *          the eigenvalues held in RITZ.  This is equal to RNORM
 *          times the last components of the eigenvectors corresponding
 *          to the eigenvalues in RITZ.
 *
 *  Q       Complex*16 N by N array.  (WORKSPACE)
 *          Workspace needed to store the eigenvectors of H.
 *
 *  LDQ     Integer.  (INPUT)
 *          Leading dimension of Q exactly as declared in the calling
 *          program.
 *
 *  WORKL   Complex*16 work array of length N**2 + 3*N.  (WORKSPACE)
 *          Private (replicated) array on each PE or array allocated on
 *          the front end.  This is needed to keep the full Schur form
 *          of H and also in the calculation of the eigenvectors of H.
 *
 *  RWORK   Double precision  work array of length N (WORKSPACE)
 *          Private (replicated) array on each PE or array allocated on
 *          the front end.
 *
 *  IERR    Integer.  (OUTPUT)
 *          Error exit flag from zlahqr or ztrevc.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Local variables:
 *     xxxxxx  Complex*16
 *
 * \Routines called:
 *     ivout   ARPACK utility routine that prints int32_ts.
 *     arscnd  ARPACK utility routine for timing.
 *     zmout   ARPACK utility routine that prints matrices
 *     zvout   ARPACK utility routine that prints vectors.
 *     dvout   ARPACK utility routine that prints vectors.
 *     zlacpy  LAPACK matrix copy routine.
 *     zlahqr  LAPACK routine to compute the Schur form of an
 *             upper Hessenberg matrix.
 *     zlaset  LAPACK matrix initialization routine.
 *     ztrevc  LAPACK routine to compute the eigenvectors of a matrix
 *             in upper triangular form
 *     zcopy   Level 1 BLAS that copies one vector to another.
 *     zdscal  Level 1 BLAS that scales a complex vector by a float number.
 *     dznrm2  Level 1 BLAS that computes the norm of a vector.
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
 * \SCCS Information: @(#)
 * FILE: neigh.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2
 *
 * \Remarks
 *     None
 *
 * \EndLib
 */

int zneigh_(double *rnorm, int32_t *n, zomplex *
	h__, int32_t *ldh, zomplex *ritz, zomplex *bounds, 
	zomplex *q, int32_t *ldq, zomplex *workl, double *
	rwork, int32_t *ierr)
{
    /* System generated locals */
    int32_t h_dim1, h_offset, q_dim1, q_offset, i__1;
    double d__1;

    /* Local variables */
    int32_t j;
    static float t0, t1;
    zomplex vl[1];
    double temp;
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
    --rwork;
    --workl;
    --bounds;
    --ritz;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;

    /* Function Body */
    arscnd_(&t0);
    msglvl = debug_1.mceigh;

    if (msglvl > 2) {
	zmout_(&debug_1.logfil, n, n, &h__[h_offset], ldh, &debug_1.ndigit, "_neigh: Entering upper Hessenberg matrix H ");
    }

     /* -------------------------------------------------------- */
     /* 1. Compute the eigenvalues, the last components of the   */
     /*    corresponding Schur vectors and the full Schur form T */
     /*    of the current upper Hessenberg matrix H.             */
     /*    zlahqr returns the full Schur form of H               */
     /*    in WORKL(1:N**2), and the Schur vectors in q.         */
     /* -------------------------------------------------------- */

    zlacpy_("All", n, n, &h__[h_offset], ldh, &workl[1], n);
    zlaset_("All", n, n, &z_zero, &z_one, &q[q_offset], ldq);
    zlahqr_(&c_true, &c_true, n, &c__1, n, &workl[1], ldh, &ritz[1], &c__1, n,
	     &q[q_offset], ldq, ierr);
    if (*ierr != 0) {
	goto L9000;
    }

    zcopy_(n, &q[*n - 1 + q_dim1], ldq, &bounds[1], &c__1);
    if (msglvl > 1) {
	zvout_(&debug_1.logfil, n, &bounds[1], &debug_1.ndigit, "_neigh: last row of the Schur matrix for H");
    }

     /* -------------------------------------------------------- */
     /* 2. Compute the eigenvectors of the full Schur form T and */
     /*    apply the Schur vectors to get the corresponding      */
     /*    eigenvectors.                                         */
     /* -------------------------------------------------------- */

    ztrevc_("Right", "Back", select, n, &workl[1], n, vl, n, &q[q_offset], 
	    ldq, n, n, &workl[*n * *n + 1], &rwork[1], ierr);

    if (*ierr != 0) {
	goto L9000;
    }

     /* ---------------------------------------------- */
     /* Scale the returning eigenvectors so that their */
     /* Euclidean norms are all one. LAPACK subroutine */
     /* ztrevc returns each eigenvector normalized so  */
     /* that the element of largest magnitude has      */
     /* magnitude 1; here the magnitude of a complex   */
     /* number (x,y) is taken to be |x| + |y|.         */
     /* ---------------------------------------------- */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	temp = dznrm2_(n, &q[j * q_dim1 + 1], &c__1);
	d__1 = 1. / temp;
	zdscal_(n, &d__1, &q[j * q_dim1 + 1], &c__1);
/* L10: */
    }

    if (msglvl > 1) {
	zcopy_(n, &q[*n + q_dim1], ldq, &workl[1], &c__1);
	zvout_(&debug_1.logfil, n, &workl[1], &debug_1.ndigit, "_neigh: Last row of the eigenvector matrix for H");
    }

     /* -------------------------- */
     /* Compute the Ritz estimates */
     /* -------------------------- */

    zcopy_(n, &q[*n + q_dim1], n, &bounds[1], &c__1);
    zdscal_(n, rnorm, &bounds[1], &c__1);

    if (msglvl > 2) {
	zvout_(&debug_1.logfil, n, &ritz[1], &debug_1.ndigit, "_neigh: The eigenvalues of H");
	zvout_(&debug_1.logfil, n, &bounds[1], &debug_1.ndigit, "_neigh: Ritz estimates for the eigenvalues of H");
    }

    arscnd_(&t1);
    timing_1.tceigh += t1 - t0;

L9000:
    return 0;

     /* ------------- */
     /* End of zneigh */
     /* ------------- */

} /* zneigh_ */

