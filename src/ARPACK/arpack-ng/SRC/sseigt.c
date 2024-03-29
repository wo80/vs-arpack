/* arpack-ng\SRC\sseigt.f -- translated by f2c (version 20100827). */

#include "arpack_internal.h"

/**
 * \BeginDoc
 *
 * \Name: sseigt
 *
 * \Description:
 *  Compute the eigenvalues of the current symmetric tridiagonal matrix
 *  and the corresponding error bounds given the current residual norm.
 *
 * \Usage:
 *  call sseigt
 *     ( RNORM, N, H, LDH, EIG, BOUNDS, WORKL, IERR )
 *
 * \Arguments
 *  RNORM   Real scalar.  (INPUT)
 *          RNORM contains the residual norm corresponding to the current
 *          symmetric tridiagonal matrix H.
 *
 *  N       Integer.  (INPUT)
 *          Size of the symmetric tridiagonal matrix H.
 *
 *  H       Real N by 2 array.  (INPUT)
 *          H contains the symmetric tridiagonal matrix with the
 *          subdiagonal in the first column starting at H(2,1) and the
 *          main diagonal in second column.
 *
 *  LDH     Integer.  (INPUT)
 *          Leading dimension of H exactly as declared in the calling
 *          program.
 *
 *  EIG     Real array of length N.  (OUTPUT)
 *          On output, EIG contains the N eigenvalues of H possibly
 *          unsorted.  The BOUNDS arrays are returned in the
 *          same sorted order as EIG.
 *
 *  BOUNDS  Real array of length N.  (OUTPUT)
 *          On output, BOUNDS contains the error estimates corresponding
 *          to the eigenvalues EIG.  This is equal to RNORM times the
 *          last components of the eigenvectors corresponding to the
 *          eigenvalues in EIG.
 *
 *  WORKL   Real work array of length 3*N.  (WORKSPACE)
 *          Private (replicated) array on each PE or array allocated on
 *          the front end.
 *
 *  IERR    Integer.  (OUTPUT)
 *          Error exit flag from sstqrb.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Local variables:
 *     xxxxxx  real
 *
 * \Routines called:
 *     sstqrb  ARPACK routine that computes the eigenvalues and the
 *             last components of the eigenvectors of a symmetric
 *             and tridiagonal matrix.
 *     arscnd  ARPACK utility routine for timing.
 *     svout   ARPACK utility routine that prints vectors.
 *     scopy   Level 1 BLAS that copies one vector to another.
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
 *     xx/xx/92: Version ' 2.4'
 *
 * \SCCS Information: @(#)
 * FILE: seigt.F   SID: 2.4   DATE OF SID: 8/27/96   RELEASE: 2
 *
 * \Remarks
 *     None
 *
 * \EndLib
 */
int sseigt_(float *rnorm, int *n, float *h, int *ldh,
            float *eig, float *bounds, float *workl, int *ierr)
{
    /* System generated locals */
    int h_dim, i__1;
    float r__1;

    /* Local variables */
    int k;
    static float t0, t1;
    int msglvl;

    /* ----------------------------- */
    /* Initialize timing statistics  */
    /* & message level for debugging */
    /* ----------------------------- */

    h_dim = *ldh;
    /* Function Body */
#ifndef NO_TIMER
    arscnd_(&t0);
#endif

    msglvl = debug_1.mseigt;

#ifndef NO_TRACE
    if (msglvl > 0)
    {
        svout_(n, &h[h_dim], &debug_1.ndigit, "_seigt: main diagonal of matrix H");
        if (*n > 1)
        {
            i__1 = *n - 1;
            svout_(&i__1, &h[1], &debug_1.ndigit, "_seigt: sub diagonal of matrix H");
        }
    }
#endif

    scopy_(n, &h[h_dim], &c__1, eig, &c__1);
    i__1 = *n - 1;
    scopy_(&i__1, &h[1], &c__1, workl, &c__1);
    sstqrb_(n, eig, workl, bounds, &workl[*n], ierr);
    if (*ierr != 0)
    {
        goto L9000;
    }
#ifndef NO_TRACE
    if (msglvl > 1)
    {
        svout_(n, bounds, &debug_1.ndigit, "_seigt: last row of the eigenvector matrix for H");
    }
#endif

    /* --------------------------------------------- */
    /* Finally determine the error bounds associated */
    /* with the n Ritz values of H.                  */
    /* --------------------------------------------- */

    i__1 = *n;
    for (k = 0; k < i__1; ++k)
    {
        bounds[k] = *rnorm * (r__1 = bounds[k], dabs(r__1));
    }

#ifndef NO_TIMER
    arscnd_(&t1);
    timing_1.tseigt += t1 - t0;
#endif

L9000:
    return 0;

    /* ------------- */
    /* End of sseigt */
    /* ------------- */

} /* sseigt_ */

