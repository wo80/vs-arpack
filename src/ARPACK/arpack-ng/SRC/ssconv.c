/* arpack-ng\SRC\ssconv.f -- translated by f2c (version 20100827). */

#include <math.h>
#include "arpack_internal.h"

/**
 * \BeginDoc
 *
 * \Name: ssconv
 *
 * \Description:
 *  Convergence testing for the symmetric Arnoldi eigenvalue routine.
 *
 * \Usage:
 *  call ssconv
 *     ( N, RITZ, BOUNDS, TOL, NCONV )
 *
 * \Arguments
 *  N       Integer.  (INPUT)
 *          Number of Ritz values to check for convergence.
 *
 *  RITZ    Real array of length N.  (INPUT)
 *          The Ritz values to be checked for convergence.
 *
 *  BOUNDS  Real array of length N.  (INPUT)
 *          Ritz estimates associated with the Ritz values in RITZ.
 *
 *  TOL     Real scalar.  (INPUT)
 *          Desired relative accuracy for a Ritz value to be considered
 *          "converged".
 *
 *  NCONV   Integer scalar.  (OUTPUT)
 *          Number of "converged" Ritz values.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Routines called:
 *     arscnd  ARPACK utility routine for timing.
 *     slamch  LAPACK routine that determines machine constants.
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
 * FILE: sconv.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2
 *
 * \Remarks
 *     1. Starting with version 2.4, this routine no longer uses the
 *        Parlett strategy using the gap conditions.
 *
 * \EndLib
 */
int ssconv_(int *n, float *ritz, float *bounds, float *tol,
            int *nconv)
{
    /* System generated locals */
    int i__1;
    float r__1, r__2, r__3;

    /* Builtin functions */
    /* Local variables */
    int i;
    static float t0, t1;
    float eps23, temp;

    /* Parameter adjustments */
    --bounds;
    --ritz;

    /* Function Body */
#ifndef NO_TIMER
    arscnd_(&t0);
#endif

    eps23 = slamch_("E");
    eps23 = pow((double)eps23, d_23);

    *nconv = 0;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
        /* --------------------------------------------------- */
        /* The i-th Ritz value is considered "converged"       */
        /* when: bounds(i) .le. TOL*max(eps23, abs(ritz(i)))   */
        /* --------------------------------------------------- */

        /* Computing MAX */
        r__2 = eps23, r__3 = (r__1 = ritz[i], dabs(r__1));
        temp = dmax(r__2,r__3);
        if (bounds[i] <= *tol * temp)
        {
            ++(*nconv);
        }
    }

#ifndef NO_TIMER
    arscnd_(&t1);
    timing_1.tsconv += t1 - t0;
#endif

    return 0;

    /* ------------- */
    /* End of ssconv */
    /* ------------- */

} /* ssconv_ */

