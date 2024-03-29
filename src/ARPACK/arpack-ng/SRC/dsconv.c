/* arpack-ng\SRC\dsconv.f -- translated by f2c (version 20100827). */

#include <math.h>
#include "arpack_internal.h"

/**
 * \BeginDoc
 *
 * \Name: dsconv
 *
 * \Description:
 *  Convergence testing for the symmetric Arnoldi eigenvalue routine.
 *
 * \Usage:
 *  call dsconv
 *     ( N, RITZ, BOUNDS, TOL, NCONV )
 *
 * \Arguments
 *  N       Integer.  (INPUT)
 *          Number of Ritz values to check for convergence.
 *
 *  RITZ    Double precision array of length N.  (INPUT)
 *          The Ritz values to be checked for convergence.
 *
 *  BOUNDS  Double precision array of length N.  (INPUT)
 *          Ritz estimates associated with the Ritz values in RITZ.
 *
 *  TOL     Double precision scalar.  (INPUT)
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
 *     dlamch  LAPACK routine that determines machine constants.
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
int dsconv_(int *n, double *ritz, double *bounds,
            double *tol, int *nconv)
{
    /* System generated locals */
    int i__1;
    double d__1, d__2, d__3;

    /* Builtin functions */
    /* Local variables */
    int i;
    static float t0, t1;
    double eps23, temp;

    /* Parameter adjustments */
    --bounds;
    --ritz;

    /* Function Body */
#ifndef NO_TIMER
    arscnd_(&t0);
#endif

    eps23 = dlamch_("E");
    eps23 = pow(eps23, d_23);

    *nconv = 0;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
        /* --------------------------------------------------- */
        /* The i-th Ritz value is considered "converged"       */
        /* when: bounds(i) .le. TOL*max(eps23, abs(ritz(i)))   */
        /* --------------------------------------------------- */

        /* Computing MAX */
        d__2 = eps23, d__3 = (d__1 = ritz[i], abs(d__1));
        temp = max(d__2,d__3);
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
    /* End of dsconv */
    /* ------------- */

} /* dsconv_ */

