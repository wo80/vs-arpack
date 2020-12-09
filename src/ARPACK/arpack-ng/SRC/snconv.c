/* arpack-ng\SRC\snconv.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \BeginDoc
 *
 * \Name: snconv
 *
 * \Description:
 *  Convergence testing for the nonsymmetric Arnoldi eigenvalue routine.
 *
 * \Usage:
 *  call snconv
 *     ( N, RITZR, RITZI, BOUNDS, TOL, NCONV )
 *
 * \Arguments
 *  N       Integer.  (INPUT)
 *          Number of Ritz values to check for convergence.
 *
 *  RITZR,  Real arrays of length N.  (INPUT)
 *  RITZI   Real and imaginary parts of the Ritz values to be checked
 *          for convergence.
 *  BOUNDS  Real array of length N.  (INPUT)
 *          Ritz estimates for the Ritz values in RITZR and RITZI.
 *
 *  TOL     Real scalar.  (INPUT)
 *          Desired backward error for a Ritz value to be considered
 *          "converged".
 *
 *  NCONV   Integer scalar.  (OUTPUT)
 *          Number of "converged" Ritz values.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Local variables:
 *     xxxxxx  real
 *
 * \Routines called:
 *     arscnd  ARPACK utility routine for timing.
 *     slamch  LAPACK routine that determines machine constants.
 *     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
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
 * FILE: nconv.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
 *
 * \Remarks
 *     1. xxxx
 *
 * \EndLib
 */
int snconv_(int *n, float *ritzr, float *ritzi, float *
            bounds, float *tol, int *nconv)
{
    /* System generated locals */
    int i__1;
    float r__1, r__2;
    double d__1;

    /* Builtin functions */
    double pow_dd(double *, double *);

    /* Local variables */
    int i;
    static float t0, t1;
    float eps23, temp;


    /* ----------------------------------------------------------- */
    /* Convergence test: unlike in the symmetric code, I am not    */
    /* using things like refined error bounds and gap condition    */
    /* because I don't know the exact equivalent concept.          */
    /*                                                             */
    /* Instead the i-th Ritz value is considered "converged" when: */
    /*                                                             */
    /*     bounds(i) .le. ( TOL * | ritz | )                       */
    /*                                                             */
    /* for some appropriate choice of norm.                        */
    /* ----------------------------------------------------------- */

    /* Parameter adjustments */
    --bounds;
    --ritzi;
    --ritzr;

    /* Function Body */
#ifndef NO_TIMER
    arscnd_(&t0);
#endif

    /* ------------------------------- */
    /* Get machine dependent constant. */
    /* ------------------------------- */

    eps23 = slamch_("E");
    d__1 = (double) eps23;
    eps23 = pow_dd(&d__1, &d_23);

    *nconv = 0;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
        /* Computing MAX */
        r__1 = eps23, r__2 = slapy2_(&ritzr[i], &ritzi[i]);
        temp = dmax(r__1,r__2);
        if (bounds[i] <= *tol * temp)
        {
            ++(*nconv);
        }
    }

#ifndef NO_TIMER
    arscnd_(&t1);
    timing_1.tnconv += t1 - t0;
#endif

    return 0;

    /* ------------- */
    /* End of snconv */
    /* ------------- */

} /* snconv_ */

