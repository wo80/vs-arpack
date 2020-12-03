/* arpack-ng\SRC\dnconv.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \BeginDoc
 *
 * \Name: dnconv
 *
 * \Description:
 *  Convergence testing for the nonsymmetric Arnoldi eigenvalue routine.
 *
 * \Usage:
 *  call dnconv
 *     ( N, RITZR, RITZI, BOUNDS, TOL, NCONV )
 *
 * \Arguments
 *  N       Integer.  (INPUT)
 *          Number of Ritz values to check for convergence.
 *
 *  RITZR,  Double precision arrays of length N.  (INPUT)
 *  RITZI   Real and imaginary parts of the Ritz values to be checked
 *          for convergence.
 *  BOUNDS  Double precision array of length N.  (INPUT)
 *          Ritz estimates for the Ritz values in RITZR and RITZI.
 *
 *  TOL     Double precision scalar.  (INPUT)
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
 *     xxxxxx  float
 *
 * \Routines called:
 *     arscnd  ARPACK utility routine for timing.
 *     dlamch  LAPACK routine that determines machine constants.
 *     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
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

int dnconv_(int32_t *n, double *ritzr, double *ritzi,
	 double *bounds, double *tol, int32_t *nconv)
{
    /* System generated locals */
    int32_t i__1;
    double d__1, d__2;

    /* Builtin functions */
    double pow_dd(double *, double *);

    /* Local variables */
    int32_t i;
    static float t0, t1;
    double eps23, temp;

     /* --------------------- */
     /* Executable Statements */
     /* --------------------- */

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

    eps23 = dlamch_("Epsilon-Machine");
    eps23 = pow_dd(&eps23, &d_23);

    *nconv = 0;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
/* Computing MAX */
	d__1 = eps23, d__2 = dlapy2_(&ritzr[i], &ritzi[i]);
	temp = max(d__1,d__2);
	if (bounds[i] <= *tol * temp) {
	    ++(*nconv);
	}
/* L20: */
    }

#ifndef NO_TIMER
    arscnd_(&t1);
    timing_1.tnconv += t1 - t0;
#endif

    return 0;

     /* ------------- */
     /* End of dnconv */
     /* ------------- */

} /* dnconv_ */

