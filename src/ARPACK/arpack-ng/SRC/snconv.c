/* D:\Projekte\ARPACK\arpack-ng\SRC\snconv.f -- translated by f2c (version 20100827). */

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
 *     xxxxxx  float
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

int snconv_(int32_t *n, float *ritzr, float *ritzi, float *
	bounds, float *tol, int32_t *nconv)
{
    /* System generated locals */
    int32_t i__1;
    float r__1, r__2;
    double d__1;

    /* Builtin functions */
    double pow_dd(double *, double *);

    /* Local variables */
    int32_t i__;
    static float t0, t1;
    float eps23, temp;

/*     %-----------------% */
/*     | Array Arguments | */
/*     %-----------------% */

/*     %--------------------% */
/*     | External Functions | */
/*     %--------------------% */

/*     %-----------------------% */
/*     | Executable Statements | */
/*     %-----------------------% */

/*     %-------------------------------------------------------------% */
/*     | Convergence test: unlike in the symmetric code, I am not    | */
/*     | using things like refined error bounds and gap condition    | */
/*     | because I don't know the exact equivalent concept.          | */
/*     |                                                             | */
/*     | Instead the i-th Ritz value is considered "converged" when: | */
/*     |                                                             | */
/*     |     bounds(i) .le. ( TOL * | ritz | )                       | */
/*     |                                                             | */
/*     | for some appropriate choice of norm.                        | */
/*     %-------------------------------------------------------------% */

    /* Parameter adjustments */
    --bounds;
    --ritzi;
    --ritzr;

    /* Function Body */
    arscnd_(&t0);

/*     %---------------------------------% */
/*     | Get machine dependent constant. | */
/*     %---------------------------------% */

    eps23 = slamch_("Epsilon-Machine");
    d__1 = (double) eps23;
    eps23 = pow_dd(&d__1, &d_23);

    *nconv = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	r__1 = eps23, r__2 = slapy2_(&ritzr[i__], &ritzi[i__]);
	temp = dmax(r__1,r__2);
	if (bounds[i__] <= *tol * temp) {
	    ++(*nconv);
	}
/* L20: */
    }

    arscnd_(&t1);
    timing_1.tnconv += t1 - t0;

    return 0;

/*     %---------------% */
/*     | End of snconv | */
/*     %---------------% */

} /* snconv_ */

