/* D:\Projekte\ARPACK\arpack-ng\SRC\dnconv.f -- translated by f2c (version 20100827). */

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
 *     xxxxxx  real
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


/* Subroutine */ int dnconv_(integer *n, doublereal *ritzr, doublereal *ritzi,
	 doublereal *bounds, doublereal *tol, integer *nconv)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    integer i__;
    static real t0, t1;
    doublereal eps23, temp;


/*     %----------------------------------------------------% */
/*     | Include files for debugging and timing information | */
/*     %----------------------------------------------------% */


/* \SCCS Information: @(#) */
/* FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 */

/*     %---------------------------------% */
/*     | See debug.doc for documentation | */
/*     %---------------------------------% */

/*     %------------------% */
/*     | Scalar Arguments | */
/*     %------------------% */

/*     %--------------------------------% */
/*     | See stat.doc for documentation | */
/*     %--------------------------------% */

/* \SCCS Information: @(#) */
/* FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 */



/*     %-----------------% */
/*     | Array Arguments | */
/*     %-----------------% */

/*     %---------------% */
/*     | Local Scalars | */
/*     %---------------% */


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

    eps23 = dlamch_("Epsilon-Machine");
    eps23 = pow_dd(&eps23, &d_23);

    *nconv = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__1 = eps23, d__2 = dlapy2_(&ritzr[i__], &ritzi[i__]);
	temp = max(d__1,d__2);
	if (bounds[i__] <= *tol * temp) {
	    ++(*nconv);
	}
/* L20: */
    }

    arscnd_(&t1);
    timing_1.tnconv += t1 - t0;

    return 0;

/*     %---------------% */
/*     | End of dnconv | */
/*     %---------------% */

} /* dnconv_ */

