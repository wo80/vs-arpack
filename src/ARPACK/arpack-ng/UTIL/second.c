/* arpack-ng\UTIL\second.f -- translated by f2c (version 20100827). */

#include "arpack.h"

int arscnd_(float *t)
{

/*  -- LAPACK auxiliary routine (preliminary version) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     July 26, 1991 */

/*  Purpose */
/*  ======= */

/*  SECOND returns the user time for a process in arscnds. */
/*  This version gets the time from the system function ETIME. */

/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*      T1 = ETIME( TARRAY ) */
/*      T  = TARRAY( 1 ) */
    *t = 0.f;
    return 0;

/*     End of ARSCND */

} /* arscnd_ */

