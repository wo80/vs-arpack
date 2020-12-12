/* arpack-ng\UTIL\second.f -- translated by f2c (version 20100827). */

#ifndef NO_TIMER
#include <time.h>
#include <stdlib.h>
#endif

#include "arpack.h"

/**
 *  Returns the user time for a process in arscnds.
 */
int arscnd_(float *t)
{

    *t = 0.0f;

#ifndef NO_TIMER
    clock_t now = clock();
    *t = (float)now / CLOCKS_PER_SEC;
#endif

    return 0;
} /* arscnd_ */

