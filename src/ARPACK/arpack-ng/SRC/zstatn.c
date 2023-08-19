/* arpack-ng\SRC\zstatn.f -- translated by f2c (version 20100827). */

#include "arpack_internal.h"

/**
 * Initialize statistic and timing information for complex nonsymmetric Arnoldi code.
 */
int zstatn_(void)
{

    /* ------------------------------ */
    /* See stat.doc for documentation */
    /* ------------------------------ */

    timing_1.nopx = 0;
    timing_1.nbx = 0;
    timing_1.nrorth = 0;
    timing_1.nitref = 0;
    timing_1.nrstrt = 0;
    timing_1.tcaupd = 0.0f;
    timing_1.tcaup2 = 0.0f;
    timing_1.tcaitr = 0.0f;
    timing_1.tceigh = 0.0f;
    timing_1.tcgets = 0.0f;
    timing_1.tcapps = 0.0f;
    timing_1.tcconv = 0.0f;
    timing_1.titref = 0.0f;
    timing_1.tgetv0 = 0.0f;
    timing_1.trvec = 0.0f;

    /* -------------------------------------------------- */
    /* User time including reverse communication overhead */
    /* -------------------------------------------------- */

    timing_1.tmvopx = 0.0f;
    timing_1.tmvbx = 0.0f;

    return 0;
} /* zstatn_ */

