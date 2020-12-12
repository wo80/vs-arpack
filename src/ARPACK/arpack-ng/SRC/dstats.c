/* arpack-ng\SRC\dstats.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * Initialize statistic and timing information for symmetric Arnoldi code.
 */
int dstats_(void)
{
    /* ------------------------------ */
    /* See stat.doc for documentation */
    /* ------------------------------ */

    timing_1.nopx = 0;
    timing_1.nbx = 0;
    timing_1.nrorth = 0;
    timing_1.nitref = 0;
    timing_1.nrstrt = 0;
    timing_1.tsaupd = 0.0f;
    timing_1.tsaup2 = 0.0f;
    timing_1.tsaitr = 0.0f;
    timing_1.tseigt = 0.0f;
    timing_1.tsgets = 0.0f;
    timing_1.tsapps = 0.0f;
    timing_1.tsconv = 0.0f;
    timing_1.titref = 0.0f;
    timing_1.tgetv0 = 0.0f;
    timing_1.trvec = 0.0f;

    /* -------------------------------------------------- */
    /* User time including reverse communication overhead */
    /* -------------------------------------------------- */

    timing_1.tmvopx = 0.0f;
    timing_1.tmvbx = 0.0f;

    return 0;
} /* dstats_ */

