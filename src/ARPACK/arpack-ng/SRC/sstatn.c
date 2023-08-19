/* arpack-ng\SRC\sstatn.f -- translated by f2c (version 20100827). */

#include "arpack_internal.h"

/**
 * Initialize statistic and timing information for nonsymmetric Arnoldi code.
 */
int sstatn_(void)
{
    /* ------------------------------ */
    /* See stat.doc for documentation */
    /* ------------------------------ */

    timing_1.nopx = 0;
    timing_1.nbx = 0;
    timing_1.nrorth = 0;
    timing_1.nitref = 0;
    timing_1.nrstrt = 0;

    timing_1.tnaupd = 0.0f;
    timing_1.tnaup2 = 0.0f;
    timing_1.tnaitr = 0.0f;
    timing_1.tneigh = 0.0f;
    timing_1.tngets = 0.0f;
    timing_1.tnapps = 0.0f;
    timing_1.tnconv = 0.0f;
    timing_1.titref = 0.0f;
    timing_1.tgetv0 = 0.0f;
    timing_1.trvec = 0.0f;

    /* -------------------------------------------------- */
    /* User time including reverse communication overhead */
    /* -------------------------------------------------- */

    timing_1.tmvopx = 0.0f;
    timing_1.tmvbx = 0.0f;

    return 0;
} /* sstatn_ */

