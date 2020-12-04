/* arpack-ng\SRC\sstatn.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/* ------------------------------------------- */
/* Initialize statistic and timing information */
/* for nonsymmetric Arnoldi code.              */
/* ------------------------------------------- */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \SCCS Information: @(#) */
/* FILE: statn.F   SID: 2.4   DATE OF SID: 4/20/96   RELEASE: 2 */

int sstatn_(void)
{

    /* ------------------------------ */
    /* See stat.doc for documentation */
    /* ------------------------------ */

    /* --------------------- */
    /* Executable Statements */
    /* --------------------- */

    /* ------------------------------ */
    /* See stat.doc for documentation */
    /* ------------------------------ */

    timing_1.nopx = 0;
    timing_1.nbx = 0;
    timing_1.nrorth = 0;
    timing_1.nitref = 0;
    timing_1.nrstrt = 0;

    timing_1.tnaupd = 0.f;
    timing_1.tnaup2 = 0.f;
    timing_1.tnaitr = 0.f;
    timing_1.tneigh = 0.f;
    timing_1.tngets = 0.f;
    timing_1.tnapps = 0.f;
    timing_1.tnconv = 0.f;
    timing_1.titref = 0.f;
    timing_1.tgetv0 = 0.f;
    timing_1.trvec = 0.f;

    /* -------------------------------------------------- */
    /* User time including reverse communication overhead */
    /* -------------------------------------------------- */

    timing_1.tmvopx = 0.f;
    timing_1.tmvbx = 0.f;

    return 0;

    /* ------------- */
    /* End of sstatn */
    /* ------------- */

} /* sstatn_ */

