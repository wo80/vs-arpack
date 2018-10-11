/* D:\Projekte\ARPACK\arpack-ng\SRC\sstats.f -- translated by f2c (version 20100827). */

#include "arpack.h"


/* \SCCS Information: @(#) */
/* FILE: stats.F   SID: 2.1   DATE OF SID: 4/19/96   RELEASE: 2 */
/*     %---------------------------------------------% */
/*     | Initialize statistic and timing information | */
/*     | for symmetric Arnoldi code.                 | */
/*     %---------------------------------------------% */
/* Subroutine */ int sstats_(void)
{
/*     %--------------------------------% */
/*     | See stat.doc for documentation | */
/*     %--------------------------------% */
/*     %-----------------------% */
/*     | Executable Statements | */
/*     %-----------------------% */
/*     %--------------------------------% */
/*     | See stat.doc for documentation | */
/*     %--------------------------------% */

/* \SCCS Information: @(#) */
/* FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 */


    timing_1.nopx = 0;
    timing_1.nbx = 0;
    timing_1.nrorth = 0;
    timing_1.nitref = 0;
    timing_1.nrstrt = 0;
    timing_1.tsaupd = 0.f;
    timing_1.tsaup2 = 0.f;
    timing_1.tsaitr = 0.f;
    timing_1.tseigt = 0.f;
    timing_1.tsgets = 0.f;
    timing_1.tsapps = 0.f;
    timing_1.tsconv = 0.f;
    timing_1.titref = 0.f;
    timing_1.tgetv0 = 0.f;
    timing_1.trvec = 0.f;
/*     %----------------------------------------------------% */
/*     | User time including reverse communication overhead | */
/*     %----------------------------------------------------% */
    timing_1.tmvopx = 0.f;
    timing_1.tmvbx = 0.f;
    return 0;

/*     End of sstats */

} /* sstats_ */

