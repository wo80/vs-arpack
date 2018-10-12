/* D:\Projekte\ARPACK\arpack-ng\SRC\zstatn.f -- translated by f2c (version 20100827). */

#include "arpack.h"


/* \SCCS Information: @(#) */
/* FILE: statn.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2 */

/*     %---------------------------------------------% */
/*     | Initialize statistic and timing information | */
/*     | for complex nonsymmetric Arnoldi code.      | */
/*     %---------------------------------------------% */
/* Subroutine */ int zstatn_(void)
{

/*     %--------------------------------% */
/*     | See stat.doc for documentation | */
/*     %--------------------------------% */


    timing_1.nopx = 0;
    timing_1.nbx = 0;
    timing_1.nrorth = 0;
    timing_1.nitref = 0;
    timing_1.nrstrt = 0;
    timing_1.tcaupd = 0.f;
    timing_1.tcaup2 = 0.f;
    timing_1.tcaitr = 0.f;
    timing_1.tceigh = 0.f;
    timing_1.tcgets = 0.f;
    timing_1.tcapps = 0.f;
    timing_1.tcconv = 0.f;
    timing_1.titref = 0.f;
    timing_1.tgetv0 = 0.f;
    timing_1.trvec = 0.f;
/*     %----------------------------------------------------% */
/*     | User time including reverse communication overhead | */
/*     %----------------------------------------------------% */
    timing_1.tmvopx = 0.f;
    timing_1.tmvbx = 0.f;
    return 0;

/*     %---------------% */
/*     | End of zstatn | */
/*     %---------------% */

} /* zstatn_ */

