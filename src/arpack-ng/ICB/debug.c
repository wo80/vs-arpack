
#include "arpack_internal.h"

void debug_c(int logfil_c, int ndigit_c, int mgetv0_c,
    int msaupd_c, int msaup2_c, int msaitr_c, int mseigt_c, int msapps_c, int msgets_c, int mseupd_c,
    int mnaupd_c, int mnaup2_c, int mnaitr_c, int mneigh_c, int mnapps_c, int mngets_c, int mneupd_c,
    int mcaupd_c, int mcaup2_c, int mcaitr_c, int mceigh_c, int mcapps_c, int mcgets_c, int mceupd_c)
{
    debug_1.logfil = logfil_c;
    debug_1.ndigit = ndigit_c;
    debug_1.mgetv0 = mgetv0_c;
    debug_1.msaupd = msaupd_c;
    debug_1.msaup2 = msaup2_c;
    debug_1.msaitr = msaitr_c;
    debug_1.mseigt = mseigt_c;
    debug_1.msapps = msapps_c;
    debug_1.msgets = msgets_c;
    debug_1.mseupd = mseupd_c;
    debug_1.mnaupd = mnaupd_c;
    debug_1.mnaup2 = mnaup2_c;
    debug_1.mnaitr = mnaitr_c;
    debug_1.mneigh = mneigh_c;
    debug_1.mnapps = mnapps_c;
    debug_1.mngets = mngets_c;
    debug_1.mneupd = mneupd_c;
    debug_1.mcaupd = mcaupd_c;
    debug_1.mcaup2 = mcaup2_c;
    debug_1.mcaitr = mcaitr_c;
    debug_1.mceigh = mceigh_c;
    debug_1.mcapps = mcapps_c;
    debug_1.mcgets = mcgets_c;
    debug_1.mceupd = mceupd_c;
}
