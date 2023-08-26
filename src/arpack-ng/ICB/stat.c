
#include "arpack_internal.h"

void sstats_c()
{
    sstats();
}

void sstatn_c()
{
    sstatn();
}

void cstatn_c()
{
    cstatn();
}

void stat_c(int* nopx_c, int* nbx_c, int* nrorth_c, int* nitref_c, int* nrstrt_c,
    float* tsaupd_c, float* tsaup2_c, float* tsaitr_c, float* tseigt_c, float* tsgets_c, float* tsapps_c, float* tsconv_c,
    float* tnaupd_c, float* tnaup2_c, float* tnaitr_c, float* tneigh_c, float* tngets_c, float* tnapps_c, float* tnconv_c,
    float* tcaupd_c, float* tcaup2_c, float* tcaitr_c, float* tceigh_c, float* tcgets_c, float* tcapps_c, float* tcconv_c,
    float* tmvopx_c, float* tmvbx_c, float* tgetv0_c, float* titref_c, float* trvec_c)
{
    *nopx_c = timing_1.nopx;
    *nbx_c = timing_1.nbx;
    *nrorth_c = timing_1.nrorth;
    *nitref_c = timing_1.nitref;
    *nrstrt_c = timing_1.nrstrt;

    *tsaupd_c = timing_1.tsaupd;
    *tsaup2_c = timing_1.tsaup2;
    *tsaitr_c = timing_1.tsaitr;
    *tseigt_c = timing_1.tseigt;
    *tsgets_c = timing_1.tsgets;
    *tsapps_c = timing_1.tsapps;
    *tsconv_c = timing_1.tsconv;
    *tnaupd_c = timing_1.tnaupd;
    *tnaup2_c = timing_1.tnaup2;
    *tnaitr_c = timing_1.tnaitr;
    *tneigh_c = timing_1.tneigh;
    *tngets_c = timing_1.tngets;
    *tnapps_c = timing_1.tnapps;
    *tnconv_c = timing_1.tnconv;
    *tcaupd_c = timing_1.tcaupd;
    *tcaup2_c = timing_1.tcaup2;
    *tcaitr_c = timing_1.tcaitr;
    *tceigh_c = timing_1.tceigh;
    *tcgets_c = timing_1.tcgets;
    *tcapps_c = timing_1.tcapps;
    *tcconv_c = timing_1.tcconv;
    *tmvopx_c = timing_1.tmvopx;
    *tmvbx_c = timing_1.tmvbx;
    *tgetv0_c = timing_1.tgetv0;
    *titref_c = timing_1.titref;
    *trvec_c = timing_1.trvec;
}
