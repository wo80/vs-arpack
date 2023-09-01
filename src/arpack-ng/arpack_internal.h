#pragma once

#include <stdio.h>
#include <string.h>

#include "arpack.h"
#include "lapack.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /* BEGIN: private interface */

    int arscnd_(float*);

    int ivout_(const int, const int*, const int, const char*);

    /* complex */

    int cgetv0_(int*, const char*, int*, logical*, int*, int*, a_fcomplex*, int*, a_fcomplex*, float*, int*, a_fcomplex*, int*);
    int cnaitr_(int*, const char*, int*, int*, int*, int*, a_fcomplex*, float*, a_fcomplex*, int*, a_fcomplex*, int*, int*, a_fcomplex*, int*);
    int cnapps_(int*, int*, int*, a_fcomplex*, a_fcomplex*, int*, a_fcomplex*, int*, a_fcomplex*, a_fcomplex*, int*, a_fcomplex*, a_fcomplex*);
    int cnaup2_(int*, const char*, int*, const char*, int*, int*, float*, a_fcomplex*, int*, int*, int*, int*, a_fcomplex*, int*, a_fcomplex*, int*, a_fcomplex*, a_fcomplex*, a_fcomplex*, int*, a_fcomplex*, int*, a_fcomplex*, float*, int*);
    int cneigh_(float*, int*, a_fcomplex*, int*, a_fcomplex*, a_fcomplex*, a_fcomplex*, int*, a_fcomplex*, float*, int*);
    int cngets_(int*, const char*, int*, int*, a_fcomplex*, a_fcomplex*);
    int csortc_(const char*, logical*, int*, a_fcomplex*, a_fcomplex*);
    int cstatn_(void);

    int cmout_(const int, const int, const a_fcomplex*, const int, const int, const char*);
    int cvout_(const int, const a_fcomplex*, const int, const char*);

    /* double */

    int dgetv0_(int*, const char*, int*, logical*, int*, int*, double*, int*, double*, double*, int*, double*, int*);
    int dnaitr_(int*, const char*, int*, int*, int*, int*, double*, double*, double*, int*, double*, int*, int*, double*, int*);
    int dnapps_(int*, int*, int*, double*, double*, double*, int*, double*, int*, double*, double*, int*, double*, double*);
    int dnaup2_(int*, const char*, int*, const char*, int*, int*, double*, double*, int*, int*, int*, int*, double*, int*, double*, int*, double*, double*, double*, double*, int*, double*, int*, double*, int*);
    int dnconv_(int*, double*, double*, double*, double*, int*);
    int dneigh_(double*, int*, double*, int*, double*, double*, double*, double*, int*, double*, int*);
    int dngets_(int*, const char*, int*, int*, double*, double*, double*, double*, double*);
    int dsaitr_(int*, const char*, int*, int*, int*, int*, double*, double*, double*, int*, double*, int*, int*, double*, int*);
    int dsapps_(int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*, double*);
    int dsaup2_(int*, const char*, int*, const char*, int*, int*, double*, double*, int*, int*, int*, int*, double*, int*, double*, int*, double*, double*, double*, int*, double*, int*, double*, int*);
    int dsconv_(int*, double*, double*, double*, int*);
    int dseigt_(double*, int*, double*, int*, double*, double*, double*, int*);
    int dsesrt_(const char*, logical*, int*, double*, int*, double*, int*);
    int dsgets_(int*, const char*, int*, int*, double*, double*, double*);
    int dsortc_(const char*, logical*, int*, double*, double*, double*);
    int dsortr_(const char*, logical*, int*, double*, double*);
    int dstatn_(void);
    int dstats_(void);
    int dstqrb_(int*, double*, double*, double*, double*, int*);

    int dmout_(const int, const int, const double*, const int, const int, const char*);
    int dvout_(const int, const double*, const int, const char*);

    /* single */

    int sgetv0_(int*, const char*, int*, logical*, int*, int*, float*, int*, float*, float*, int*, float*, int*);
    int snaitr_(int*, const char*, int*, int*, int*, int*, float*, float*, float*, int*, float*, int*, int*, float*, int*);
    int snapps_(int*, int*, int*, float*, float*, float*, int*, float*, int*, float*, float*, int*, float*, float*);
    int snaup2_(int*, const char*, int*, const char*, int*, int*, float*, float*, int*, int*, int*, int*, float*, int*, float*, int*, float*, float*, float*, float*, int*, float*, int*, float*, int*);
    int snconv_(int*, float*, float*, float*, float*, int*);
    int sneigh_(float*, int*, float*, int*, float*, float*, float*, float*, int*, float*, int*);
    int sngets_(int*, const char*, int*, int*, float*, float*, float*, float*, float*);
    int ssaitr_(int*, const char*, int*, int*, int*, int*, float*, float*, float*, int*, float*, int*, int*, float*, int*);
    int ssapps_(int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*, float*);
    int ssaup2_(int*, const char*, int*, const char*, int*, int*, float*, float*, int*, int*, int*, int*, float*, int*, float*, int*, float*, float*, float*, int*, float*, int*, float*, int*);
    int ssconv_(int*, float*, float*, float*, int*);
    int sseigt_(float*, int*, float*, int*, float*, float*, float*, int*);
    int ssesrt_(const char*, logical*, int*, float*, int*, float*, int*);
    int ssgets_(int*, const char*, int*, int*, float*, float*, float*);
    int ssortc_(const char*, logical*, int*, float*, float*, float*);
    int ssortr_(const char*, logical*, int*, float*, float*);
    int sstatn_(void);
    int sstats_(void);
    int sstqrb_(int*, float*, float*, float*, float*, int*);

    int smout_(const int, const int, const float*, const int, const int, const char*);
    int svout_(const int, const float*, const int, const char*);

    /* a_dcomplex */

    int zgetv0_(int*, const char*, int*, logical*, int*, int*, a_dcomplex*, int*, a_dcomplex*, double*, int*, a_dcomplex*, int*);
    int znaitr_(int*, const char*, int*, int*, int*, int*, a_dcomplex*, double*, a_dcomplex*, int*, a_dcomplex*, int*, int*, a_dcomplex*, int*);
    int znapps_(int*, int*, int*, a_dcomplex*, a_dcomplex*, int*, a_dcomplex*, int*, a_dcomplex*, a_dcomplex*, int*, a_dcomplex*, a_dcomplex*);
    int znaup2_(int*, const char*, int*, const char*, int*, int*, double*, a_dcomplex*, int*, int*, int*, int*, a_dcomplex*, int*, a_dcomplex*, int*, a_dcomplex*, a_dcomplex*, a_dcomplex*, int*, a_dcomplex*, int*, a_dcomplex*, double*, int*);
    int zneigh_(double*, int*, a_dcomplex*, int*, a_dcomplex*, a_dcomplex*, a_dcomplex*, int*, a_dcomplex*, double*, int*);
    int zngets_(int*, const char*, int*, int*, a_dcomplex*, a_dcomplex*);
    int zsortc_(const char*, logical*, int*, a_dcomplex*, a_dcomplex*);
    int zstatn_(void);

    int zmout_(const int, const int, const a_dcomplex*, const int, const int, const char*);
    int zvout_(const int, const a_dcomplex*, const int, const char*);

    /* END: private interface */

#ifdef __cplusplus
}
#endif
