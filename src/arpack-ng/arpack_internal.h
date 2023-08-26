#pragma once

#include <stdio.h>
#include <string.h>

#include "arpack.h"
#include "arpack_constants.h"
#include "lapack.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /* BEGIN: private interface */

    int arscnd_(float*);

    int ivout_(int*, int*, int*, char*);

    /* complex */

    int cgetv0_(int*, char*, int*, logical*, int*, int*, complex*, int*, complex*, float*, int*, complex*, int*);
    int cnaitr_(int*, char*, int*, int*, int*, int*, complex*, float*, complex*, int*, complex*, int*, int*, complex*, int*);
    int cnapps_(int*, int*, int*, complex*, complex*, int*, complex*, int*, complex*, complex*, int*, complex*, complex*);
    int cnaup2_(int*, char*, int*, char*, int*, int*, float*, complex*, int*, int*, int*, int*, complex*, int*, complex*, int*, complex*, complex*, complex*, int*, complex*, int*, complex*, float*, int*);
    int cneigh_(float*, int*, complex*, int*, complex*, complex*, complex*, int*, complex*, float*, int*);
    int cngets_(int*, char*, int*, int*, complex*, complex*);
    int csortc_(char*, logical*, int*, complex*, complex*);
    int cstatn_(void);

    int cmout_(int*, int*, complex*, int*, int*, char*);
    int cvout_(int*, complex*, int*, char*);

    /* double */

    int dgetv0_(int*, char*, int*, logical*, int*, int*, double*, int*, double*, double*, int*, double*, int*);
    int dnaitr_(int*, char*, int*, int*, int*, int*, double*, double*, double*, int*, double*, int*, int*, double*, int*);
    int dnapps_(int*, int*, int*, double*, double*, double*, int*, double*, int*, double*, double*, int*, double*, double*);
    int dnaup2_(int*, char*, int*, char*, int*, int*, double*, double*, int*, int*, int*, int*, double*, int*, double*, int*, double*, double*, double*, double*, int*, double*, int*, double*, int*);
    int dnconv_(int*, double*, double*, double*, double*, int*);
    int dneigh_(double*, int*, double*, int*, double*, double*, double*, double*, int*, double*, int*);
    int dngets_(int*, char*, int*, int*, double*, double*, double*, double*, double*);
    int dsaitr_(int*, char*, int*, int*, int*, int*, double*, double*, double*, int*, double*, int*, int*, double*, int*);
    int dsapps_(int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*, double*);
    int dsaup2_(int*, char*, int*, char*, int*, int*, double*, double*, int*, int*, int*, int*, double*, int*, double*, int*, double*, double*, double*, int*, double*, int*, double*, int*);
    int dsconv_(int*, double*, double*, double*, int*);
    int dseigt_(double*, int*, double*, int*, double*, double*, double*, int*);
    int dsesrt_(char*, logical*, int*, double*, int*, double*, int*);
    int dsgets_(int*, char*, int*, int*, double*, double*, double*);
    int dsortc_(char*, logical*, int*, double*, double*, double*);
    int dsortr_(char*, logical*, int*, double*, double*);
    int dstatn_(void);
    int dstats_(void);
    int dstqrb_(int*, double*, double*, double*, double*, int*);

    int dmout_(int*, int*, double*, int*, int*, char*);
    int dvout_(int*, double*, int*, char*);

    /* single */

    int sgetv0_(int*, char*, int*, logical*, int*, int*, float*, int*, float*, float*, int*, float*, int*);
    int snaitr_(int*, char*, int*, int*, int*, int*, float*, float*, float*, int*, float*, int*, int*, float*, int*);
    int snapps_(int*, int*, int*, float*, float*, float*, int*, float*, int*, float*, float*, int*, float*, float*);
    int snaup2_(int*, char*, int*, char*, int*, int*, float*, float*, int*, int*, int*, int*, float*, int*, float*, int*, float*, float*, float*, float*, int*, float*, int*, float*, int*);
    int snconv_(int*, float*, float*, float*, float*, int*);
    int sneigh_(float*, int*, float*, int*, float*, float*, float*, float*, int*, float*, int*);
    int sngets_(int*, char*, int*, int*, float*, float*, float*, float*, float*);
    int ssaitr_(int*, char*, int*, int*, int*, int*, float*, float*, float*, int*, float*, int*, int*, float*, int*);
    int ssapps_(int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*, float*);
    int ssaup2_(int*, char*, int*, char*, int*, int*, float*, float*, int*, int*, int*, int*, float*, int*, float*, int*, float*, float*, float*, int*, float*, int*, float*, int*);
    int ssconv_(int*, float*, float*, float*, int*);
    int sseigt_(float*, int*, float*, int*, float*, float*, float*, int*);
    int ssesrt_(char*, logical*, int*, float*, int*, float*, int*);
    int ssgets_(int*, char*, int*, int*, float*, float*, float*);
    int ssortc_(char*, logical*, int*, float*, float*, float*);
    int ssortr_(char*, logical*, int*, float*, float*);
    int sstatn_(void);
    int sstats_(void);
    int sstqrb_(int*, float*, float*, float*, float*, int*);

    int smout_(int*, int*, float*, int*, int*, char*);
    int svout_(int*, float*, int*, char*);

    /* zomplex */

    int zgetv0_(int*, char*, int*, logical*, int*, int*, zomplex*, int*, zomplex*, double*, int*, zomplex*, int*);
    int znaitr_(int*, char*, int*, int*, int*, int*, zomplex*, double*, zomplex*, int*, zomplex*, int*, int*, zomplex*, int*);
    int znapps_(int*, int*, int*, zomplex*, zomplex*, int*, zomplex*, int*, zomplex*, zomplex*, int*, zomplex*, zomplex*);
    int znaup2_(int*, char*, int*, char*, int*, int*, double*, zomplex*, int*, int*, int*, int*, zomplex*, int*, zomplex*, int*, zomplex*, zomplex*, zomplex*, int*, zomplex*, int*, zomplex*, double*, int*);
    int zneigh_(double*, int*, zomplex*, int*, zomplex*, zomplex*, zomplex*, int*, zomplex*, double*, int*);
    int zngets_(int*, char*, int*, int*, zomplex*, zomplex*);
    int zsortc_(char*, logical*, int*, zomplex*, zomplex*);
    int zstatn_(void);

    int zmout_(int*, int*, zomplex*, int*, int*, char*);
    int zvout_(int*, zomplex*, int*, char*);

    /* END: private interface */

#ifdef __cplusplus
}
#endif
