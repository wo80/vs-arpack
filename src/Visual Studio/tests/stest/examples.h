#pragma once

#include "arpack.h"

int sndrv1();
int sndrv1_av_(int32_t* nx, float* v, float* w);
int sndrv1_tv_(int32_t* nx, float* x, float* y);

int sndrv2();
int sndrv2_av_(int32_t* n, float* v, float* w);

int sndrv3();
int sndrv3_av_(int32_t* n, float* v, float* w);
int sndrv3_mv_(int32_t* n, float* v, float* w);

int sndrv4();
int sndrv4_av_(int32_t* n, float* v, float* w);
int sndrv4_mv_(int32_t* n, float* v, float* w);

int sndrv5();
int sndrv5_av_(int32_t* n, float* v, float* w);
int sndrv5_mv_(int32_t* n, float* v, float* w);

int sndrv6();
int sndrv6_av_(int32_t* n, float* v, float* w);
int sndrv6_mv_(int32_t* n, float* v, float* w);

int snsimp();
int snsimp_av_(int32_t* nx, float* v, float* w);
int snsimp_tv_(int32_t* nx, float* x, float* y);

int ssdrv1();
int ssdrv1_av_(int32_t* nx, float* v, float* w);

int ssdrv2();
int ssdrv2_av_(int32_t* nx, float* v, float* w);

int ssdrv3();
int ssdrv3_av_(int32_t* nx, float* v, float* w);
int ssdrv3_mv_(int32_t* n, float* v, float* w);

int ssdrv4();
int ssdrv4_av_(int32_t* nx, float* v, float* w);
int ssdrv4_mv_(int32_t* n, float* v, float* w);

int ssdrv5();
int ssdrv5_av_(int32_t* nx, float* v, float* w);
int ssdrv5_mv_(int32_t* n, float* v, float* w);

int ssdrv6();
int ssdrv6_av_(int32_t* nx, float* v, float* w);
int ssdrv6_mv_(int32_t* n, float* v, float* w);

int sssimp();
int sssimp_av_(int32_t* nx, float* v, float* w);
int sssimp_tv_(int32_t* nx, float* x, float* y);

int ssvd();
int ssvd_av_(int32_t* m, int32_t* n, float* x, float* w);
int ssvd_atv_(int32_t* m, int32_t* n, float* w, float* y);
