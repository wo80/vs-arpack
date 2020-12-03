#pragma once

#include "arpack.h"

int dndrv1();
int dndrv1_av_(int32_t* nx, double* v, double* w);
int dndrv1_tv_(int32_t* nx, double* x, double* y);

int dndrv2();
int dndrv2_av_(int32_t* n, double* v, double* w);

int dndrv3();
int dndrv3_av_(int32_t* n, double* v, double* w);
int dndrv3_mv_(int32_t* n, double* v, double* w);

int dndrv4();
int dndrv4_av_(int32_t* n, double* v, double* w);
int dndrv4_mv_(int32_t* n, double* v, double* w);

int dndrv5();
int dndrv5_av_(int32_t* n, double* v, double* w);
int dndrv5_mv_(int32_t* n, double* v, double* w);

int dndrv6();
int dndrv6_av_(int32_t* n, double* v, double* w);
int dndrv6_mv_(int32_t* n, double* v, double* w);

int dnsimp();
int dnsimp_av_(int32_t* nx, double* v, double* w);
int dnsimp_tv_(int32_t* nx, double* x, double* y);

int dsdrv1();
int dsdrv1_av_(int32_t* nx, double* v, double* w);

int dsdrv2();
int dsdrv2_av_(int32_t* nx, double* v, double* w);

int dsdrv3();
int dsdrv3_av_(int32_t* nx, double* v, double* w);
int dsdrv3_mv_(int32_t* n, double* v, double* w);

int dsdrv4();
int dsdrv4_av_(int32_t* nx, double* v, double* w);
int dsdrv4_mv_(int32_t* n, double* v, double* w);

int dsdrv5();
int dsdrv5_av_(int32_t* nx, double* v, double* w);
int dsdrv5_mv_(int32_t* n, double* v, double* w);

int dsdrv6();
int dsdrv6_av_(int32_t* nx, double* v, double* w);
int dsdrv6_mv_(int32_t* n, double* v, double* w);

int dssimp();
int dssimp_av_(int32_t* nx, double* v, double* w);
int dssimp_tv_(int32_t* nx, double* x, double* y);

int dsvd();
int dsvd_av_(int32_t* m, int32_t* n, double* x, double* w);
int dsvd_atv_(int32_t* m, int32_t* n, double* w, double* y);
