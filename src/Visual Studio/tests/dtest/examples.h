#pragma once

#include "arpack.h"

int dndrv1();
int dndrv1_av_(integer* nx, doublereal* v, doublereal* w);
int dndrv1_tv_(integer* nx, doublereal* x, doublereal* y);

int dndrv2();
int dndrv2_av_(integer* n, doublereal* v, doublereal* w);

int dndrv3();
int dndrv3_av_(integer* n, doublereal* v, doublereal* w);
int dndrv3_mv_(integer* n, doublereal* v, doublereal* w);

int dndrv4();
int dndrv4_av_(integer* n, doublereal* v, doublereal* w);
int dndrv4_mv_(integer* n, doublereal* v, doublereal* w);

int dndrv5();
int dndrv5_av_(integer* n, doublereal* v, doublereal* w);
int dndrv5_mv_(integer* n, doublereal* v, doublereal* w);

int dndrv6();
int dndrv6_av_(integer* n, doublereal* v, doublereal* w);
int dndrv6_mv_(integer* n, doublereal* v, doublereal* w);

int dnsimp();
int dnsimp_av_(integer* nx, doublereal* v, doublereal* w);
int dnsimp_tv_(integer* nx, doublereal* x, doublereal* y);

int dsdrv1();
int dsdrv1_av_(integer* nx, doublereal* v, doublereal* w);

int dsdrv2();
int dsdrv2_av_(integer* nx, doublereal* v, doublereal* w);

int dsdrv3();
int dsdrv3_av_(integer* nx, doublereal* v, doublereal* w);
int dsdrv3_mv_(integer* n, doublereal* v, doublereal* w);

int dsdrv4();
int dsdrv4_av_(integer* nx, doublereal* v, doublereal* w);
int dsdrv4_mv_(integer* n, doublereal* v, doublereal* w);

int dsdrv5();
int dsdrv5_av_(integer* nx, doublereal* v, doublereal* w);
int dsdrv5_mv_(integer* n, doublereal* v, doublereal* w);

int dsdrv6();
int dsdrv6_av_(integer* nx, doublereal* v, doublereal* w);
int dsdrv6_mv_(integer* n, doublereal* v, doublereal* w);

int dssimp();
int dssimp_av_(integer* nx, doublereal* v, doublereal* w);
int dssimp_tv_(integer* nx, doublereal* x, doublereal* y);

int dsvd();
int dsvd_av_(integer* m, integer* n, doublereal* x, doublereal* w);
int dsvd_atv_(integer* m, integer* n, doublereal* w, doublereal* y);
