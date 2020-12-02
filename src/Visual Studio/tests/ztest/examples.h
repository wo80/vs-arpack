#pragma once

#include "arpack.h"

int znsimp();
int znsimp_av_(integer* nx, doublecomplex* v, doublecomplex* w);
int znsimp_tv_(integer* nx, doublecomplex* x, doublecomplex* y);

int zndrv1();
int zndrv1_av_(integer* nx, doublecomplex* v, doublecomplex* w);
int zndrv1_tv_(integer* nx, doublecomplex* x, doublecomplex* y);

int zndrv2();
int zndrv2_av_(integer* n, doublecomplex* v, doublecomplex* w);

int zndrv3();
int zndrv3_av_(integer* n, doublecomplex* v, doublecomplex* w);
int zndrv3_mv_(integer* n, doublecomplex* v, doublecomplex* w);

int zndrv4();
int zndrv4_av_(integer* n, doublecomplex* v, doublecomplex* w);
int zndrv4_mv_(integer* n, doublecomplex* v, doublecomplex* w);