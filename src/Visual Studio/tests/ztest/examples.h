#pragma once

#include "arpack.h"

int znsimp();
int znsimp_av_(int* nx, zomplex* v, zomplex* w);
int znsimp_tv_(int* nx, zomplex* x, zomplex* y);

int zndrv1();
int zndrv1_av_(int* nx, zomplex* v, zomplex* w);
int zndrv1_tv_(int* nx, zomplex* x, zomplex* y);

int zndrv2();
int zndrv2_av_(int* n, zomplex* v, zomplex* w);

int zndrv3();
int zndrv3_av_(int* n, zomplex* v, zomplex* w);
int zndrv3_mv_(int* n, zomplex* v, zomplex* w);

int zndrv4();
int zndrv4_av_(int* n, zomplex* v, zomplex* w);
int zndrv4_mv_(int* n, zomplex* v, zomplex* w);