#pragma once

#include "arpack_internal.h"

int znsimp();
int znsimp_av_(const int nx, zomplex* v, zomplex* w);
int znsimp_tv_(const int nx, zomplex* x, zomplex* y);

int zndrv1();
int zndrv1_av_(const int nx, zomplex* v, zomplex* w);
int zndrv1_tv_(const int nx, zomplex* x, zomplex* y);

int zndrv2();
int zndrv2_av_(const int n, zomplex* v, zomplex* w);

int zndrv3();
int zndrv3_av_(const int n, zomplex* v, zomplex* w);
int zndrv3_mv_(const int n, zomplex* v, zomplex* w);

int zndrv4();
int zndrv4_av_(const int n, zomplex* v, zomplex* w);
int zndrv4_mv_(const int n, zomplex* v, zomplex* w);