#pragma once

#include "arpack.h"

int cnsimp();
int cnsimp_av_(int* nx, complex* v, complex* w);
int cnsimp_tv_(int* nx, complex* x, complex* y);

int cndrv1();
int cndrv1_av_(int* nx, complex* v, complex* w);
int cndrv1_tv_(int* nx, complex* x, complex* y);

int cndrv2();
int cndrv2_av_(int* n, complex* v, complex* w);

int cndrv3();
int cndrv3_av_(int* n, complex* v, complex* w);
int cndrv3_mv_(int* n, complex* v, complex* w);

int cndrv4();
int cndrv4_av_(int* n, complex* v, complex* w);
int cndrv4_mv_(int* n, complex* v, complex* w);