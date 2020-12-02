#pragma once

#include "arpack.h"

int sndrv1();
int sndrv1_av_(integer* nx, real* v, real* w);
int sndrv1_tv_(integer* nx, real* x, real* y);

int sndrv2();
int sndrv2_av_(integer* n, real* v, real* w);

int sndrv3();
int sndrv3_av_(integer* n, real* v, real* w);
int sndrv3_mv_(integer* n, real* v, real* w);

int sndrv4();
int sndrv4_av_(integer* n, real* v, real* w);
int sndrv4_mv_(integer* n, real* v, real* w);

int sndrv5();
int sndrv5_av_(integer* n, real* v, real* w);
int sndrv5_mv_(integer* n, real* v, real* w);

int sndrv6();
int sndrv6_av_(integer* n, real* v, real* w);
int sndrv6_mv_(integer* n, real* v, real* w);

int snsimp();
int snsimp_av_(integer* nx, real* v, real* w);
int snsimp_tv_(integer* nx, real* x, real* y);

int ssdrv1();
int ssdrv1_av_(integer* nx, real* v, real* w);

int ssdrv2();
int ssdrv2_av_(integer* nx, real* v, real* w);

int ssdrv3();
int ssdrv3_av_(integer* nx, real* v, real* w);
int ssdrv3_mv_(integer* n, real* v, real* w);

int ssdrv4();
int ssdrv4_av_(integer* nx, real* v, real* w);
int ssdrv4_mv_(integer* n, real* v, real* w);

int ssdrv5();
int ssdrv5_av_(integer* nx, real* v, real* w);
int ssdrv5_mv_(integer* n, real* v, real* w);

int ssdrv6();
int ssdrv6_av_(integer* nx, real* v, real* w);
int ssdrv6_mv_(integer* n, real* v, real* w);

int sssimp();
int sssimp_av_(integer* nx, real* v, real* w);
int sssimp_tv_(integer* nx, real* x, real* y);

int ssvd();
int ssvd_av_(integer* m, integer* n, real* x, real* w);
int ssvd_atv_(integer* m, integer* n, real* w, real* y);
