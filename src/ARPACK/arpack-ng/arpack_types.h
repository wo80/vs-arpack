#pragma once

#include <stdint.h>
#include <stdbool.h>
#include "f2c.h"

#define zomplex doublecomplex

/* Common blocks */

struct {
	int32_t logfil, ndigit, mgetv0, msaupd, msaup2, msaitr, mseigt, msapps,
		msgets, mseupd, mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets,
		mneupd, mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd;
} debug_;

#define debug_1 debug_

struct {
	int32_t nopx, nbx, nrorth, nitref, nrstrt;
	float tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv, tnaupd,
		tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv, tcaupd, tcaup2,
		tcaitr, tceigh, tcgets, tcapps, tcconv, tmvopx, tmvbx, tgetv0,
		titref, trvec;
} timing_;

#define timing_1 timing_
