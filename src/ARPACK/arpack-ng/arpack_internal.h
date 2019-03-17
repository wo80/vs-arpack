#pragma once

#include "f2c.h"

#ifdef __cplusplus
extern "C"
{
#endif

	int arscnd_(real *);

	int ivout_(integer *, integer *, integer *, integer *, char *, ftnlen);

	/* complex */
	
	int cgetv0_(integer *, char *, integer *, logical *, integer *, integer *, complex *, integer *, complex *, real *, integer *, complex *, integer *);
	int cnaitr_(integer *, char *, integer *, integer *, integer *, integer *, complex *, real *, complex *, integer *, complex *, integer *, integer *, complex *, integer *);
	int cnapps_(integer *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *, complex *, complex *);
	int cnaup2_(integer *, char *, integer *, char *, integer *, integer *, real *, complex *, integer *, integer *, integer *, integer *, complex *, integer *, complex *, integer *, complex *, complex *, complex *, integer *, complex *, integer *, complex *, real *, integer *);
	int cneigh_(real *, integer *, complex *, integer *, complex *, complex *, complex *, integer *, complex *, real *, integer *);
	int cngets_(integer *, char *, integer *, integer *, complex *, complex *);
	int csortc_(char *, logical *, integer *, complex *, complex *);
	int cstatn_(void);

	int cmout_(integer *, integer *, integer *, complex *, integer *, integer *, char *, ftnlen);
	int cvout_(integer *, integer *, complex *, integer *, char *, ftnlen);

	/* double */

	int dgetv0_(integer *, char *, integer *, logical *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
	int dnaitr_(integer *, char *, integer *, integer *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *, doublereal *, integer *);
	int dnapps_(integer *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *, doublereal *);
	int dnaup2_(integer *, char *, integer *, char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *);
	int dnconv_(integer *, doublereal *, doublereal *, doublereal *, doublereal *, integer *);
	int dneigh_(doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *);
	int dngets_(integer *, char *, integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
	int dsaitr_(integer *, char *, integer *, integer *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *, doublereal *, integer *);
	int dsapps_(integer *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *);
	int dsaup2_(integer *, char *, integer *, char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *);
	int dsconv_(integer *, doublereal *, doublereal *, doublereal *, integer *);
	int dseigt_(doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *, integer *);
	int dsesrt_(char *, logical *, integer *, doublereal *, integer *, doublereal *, integer *);
	int dsgets_(integer *, char *, integer *, integer *, doublereal *, doublereal *, doublereal *);
	int dsortc_(char *, logical *, integer *, doublereal *, doublereal *, doublereal *);
	int dsortr_(char *, logical *, integer *, doublereal *, doublereal *);
	int dstatn_(void);
	int dstats_(void);
	int dstqrb_(integer *, doublereal *, doublereal *, doublereal *, doublereal *, integer *);

	int dmout_(integer *, integer *, integer *, doublereal *, integer *, integer *, char *, ftnlen);
	int dvout_(integer *, integer *, doublereal *, integer *, char *, ftnlen);

	/* single */

	int sgetv0_(integer *, char *, integer *, logical *, integer *, integer *, real *, integer *, real *, real *, integer *, real *, integer *);
	int snaitr_(integer *, char *, integer *, integer *, integer *, integer *, real *, real *, real *, integer *, real *, integer *, integer *, real *, integer *);
	int snapps_(integer *, integer *, integer *, real *, real *, real *, integer *, real *, integer *, real *, real *, integer *, real *, real *);
	int snaup2_(integer *, char *, integer *, char *, integer *, integer *, real *, real *, integer *, integer *, integer *, integer *, real *, integer *, real *, integer *, real *, real *, real *, real *, integer *, real *, integer *, real *, integer *);
	int snconv_(integer *, real *, real *, real *, real *, integer *);
	int sneigh_(real *, integer *, real *, integer *, real *, real *, real *, real *, integer *, real *, integer *);
	int sngets_(integer *, char *, integer *, integer *, real *, real *, real *, real *, real *);
	int ssaitr_(integer *, char *, integer *, integer *, integer *, integer *, real *, real *, real *, integer *, real *, integer *, integer *, real *, integer *);
	int ssapps_(integer *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *, real *);
	int ssaup2_(integer *, char *, integer *, char *, integer *, integer *, real *, real *, integer *, integer *, integer *, integer *, real *, integer *, real *, integer *, real *, real *, real *, integer *, real *, integer *, real *, integer *);
	int ssconv_(integer *, real *, real *, real *, integer *);
	int sseigt_(real *, integer *, real *, integer *, real *, real *, real *, integer *);
	int ssesrt_(char *, logical *, integer *, real *, integer *, real *, integer *);
	int ssgets_(integer *, char *, integer *, integer *, real *, real *, real *);
	int ssortc_(char *, logical *, integer *, real *, real *, real *);
	int ssortr_(char *, logical *, integer *, real *, real *);
	int sstatn_(void);
	int sstats_(void);
	int sstqrb_(integer *, real *, real *, real *, real *, integer *);

	int smout_(integer *, integer *, integer *, real *, integer *, integer *, char *, ftnlen);
	int svout_(integer *, integer *, real *, integer *, char *, ftnlen);

	/* doublecomplex */

	int zgetv0_(integer *, char *, integer *, logical *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublereal *, integer *, doublecomplex *, integer *);
	int znaitr_(integer *, char *, integer *, integer *, integer *, integer *, doublecomplex *, doublereal *, doublecomplex *, integer *, doublecomplex *, integer *, integer *, doublecomplex *, integer *);
	int znapps_(integer *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, doublecomplex *);
	int znaup2_(integer *, char *, integer *, char *, integer *, integer *, doublereal *, doublecomplex *, integer *, integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublereal *, integer *);
	int zneigh_(doublereal *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, doublereal *, integer *);
	int zngets_(integer *, char *, integer *, integer *, doublecomplex *, doublecomplex *);
	int zsortc_(char *, logical *, integer *, doublecomplex *, doublecomplex *);
	int zstatn_(void);

	int zmout_(integer *, integer *, integer *, doublecomplex *, integer *, integer *, char *, ftnlen);
	int zvout_(integer *, integer *, doublecomplex *, integer *, char *, ftnlen);

#ifdef __cplusplus
}
#endif