#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

	int arscnd_(float *);

	int ivout_(int32_t *, int32_t *, int32_t *, char *);

	/* complex */
	
	int cgetv0_(int32_t *, char *, int32_t *, bool *, int32_t *, int32_t *, complex *, int32_t *, complex *, float *, int32_t *, complex *, int32_t *);
	int cnaitr_(int32_t *, char *, int32_t *, int32_t *, int32_t *, int32_t *, complex *, float *, complex *, int32_t *, complex *, int32_t *, int32_t *, complex *, int32_t *);
	int cnapps_(int32_t *, int32_t *, int32_t *, complex *, complex *, int32_t *, complex *, int32_t *, complex *, complex *, int32_t *, complex *, complex *);
	int cnaup2_(int32_t *, char *, int32_t *, char *, int32_t *, int32_t *, float *, complex *, int32_t *, int32_t *, int32_t *, int32_t *, complex *, int32_t *, complex *, int32_t *, complex *, complex *, complex *, int32_t *, complex *, int32_t *, complex *, float *, int32_t *);
	int cneigh_(float *, int32_t *, complex *, int32_t *, complex *, complex *, complex *, int32_t *, complex *, float *, int32_t *);
	int cngets_(int32_t *, char *, int32_t *, int32_t *, complex *, complex *);
	int csortc_(char *, bool *, int32_t *, complex *, complex *);
	int cstatn_(void);

	int cmout_(int32_t *, int32_t *, complex *, int32_t *, int32_t *, char *);
	int cvout_(int32_t *, complex *, int32_t *, char *);

	/* double */

	int dgetv0_(int32_t *, char *, int32_t *, bool *, int32_t *, int32_t *, double *, int32_t *, double *, double *, int32_t *, double *, int32_t *);
	int dnaitr_(int32_t *, char *, int32_t *, int32_t *, int32_t *, int32_t *, double *, double *, double *, int32_t *, double *, int32_t *, int32_t *, double *, int32_t *);
	int dnapps_(int32_t *, int32_t *, int32_t *, double *, double *, double *, int32_t *, double *, int32_t *, double *, double *, int32_t *, double *, double *);
	int dnaup2_(int32_t *, char *, int32_t *, char *, int32_t *, int32_t *, double *, double *, int32_t *, int32_t *, int32_t *, int32_t *, double *, int32_t *, double *, int32_t *, double *, double *, double *, double *, int32_t *, double *, int32_t *, double *, int32_t *);
	int dnconv_(int32_t *, double *, double *, double *, double *, int32_t *);
	int dneigh_(double *, int32_t *, double *, int32_t *, double *, double *, double *, double *, int32_t *, double *, int32_t *);
	int dngets_(int32_t *, char *, int32_t *, int32_t *, double *, double *, double *, double *, double *);
	int dsaitr_(int32_t *, char *, int32_t *, int32_t *, int32_t *, int32_t *, double *, double *, double *, int32_t *, double *, int32_t *, int32_t *, double *, int32_t *);
	int dsapps_(int32_t *, int32_t *, int32_t *, double *, double *, int32_t *, double *, int32_t *, double *, double *, int32_t *, double *);
	int dsaup2_(int32_t *, char *, int32_t *, char *, int32_t *, int32_t *, double *, double *, int32_t *, int32_t *, int32_t *, int32_t *, double *, int32_t *, double *, int32_t *, double *, double *, double *, int32_t *, double *, int32_t *, double *, int32_t *);
	int dsconv_(int32_t *, double *, double *, double *, int32_t *);
	int dseigt_(double *, int32_t *, double *, int32_t *, double *, double *, double *, int32_t *);
	int dsesrt_(char *, bool *, int32_t *, double *, int32_t *, double *, int32_t *);
	int dsgets_(int32_t *, char *, int32_t *, int32_t *, double *, double *, double *);
	int dsortc_(char *, bool *, int32_t *, double *, double *, double *);
	int dsortr_(char *, bool *, int32_t *, double *, double *);
	int dstatn_(void);
	int dstats_(void);
	int dstqrb_(int32_t *, double *, double *, double *, double *, int32_t *);

	int dmout_(int32_t *, int32_t *, double *, int32_t *, int32_t *, char *);
	int dvout_(int32_t *, double *, int32_t *, char *);

	/* single */

	int sgetv0_(int32_t *, char *, int32_t *, bool *, int32_t *, int32_t *, float *, int32_t *, float *, float *, int32_t *, float *, int32_t *);
	int snaitr_(int32_t *, char *, int32_t *, int32_t *, int32_t *, int32_t *, float *, float *, float *, int32_t *, float *, int32_t *, int32_t *, float *, int32_t *);
	int snapps_(int32_t *, int32_t *, int32_t *, float *, float *, float *, int32_t *, float *, int32_t *, float *, float *, int32_t *, float *, float *);
	int snaup2_(int32_t *, char *, int32_t *, char *, int32_t *, int32_t *, float *, float *, int32_t *, int32_t *, int32_t *, int32_t *, float *, int32_t *, float *, int32_t *, float *, float *, float *, float *, int32_t *, float *, int32_t *, float *, int32_t *);
	int snconv_(int32_t *, float *, float *, float *, float *, int32_t *);
	int sneigh_(float *, int32_t *, float *, int32_t *, float *, float *, float *, float *, int32_t *, float *, int32_t *);
	int sngets_(int32_t *, char *, int32_t *, int32_t *, float *, float *, float *, float *, float *);
	int ssaitr_(int32_t *, char *, int32_t *, int32_t *, int32_t *, int32_t *, float *, float *, float *, int32_t *, float *, int32_t *, int32_t *, float *, int32_t *);
	int ssapps_(int32_t *, int32_t *, int32_t *, float *, float *, int32_t *, float *, int32_t *, float *, float *, int32_t *, float *);
	int ssaup2_(int32_t *, char *, int32_t *, char *, int32_t *, int32_t *, float *, float *, int32_t *, int32_t *, int32_t *, int32_t *, float *, int32_t *, float *, int32_t *, float *, float *, float *, int32_t *, float *, int32_t *, float *, int32_t *);
	int ssconv_(int32_t *, float *, float *, float *, int32_t *);
	int sseigt_(float *, int32_t *, float *, int32_t *, float *, float *, float *, int32_t *);
	int ssesrt_(char *, bool *, int32_t *, float *, int32_t *, float *, int32_t *);
	int ssgets_(int32_t *, char *, int32_t *, int32_t *, float *, float *, float *);
	int ssortc_(char *, bool *, int32_t *, float *, float *, float *);
	int ssortr_(char *, bool *, int32_t *, float *, float *);
	int sstatn_(void);
	int sstats_(void);
	int sstqrb_(int32_t *, float *, float *, float *, float *, int32_t *);

	int smout_(int32_t *, int32_t *, float *, int32_t *, int32_t *, char *);
	int svout_(int32_t *, float *, int32_t *, char *);

	/* zomplex */

	int zgetv0_(int32_t *, char *, int32_t *, bool *, int32_t *, int32_t *, zomplex *, int32_t *, zomplex *, double *, int32_t *, zomplex *, int32_t *);
	int znaitr_(int32_t *, char *, int32_t *, int32_t *, int32_t *, int32_t *, zomplex *, double *, zomplex *, int32_t *, zomplex *, int32_t *, int32_t *, zomplex *, int32_t *);
	int znapps_(int32_t *, int32_t *, int32_t *, zomplex *, zomplex *, int32_t *, zomplex *, int32_t *, zomplex *, zomplex *, int32_t *, zomplex *, zomplex *);
	int znaup2_(int32_t *, char *, int32_t *, char *, int32_t *, int32_t *, double *, zomplex *, int32_t *, int32_t *, int32_t *, int32_t *, zomplex *, int32_t *, zomplex *, int32_t *, zomplex *, zomplex *, zomplex *, int32_t *, zomplex *, int32_t *, zomplex *, double *, int32_t *);
	int zneigh_(double *, int32_t *, zomplex *, int32_t *, zomplex *, zomplex *, zomplex *, int32_t *, zomplex *, double *, int32_t *);
	int zngets_(int32_t *, char *, int32_t *, int32_t *, zomplex *, zomplex *);
	int zsortc_(char *, bool *, int32_t *, zomplex *, zomplex *);
	int zstatn_(void);

	int zmout_(int32_t *, int32_t *, zomplex *, int32_t *, int32_t *, char *);
	int zvout_(int32_t *, zomplex *, int32_t *, char *);

#ifdef __cplusplus
}
#endif