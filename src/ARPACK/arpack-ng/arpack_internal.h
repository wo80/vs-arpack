#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

	/* BEGIN: public interface */

	int cnaupd_(int32_t* ido, char* bmat, int32_t* n, char* which, int32_t* nev, float* tol,
		complex* resid, int32_t* ncv, complex* v, int32_t* ldv, int32_t* iparam, int32_t* ipntr,
		complex* workd, complex* workl, int32_t* lworkl, float* rwork, int32_t* info);

	int cneupd_(bool* rvec, char* howmny, bool* select, complex* d, complex* z, int32_t* ldz,
		complex* sigma, complex* workev, char* bmat, int32_t* n, char* which, int32_t* nev,
		float* tol, complex* resid, int32_t* ncv, complex* v, int32_t* ldv, int32_t* iparam,
		int32_t* ipntr, complex* workd, complex* workl, int32_t* lworkl, float* rwork,
		int32_t* info);

	int dnaupd_(int32_t* ido, char* bmat, int32_t* n, char* which, int32_t* nev, double* tol,
		double* resid, int32_t* ncv, double* v, int32_t* ldv, int32_t* iparam, int32_t* ipntr,
		double* workd, double* workl, int32_t* lworkl, int32_t* info);

	int dneupd_(bool* rvec, char* howmny, bool* select, double* dr, double* di, double* z,
		int32_t* ldz, double* sigmar, double* sigmai, double* workev, char* bmat, int32_t* n,
		char* which, int32_t* nev, double* tol, double* resid, int32_t* ncv, double* v,
		int32_t* ldv, int32_t* iparam, int32_t* ipntr, double* workd, double* workl, int32_t* lworkl,
		int32_t* info);

	int dsaupd_(int32_t* ido, char* bmat, int32_t* n, char* which, int32_t* nev, double* tol,
		double* resid, int32_t* ncv, double* v, int32_t* ldv, int32_t* iparam, int32_t* ipntr,
		double* workd, double* workl, int32_t* lworkl, int32_t* info);

	int dseupd_(bool* rvec, char* howmny, bool* select, double* d, double* z, int32_t* ldz,
		double* sigma, char* bmat, int32_t* n, char* which, int32_t* nev, double* tol,
		double* resid, int32_t* ncv, double* v, int32_t* ldv, int32_t* iparam, int32_t* ipntr,
		double* workd, double* workl, int32_t* lworkl, int32_t* info);

	int snaupd_(int32_t* ido, char* bmat, int32_t* n, char* which, int32_t* nev, float* tol,
		float* resid, int32_t* ncv, float* v, int32_t* ldv, int32_t* iparam, int32_t* ipntr,
		float* workd, float* workl, int32_t* lworkl, int32_t* info);

	int sneupd_(bool* rvec, char* howmny, bool* select, float* dr, float* di, float* z,
		int32_t* ldz, float* sigmar, float* sigmai, float* workev, char* bmat, int32_t* n,
		char* which, int32_t* nev, float* tol, float* resid, int32_t* ncv, float* v,
		int32_t* ldv, int32_t* iparam, int32_t* ipntr, float* workd, float* workl, int32_t* lworkl,
		int32_t* info);

	int ssaupd_(int32_t* ido, char* bmat, int32_t* n, char* which, int32_t* nev, float* tol,
		float* resid, int32_t* ncv, float* v, int32_t* ldv, int32_t* iparam, int32_t* ipntr,
		float* workd, float* workl, int32_t* lworkl, int32_t* info);

	int sseupd_(bool* rvec, char* howmny, bool* select, float* d, float* z, int32_t* ldz,
		float* sigma, char* bmat, int32_t* n, char* which, int32_t* nev, float* tol,
		float* resid, int32_t* ncv, float* v, int32_t* ldv, int32_t* iparam, int32_t* ipntr,
		float* workd, float* workl, int32_t* lworkl, int32_t* info);

	int znaupd_(int32_t* ido, char* bmat, int32_t* n, char* which, int32_t* nev, double* tol,
		zomplex* resid, int32_t* ncv, zomplex* v, int32_t* ldv, int32_t* iparam, int32_t* ipntr,
		zomplex* workd, zomplex* workl, int32_t* lworkl, double* rwork, int32_t* info);

	int zneupd_(bool* rvec, char* howmny, bool* select, zomplex* d, zomplex* z, int32_t* ldz,
		zomplex* sigma, zomplex* workev, char* bmat, int32_t* n, char* which, int32_t* nev,
		double* tol, zomplex* resid, int32_t* ncv, zomplex* v, int32_t* ldv, int32_t* iparam,
		int32_t* ipntr, zomplex* workd, zomplex* workl, int32_t* lworkl, double* rwork,
		int32_t* info);

	/* END: public interface */

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