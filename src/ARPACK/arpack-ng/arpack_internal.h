#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

	/* BEGIN: public interface */

	int cnaupd_(int* ido, char* bmat, int* n, char* which, int* nev, float* tol,
		complex* resid, int* ncv, complex* v, int* ldv, int* iparam, int* ipntr,
		complex* workd, complex* workl, int* lworkl, float* rwork, int* info);

	int cneupd_(bool* rvec, char* howmny, bool* select, complex* d, complex* z, int* ldz,
		complex* sigma, complex* workev, char* bmat, int* n, char* which, int* nev,
		float* tol, complex* resid, int* ncv, complex* v, int* ldv, int* iparam,
		int* ipntr, complex* workd, complex* workl, int* lworkl, float* rwork,
		int* info);

	int dnaupd_(int* ido, char* bmat, int* n, char* which, int* nev, double* tol,
		double* resid, int* ncv, double* v, int* ldv, int* iparam, int* ipntr,
		double* workd, double* workl, int* lworkl, int* info);

	int dneupd_(bool* rvec, char* howmny, bool* select, double* dr, double* di, double* z,
		int* ldz, double* sigmar, double* sigmai, double* workev, char* bmat, int* n,
		char* which, int* nev, double* tol, double* resid, int* ncv, double* v,
		int* ldv, int* iparam, int* ipntr, double* workd, double* workl, int* lworkl,
		int* info);

	int dsaupd_(int* ido, char* bmat, int* n, char* which, int* nev, double* tol,
		double* resid, int* ncv, double* v, int* ldv, int* iparam, int* ipntr,
		double* workd, double* workl, int* lworkl, int* info);

	int dseupd_(bool* rvec, char* howmny, bool* select, double* d, double* z, int* ldz,
		double* sigma, char* bmat, int* n, char* which, int* nev, double* tol,
		double* resid, int* ncv, double* v, int* ldv, int* iparam, int* ipntr,
		double* workd, double* workl, int* lworkl, int* info);

	int snaupd_(int* ido, char* bmat, int* n, char* which, int* nev, float* tol,
		float* resid, int* ncv, float* v, int* ldv, int* iparam, int* ipntr,
		float* workd, float* workl, int* lworkl, int* info);

	int sneupd_(bool* rvec, char* howmny, bool* select, float* dr, float* di, float* z,
		int* ldz, float* sigmar, float* sigmai, float* workev, char* bmat, int* n,
		char* which, int* nev, float* tol, float* resid, int* ncv, float* v,
		int* ldv, int* iparam, int* ipntr, float* workd, float* workl, int* lworkl,
		int* info);

	int ssaupd_(int* ido, char* bmat, int* n, char* which, int* nev, float* tol,
		float* resid, int* ncv, float* v, int* ldv, int* iparam, int* ipntr,
		float* workd, float* workl, int* lworkl, int* info);

	int sseupd_(bool* rvec, char* howmny, bool* select, float* d, float* z, int* ldz,
		float* sigma, char* bmat, int* n, char* which, int* nev, float* tol,
		float* resid, int* ncv, float* v, int* ldv, int* iparam, int* ipntr,
		float* workd, float* workl, int* lworkl, int* info);

	int znaupd_(int* ido, char* bmat, int* n, char* which, int* nev, double* tol,
		zomplex* resid, int* ncv, zomplex* v, int* ldv, int* iparam, int* ipntr,
		zomplex* workd, zomplex* workl, int* lworkl, double* rwork, int* info);

	int zneupd_(bool* rvec, char* howmny, bool* select, zomplex* d, zomplex* z, int* ldz,
		zomplex* sigma, zomplex* workev, char* bmat, int* n, char* which, int* nev,
		double* tol, zomplex* resid, int* ncv, zomplex* v, int* ldv, int* iparam,
		int* ipntr, zomplex* workd, zomplex* workl, int* lworkl, double* rwork,
		int* info);

	/* END: public interface */

	int arscnd_(float *);

	int ivout_(int *, int *, int *, char *);

	/* complex */
	
	int cgetv0_(int *, char *, int *, bool *, int *, int *, complex *, int *, complex *, float *, int *, complex *, int *);
	int cnaitr_(int *, char *, int *, int *, int *, int *, complex *, float *, complex *, int *, complex *, int *, int *, complex *, int *);
	int cnapps_(int *, int *, int *, complex *, complex *, int *, complex *, int *, complex *, complex *, int *, complex *, complex *);
	int cnaup2_(int *, char *, int *, char *, int *, int *, float *, complex *, int *, int *, int *, int *, complex *, int *, complex *, int *, complex *, complex *, complex *, int *, complex *, int *, complex *, float *, int *);
	int cneigh_(float *, int *, complex *, int *, complex *, complex *, complex *, int *, complex *, float *, int *);
	int cngets_(int *, char *, int *, int *, complex *, complex *);
	int csortc_(char *, bool *, int *, complex *, complex *);
	int cstatn_(void);

	int cmout_(int *, int *, complex *, int *, int *, char *);
	int cvout_(int *, complex *, int *, char *);

	/* double */

	int dgetv0_(int *, char *, int *, bool *, int *, int *, double *, int *, double *, double *, int *, double *, int *);
	int dnaitr_(int *, char *, int *, int *, int *, int *, double *, double *, double *, int *, double *, int *, int *, double *, int *);
	int dnapps_(int *, int *, int *, double *, double *, double *, int *, double *, int *, double *, double *, int *, double *, double *);
	int dnaup2_(int *, char *, int *, char *, int *, int *, double *, double *, int *, int *, int *, int *, double *, int *, double *, int *, double *, double *, double *, double *, int *, double *, int *, double *, int *);
	int dnconv_(int *, double *, double *, double *, double *, int *);
	int dneigh_(double *, int *, double *, int *, double *, double *, double *, double *, int *, double *, int *);
	int dngets_(int *, char *, int *, int *, double *, double *, double *, double *, double *);
	int dsaitr_(int *, char *, int *, int *, int *, int *, double *, double *, double *, int *, double *, int *, int *, double *, int *);
	int dsapps_(int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *, double *);
	int dsaup2_(int *, char *, int *, char *, int *, int *, double *, double *, int *, int *, int *, int *, double *, int *, double *, int *, double *, double *, double *, int *, double *, int *, double *, int *);
	int dsconv_(int *, double *, double *, double *, int *);
	int dseigt_(double *, int *, double *, int *, double *, double *, double *, int *);
	int dsesrt_(char *, bool *, int *, double *, int *, double *, int *);
	int dsgets_(int *, char *, int *, int *, double *, double *, double *);
	int dsortc_(char *, bool *, int *, double *, double *, double *);
	int dsortr_(char *, bool *, int *, double *, double *);
	int dstatn_(void);
	int dstats_(void);
	int dstqrb_(int *, double *, double *, double *, double *, int *);

	int dmout_(int *, int *, double *, int *, int *, char *);
	int dvout_(int *, double *, int *, char *);

	/* single */

	int sgetv0_(int *, char *, int *, bool *, int *, int *, float *, int *, float *, float *, int *, float *, int *);
	int snaitr_(int *, char *, int *, int *, int *, int *, float *, float *, float *, int *, float *, int *, int *, float *, int *);
	int snapps_(int *, int *, int *, float *, float *, float *, int *, float *, int *, float *, float *, int *, float *, float *);
	int snaup2_(int *, char *, int *, char *, int *, int *, float *, float *, int *, int *, int *, int *, float *, int *, float *, int *, float *, float *, float *, float *, int *, float *, int *, float *, int *);
	int snconv_(int *, float *, float *, float *, float *, int *);
	int sneigh_(float *, int *, float *, int *, float *, float *, float *, float *, int *, float *, int *);
	int sngets_(int *, char *, int *, int *, float *, float *, float *, float *, float *);
	int ssaitr_(int *, char *, int *, int *, int *, int *, float *, float *, float *, int *, float *, int *, int *, float *, int *);
	int ssapps_(int *, int *, int *, float *, float *, int *, float *, int *, float *, float *, int *, float *);
	int ssaup2_(int *, char *, int *, char *, int *, int *, float *, float *, int *, int *, int *, int *, float *, int *, float *, int *, float *, float *, float *, int *, float *, int *, float *, int *);
	int ssconv_(int *, float *, float *, float *, int *);
	int sseigt_(float *, int *, float *, int *, float *, float *, float *, int *);
	int ssesrt_(char *, bool *, int *, float *, int *, float *, int *);
	int ssgets_(int *, char *, int *, int *, float *, float *, float *);
	int ssortc_(char *, bool *, int *, float *, float *, float *);
	int ssortr_(char *, bool *, int *, float *, float *);
	int sstatn_(void);
	int sstats_(void);
	int sstqrb_(int *, float *, float *, float *, float *, int *);

	int smout_(int *, int *, float *, int *, int *, char *);
	int svout_(int *, float *, int *, char *);

	/* zomplex */

	int zgetv0_(int *, char *, int *, bool *, int *, int *, zomplex *, int *, zomplex *, double *, int *, zomplex *, int *);
	int znaitr_(int *, char *, int *, int *, int *, int *, zomplex *, double *, zomplex *, int *, zomplex *, int *, int *, zomplex *, int *);
	int znapps_(int *, int *, int *, zomplex *, zomplex *, int *, zomplex *, int *, zomplex *, zomplex *, int *, zomplex *, zomplex *);
	int znaup2_(int *, char *, int *, char *, int *, int *, double *, zomplex *, int *, int *, int *, int *, zomplex *, int *, zomplex *, int *, zomplex *, zomplex *, zomplex *, int *, zomplex *, int *, zomplex *, double *, int *);
	int zneigh_(double *, int *, zomplex *, int *, zomplex *, zomplex *, zomplex *, int *, zomplex *, double *, int *);
	int zngets_(int *, char *, int *, int *, zomplex *, zomplex *);
	int zsortc_(char *, bool *, int *, zomplex *, zomplex *);
	int zstatn_(void);

	int zmout_(int *, int *, zomplex *, int *, int *, char *);
	int zvout_(int *, zomplex *, int *, char *);

#ifdef __cplusplus
}
#endif