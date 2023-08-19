#pragma once

#include "arpack_types.h"

/* Legacy ARPACK version */
#define ARPACK_VERSION "2.4.0"
#define ARPACK_DATE "07/31/1996"

#define ARPACK_VERSION_MAJOR 2
#define ARPACK_VERSION_MINOR 4
#define ARPACK_VERSION_PATCH 0

/* ARPACK NG version */
#define arpack_ng_VERSION "3.9.1"

#define arpack_ng_VERSION_MAJOR 3
#define arpack_ng_VERSION_MINOR 9
#define arpack_ng_VERSION_PATCH 1

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

#ifdef __cplusplus
}
#endif
