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
        a_fcomplex* resid, int* ncv, a_fcomplex* v, int* ldv, int* iparam, int* ipntr,
        a_fcomplex* workd, a_fcomplex* workl, int* lworkl, float* rwork, int* info);

    int cneupd_(logical* rvec, char* howmny, logical* select, a_fcomplex* d, a_fcomplex* z, int* ldz,
        a_fcomplex* sigma, a_fcomplex* workev, char* bmat, int* n, char* which, int* nev,
        float* tol, a_fcomplex* resid, int* ncv, a_fcomplex* v, int* ldv, int* iparam,
        int* ipntr, a_fcomplex* workd, a_fcomplex* workl, int* lworkl, float* rwork,
        int* info);

    int dnaupd_(int* ido, char* bmat, int* n, char* which, int* nev, double* tol,
        double* resid, int* ncv, double* v, int* ldv, int* iparam, int* ipntr,
        double* workd, double* workl, int* lworkl, int* info);

    int dneupd_(logical* rvec, char* howmny, logical* select, double* dr, double* di, double* z,
        int* ldz, double* sigmar, double* sigmai, double* workev, char* bmat, int* n,
        char* which, int* nev, double* tol, double* resid, int* ncv, double* v,
        int* ldv, int* iparam, int* ipntr, double* workd, double* workl, int* lworkl,
        int* info);

    int dsaupd_(int* ido, char* bmat, int* n, char* which, int* nev, double* tol,
        double* resid, int* ncv, double* v, int* ldv, int* iparam, int* ipntr,
        double* workd, double* workl, int* lworkl, int* info);

    int dseupd_(logical* rvec, char* howmny, logical* select, double* d, double* z, int* ldz,
        double* sigma, char* bmat, int* n, char* which, int* nev, double* tol,
        double* resid, int* ncv, double* v, int* ldv, int* iparam, int* ipntr,
        double* workd, double* workl, int* lworkl, int* info);

    int snaupd_(int* ido, char* bmat, int* n, char* which, int* nev, float* tol,
        float* resid, int* ncv, float* v, int* ldv, int* iparam, int* ipntr,
        float* workd, float* workl, int* lworkl, int* info);

    int sneupd_(logical* rvec, char* howmny, logical* select, float* dr, float* di, float* z,
        int* ldz, float* sigmar, float* sigmai, float* workev, char* bmat, int* n,
        char* which, int* nev, float* tol, float* resid, int* ncv, float* v,
        int* ldv, int* iparam, int* ipntr, float* workd, float* workl, int* lworkl,
        int* info);

    int ssaupd_(int* ido, char* bmat, int* n, char* which, int* nev, float* tol,
        float* resid, int* ncv, float* v, int* ldv, int* iparam, int* ipntr,
        float* workd, float* workl, int* lworkl, int* info);

    int sseupd_(logical* rvec, char* howmny, logical* select, float* d, float* z, int* ldz,
        float* sigma, char* bmat, int* n, char* which, int* nev, float* tol,
        float* resid, int* ncv, float* v, int* ldv, int* iparam, int* ipntr,
        float* workd, float* workl, int* lworkl, int* info);

    int znaupd_(int* ido, char* bmat, int* n, char* which, int* nev, double* tol,
        a_dcomplex* resid, int* ncv, a_dcomplex* v, int* ldv, int* iparam, int* ipntr,
        a_dcomplex* workd, a_dcomplex* workl, int* lworkl, double* rwork, int* info);

    int zneupd_(logical* rvec, char* howmny, logical* select, a_dcomplex* d, a_dcomplex* z, int* ldz,
        a_dcomplex* sigma, a_dcomplex* workev, char* bmat, int* n, char* which, int* nev,
        double* tol, a_dcomplex* resid, int* ncv, a_dcomplex* v, int* ldv, int* iparam,
        int* ipntr, a_dcomplex* workd, a_dcomplex* workl, int* lworkl, double* rwork,
        int* info);

    /* END: public interface */

#ifdef __cplusplus
}
#endif
