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

    int cnaupd_(int* ido, const char* bmat, int* n, const char* which, int* nev, float* tol,
        a_fcomplex* resid, int* ncv, a_fcomplex* v, int* ldv, int* iparam, int* ipntr,
        a_fcomplex* workd, a_fcomplex* workl, int* lworkl, float* rwork, int* info);

    int cneupd_(logical* rvec, const char* howmny, logical* select, a_fcomplex* d, a_fcomplex* z, int* ldz,
        a_fcomplex* sigma, a_fcomplex* workev, const char* bmat, int* n, const char* which, int* nev,
        float* tol, a_fcomplex* resid, int* ncv, a_fcomplex* v, int* ldv, int* iparam,
        int* ipntr, a_fcomplex* workd, a_fcomplex* workl, int* lworkl, float* rwork,
        int* info);

    int dnaupd_(int* ido, const char* bmat, int* n, const char* which, int* nev, double* tol,
        double* resid, int* ncv, double* v, int* ldv, int* iparam, int* ipntr,
        double* workd, double* workl, int* lworkl, int* info);

    int dneupd_(logical* rvec, const char* howmny, logical* select, double* dr, double* di, double* z,
        int* ldz, double* sigmar, double* sigmai, double* workev, const char* bmat, int* n,
        char* which, int* nev, double* tol, double* resid, int* ncv, double* v,
        int* ldv, int* iparam, int* ipntr, double* workd, double* workl, int* lworkl,
        int* info);

    int dsaupd_(int* ido, const char* bmat, int* n, const char* which, int* nev, double* tol,
        double* resid, int* ncv, double* v, int* ldv, int* iparam, int* ipntr,
        double* workd, double* workl, int* lworkl, int* info);

    int dseupd_(logical* rvec, const char* howmny, logical* select, double* d, double* z, int* ldz,
        double* sigma, const char* bmat, int* n, const char* which, int* nev, double* tol,
        double* resid, int* ncv, double* v, int* ldv, int* iparam, int* ipntr,
        double* workd, double* workl, int* lworkl, int* info);

    int snaupd_(int* ido, const char* bmat, int* n, const char* which, int* nev, float* tol,
        float* resid, int* ncv, float* v, int* ldv, int* iparam, int* ipntr,
        float* workd, float* workl, int* lworkl, int* info);

    int sneupd_(logical* rvec, const char* howmny, logical* select, float* dr, float* di, float* z,
        int* ldz, float* sigmar, float* sigmai, float* workev, const char* bmat, int* n,
        char* which, int* nev, float* tol, float* resid, int* ncv, float* v,
        int* ldv, int* iparam, int* ipntr, float* workd, float* workl, int* lworkl,
        int* info);

    int ssaupd_(int* ido, const char* bmat, int* n, const char* which, int* nev, float* tol,
        float* resid, int* ncv, float* v, int* ldv, int* iparam, int* ipntr,
        float* workd, float* workl, int* lworkl, int* info);

    int sseupd_(logical* rvec, const char* howmny, logical* select, float* d, float* z, int* ldz,
        float* sigma, const char* bmat, int* n, const char* which, int* nev, float* tol,
        float* resid, int* ncv, float* v, int* ldv, int* iparam, int* ipntr,
        float* workd, float* workl, int* lworkl, int* info);

    int znaupd_(int* ido, const char* bmat, int* n, const char* which, int* nev, double* tol,
        a_dcomplex* resid, int* ncv, a_dcomplex* v, int* ldv, int* iparam, int* ipntr,
        a_dcomplex* workd, a_dcomplex* workl, int* lworkl, double* rwork, int* info);

    int zneupd_(logical* rvec, const char* howmny, logical* select, a_dcomplex* d, a_dcomplex* z, int* ldz,
        a_dcomplex* sigma, a_dcomplex* workev, const char* bmat, int* n, const char* which, int* nev,
        double* tol, a_dcomplex* resid, int* ncv, a_dcomplex* v, int* ldv, int* iparam,
        int* ipntr, a_dcomplex* workd, a_dcomplex* workl, int* lworkl, double* rwork,
        int* info);

#ifdef USE_ICB

    void cnaupd_c(int* ido, char const* bmat, int n, char const* which, int nev, float tol, a_fcomplex* resid, int ncv, a_fcomplex* v, int ldv, int* iparam, int* ipntr, a_fcomplex* workd, a_fcomplex* workl, int lworkl, float* rwork, int* info);
    void cneupd_c(int rvec, char const* howmny, int const* select, a_fcomplex* d, a_fcomplex* z, int ldz, a_fcomplex sigma, a_fcomplex* workev, char const* bmat, int n, char const* which, int nev, float tol, a_fcomplex* resid, int ncv, a_fcomplex* v, int ldv, int* iparam, int* ipntr, a_fcomplex* workd, a_fcomplex* workl, int lworkl, float* rwork, int* info);
    void dnaupd_c(int* ido, char const* bmat, int n, char const* which, int nev, double tol, double* resid, int ncv, double* v, int ldv, int* iparam, int* ipntr, double* workd, double* workl, int lworkl, int* info);
    void dneupd_c(int rvec, char const* howmny, int const* select, double* dr, double* di, double* z, int ldz, double sigmar, double sigmai, double* workev, char const* bmat, int n, char const* which, int nev, double tol, double* resid, int ncv, double* v, int ldv, int* iparam, int* ipntr, double* workd, double* workl, int lworkl, int* info);
    void dsaupd_c(int* ido, char const* bmat, int n, char const* which, int nev, double tol, double* resid, int ncv, double* v, int ldv, int* iparam, int* ipntr, double* workd, double* workl, int lworkl, int* info);
    void dseupd_c(int rvec, char const* howmny, int const* select, double* d, double* z, int ldz, double sigma, char const* bmat, int n, char const* which, int nev, double tol, double* resid, int ncv, double* v, int ldv, int* iparam, int* ipntr, double* workd, double* workl, int lworkl, int* info);
    void snaupd_c(int* ido, char const* bmat, int n, char const* which, int nev, float tol, float* resid, int ncv, float* v, int ldv, int* iparam, int* ipntr, float* workd, float* workl, int lworkl, int* info);
    void sneupd_c(int rvec, char const* howmny, int const* select, float* dr, float* di, float* z, int ldz, float sigmar, float sigmai, float* workev, char const* bmat, int n, char const* which, int nev, float tol, float* resid, int ncv, float* v, int ldv, int* iparam, int* ipntr, float* workd, float* workl, int lworkl, int* info);
    void ssaupd_c(int* ido, char const* bmat, int n, char const* which, int nev, float tol, float* resid, int ncv, float* v, int ldv, int* iparam, int* ipntr, float* workd, float* workl, int lworkl, int* info);
    void sseupd_c(int rvec, char const* howmny, int const* select, float* d, float* z, int ldz, float sigma, char const* bmat, int n, char const* which, int nev, float tol, float* resid, int ncv, float* v, int ldv, int* iparam, int* ipntr, float* workd, float* workl, int lworkl, int* info);
    void znaupd_c(int* ido, char const* bmat, int n, char const* which, int nev, double tol, a_dcomplex* resid, int ncv, a_dcomplex* v, int ldv, int* iparam, int* ipntr, a_dcomplex* workd, a_dcomplex* workl, int lworkl, double* rwork, int* info);
    void zneupd_c(int rvec, char const* howmny, int const* select, a_dcomplex* d, a_dcomplex* z, int ldz, a_dcomplex sigma, a_dcomplex* workev, char const* bmat, int n, char const* which, int nev, double tol, a_dcomplex* resid, int ncv, a_dcomplex* v, int ldv, int* iparam, int* ipntr, a_dcomplex* workd, a_dcomplex* workl, int lworkl, double* rwork, int* info);

    /* Set debug levels */

    void debug_c(int logfil_c, int ndigit_c, int mgetv0_c,
        int msaupd_c, int msaup2_c, int msaitr_c, int mseigt_c, int msapps_c, int msgets_c, int mseupd_c,
        int mnaupd_c, int mnaup2_c, int mnaitr_c, int mneigh_c, int mnapps_c, int mngets_c, int mneupd_c,
        int mcaupd_c, int mcaup2_c, int mcaitr_c, int mceigh_c, int mcapps_c, int mcgets_c, int mceupd_c);

    /* Reset timers */

    void sstats_c();
    void sstatn_c();
    void cstatn_c();

    /* Get timers */

    void stat_c(int* nopx_c, int* nbx_c, int* nrorth_c, int* nitref_c, int* nrstrt_c,
        float* tsaupd_c, float* tsaup2_c, float* tsaitr_c, float* tseigt_c, float* tsgets_c, float* tsapps_c, float* tsconv_c,
        float* tnaupd_c, float* tnaup2_c, float* tnaitr_c, float* tneigh_c, float* tngets_c, float* tnapps_c, float* tnconv_c,
        float* tcaupd_c, float* tcaup2_c, float* tcaitr_c, float* tceigh_c, float* tcgets_c, float* tcapps_c, float* tcconv_c,
        float* tmvopx_c, float* tmvbx_c, float* tgetv0_c, float* titref_c, float* trvec_c);

#endif

    /* END: public interface */

#ifdef __cplusplus
}
#endif
