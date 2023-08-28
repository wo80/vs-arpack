
#include <stdlib.h>
#include "arpack_internal.h"

void cnaupd_c(int* ido, char const* bmat, int n, char const* which, int nev, float tol, a_fcomplex* resid, int ncv, a_fcomplex* v,int ldv, int* iparam, int* ipntr, a_fcomplex* workd, a_fcomplex* workl, int lworkl, float* rwork, int* info)
{
    char w[2];

    w[0] = which[0];
    w[1] = which[1];

    cnaupd_(ido, bmat, &n, w, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, info);
}

void cneupd_c(int rvec, char const* howmny, int const* select, a_fcomplex* d, a_fcomplex* z, int ldz, a_fcomplex sigma, a_fcomplex* workev, char const* bmat, int n, char const* which, int nev, float tol, a_fcomplex* resid, int ncv, a_fcomplex* v, int ldv, int* iparam, int* ipntr, a_fcomplex* workd, a_fcomplex* workl, int lworkl, float* rwork, int* info)
{
    char w[2];

    logical rv;
    logical* slt = (logical*)malloc(ncv * sizeof(logical));

    rv = (rvec != 0) ? TRUE_ : FALSE_;

    for (int i = 0; i < ncv; i++)
    {
        slt[i] = (select[i] != 0) ? TRUE_ : FALSE_;
    }

    w[0] = which[0];
    w[1] = which[1];

    cneupd_(&rv, howmny, slt, d, z, &ldz, &sigma, workev, bmat, &n, w, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, info);
}

void dnaupd_c(int* ido, char const* bmat, int n, char const* which, int nev, double tol, double* resid, int ncv, double* v, int ldv, int* iparam, int* ipntr, double* workd, double* workl, int lworkl, int* info)
{
    char w[2];

    w[0] = which[0];
    w[1] = which[1];

    dnaupd_(ido, bmat, &n, w, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, info);
}

void dneupd_c(int rvec, char const* howmny, int const* select, double* dr, double* di, double* z, int ldz, double sigmar, double sigmai, double* workev, char const* bmat, int n, char const* which, int nev, double tol, double* resid, int ncv, double* v, int ldv, int* iparam, int* ipntr, double* workd, double* workl, int lworkl, int* info)
{
    char w[2];
    logical rv;
    logical* slt = (logical*)malloc(ncv * sizeof(logical));

    rv = (rvec != 0) ? TRUE_ : FALSE_;

    for (int i = 0; i < ncv; i++)
    {
        slt[i] = (select[i] != 0) ? TRUE_ : FALSE_;
    }

    w[0] = which[0];
    w[1] = which[1];

    dneupd_(&rv, howmny, slt, dr, di, z, &ldz, &sigmar, &sigmai, workev, bmat, &n, w, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, info);
}

void dsaupd_c(int* ido, char const* bmat, int n, char const* which, int nev, double tol, double* resid, int ncv, double* v, int ldv, int* iparam, int* ipntr, double* workd, double* workl, int lworkl, int* info)
{
    char w[2];

    w[0] = which[0];
    w[1] = which[1];

    dsaupd_(ido, bmat, &n, w, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, info);
}

void dseupd_c(int rvec, char const* howmny, int const* select, double* d, double* z, int ldz, double sigma, char const* bmat, int n, char const* which, int nev, double tol, double* resid, int ncv, double* v, int ldv, int* iparam, int* ipntr, double* workd, double* workl, int lworkl, int* info)
{
    char w[2];
    logical rv;
    logical* slt = (logical*)malloc(ncv * sizeof(logical));

    rv = (rvec != 0) ? TRUE_ : FALSE_;

    for (int i = 0; i < ncv; i++)
    {
        slt[i] = (select[i] != 0) ? TRUE_ : FALSE_;
    }

    w[0] = which[0];
    w[1] = which[1];

    dseupd_(&rv, howmny, slt, d, z, &ldz, &sigma, bmat, &n, w, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, info);
}

void snaupd_c(int* ido, char const* bmat, int n, char const* which, int nev, float tol, float* resid, int ncv, float* v, int ldv, int* iparam, int* ipntr, float* workd, float* workl, int lworkl, int* info)
{
    char w[2];

    w[0] = which[0];
    w[1] = which[1];

    snaupd_(ido, bmat, &n, w, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, info);
}

void sneupd_c(int rvec, char const* howmny, int const* select, float* dr, float* di, float* z, int ldz, float sigmar, float sigmai, float* workev, char const* bmat, int n, char const* which, int nev, float tol, float* resid, int ncv, float* v, int ldv, int* iparam, int* ipntr, float* workd, float* workl, int lworkl, int* info)
{
    char w[2];
    logical rv;
    logical* slt = (logical*)malloc(ncv * sizeof(logical));

    rv = (rvec != 0) ? TRUE_ : FALSE_;

    for (int i = 0; i < ncv; i++)
    {
        slt[i] = (select[i] != 0) ? TRUE_ : FALSE_;
    }

    w[0] = which[0];
    w[1] = which[1];

    sneupd_(&rv, howmny, slt, dr, di, z, &ldz, &sigmar, &sigmai, workev, bmat, &n, w, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, info);
}

void ssaupd_c(int* ido, char const* bmat, int n, char const* which, int nev, float tol, float* resid, int ncv, float* v, int ldv, int* iparam, int* ipntr, float* workd, float* workl, int lworkl, int* info)
{
    char w[2];

    w[0] = which[0];
    w[1] = which[1];

    ssaupd_(ido, bmat, &n, w, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, info);
}

void sseupd_c(int rvec, char const* howmny, int const* select, float* d, float* z, int ldz, float sigma, char const* bmat, int n, char const* which, int nev, float tol, float* resid, int ncv, float* v, int ldv, int* iparam, int* ipntr, float* workd, float* workl, int lworkl, int* info)
{
    char w[2];
    logical rv;
    logical* slt = (logical*)malloc(ncv * sizeof(logical));

    rv = (rvec != 0) ? TRUE_ : FALSE_;

    for (int i = 0; i < ncv; i++)
    {
        slt[i] = (select[i] != 0) ? TRUE_ : FALSE_;
    }

    w[0] = which[0];
    w[1] = which[1];

    sseupd_(&rv, howmny, slt, d, z, &ldz, &sigma, bmat, &n, w, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, info);
}


void znaupd_c(int* ido, char const* bmat, int n, char const* which, int nev, double tol, a_dcomplex* resid, int ncv, a_dcomplex* v, int ldv, int* iparam, int* ipntr, a_dcomplex* workd, a_dcomplex* workl, int lworkl, double* rwork, int* info)
{
    char w[2];

    w[0] = which[0];
    w[1] = which[1];

    znaupd_(ido, bmat, &n, w, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, info);
}

void zneupd_c(int rvec, char const* howmny, int const* select, a_dcomplex* d, a_dcomplex* z, int ldz, a_dcomplex sigma, a_dcomplex* workev, char const* bmat, int n, char const* which, int nev, double tol, a_dcomplex* resid, int ncv, a_dcomplex* v, int ldv, int* iparam, int* ipntr, a_dcomplex* workd, a_dcomplex* workl, int lworkl, double* rwork, int* info)
{
    char w[2];
    logical rv;
    logical* slt = (logical*)malloc(ncv * sizeof(logical));

    rv = (rvec != 0) ? TRUE_ : FALSE_;

    for (int i = 0; i < ncv; i++)
    {
        slt[i] = (select[i] != 0) ? TRUE_ : FALSE_;
    }

    w[0] = which[0];
    w[1] = which[1];

    zneupd_(&rv, howmny, slt, d, z, &ldz, &sigma, workev, bmat, &n, w, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, info);
}
