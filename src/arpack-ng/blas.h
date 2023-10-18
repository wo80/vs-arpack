#pragma once

#include "f2c.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /* CBLAS */

    extern void cblas_cdotc_sub(const int n, const void *x, const int incx, const void *y, const int incy, void *ret);
    extern void cblas_zdotc_sub(const int n, const void *x, const int incx, const void *y, const int incy, void *ret);

    /* BLAS (complex) */

    extern float scnrm2_(int*, a_fcomplex*, int*);
    extern void cdotc_(a_fcomplex*, int*, a_fcomplex*, int*, a_fcomplex*, int*);
    extern void caxpy_(const int*, a_fcomplex*, a_fcomplex*, int*, a_fcomplex*, int*);
    extern void ccopy_(int*, a_fcomplex*, int*, a_fcomplex*, int*);
    extern void cgemv_(char*, int*, int*, a_fcomplex*, a_fcomplex*, int*, a_fcomplex*, int*, a_fcomplex*, a_fcomplex*, int*);
    extern void cgeru_(int*, int*, a_fcomplex*, a_fcomplex*, int*, a_fcomplex*, int*, a_fcomplex*, int*);
    extern void cscal_(const int*, a_fcomplex*, a_fcomplex*, int*);
    extern void csscal_(int*, float*, a_fcomplex*, int*);
    extern void ctrmm_(char*, char*, char*, char*, int*, int*, a_fcomplex*, a_fcomplex*, int*, a_fcomplex*, int*);

    /* BLAS (double) */

    extern void daxpy_(const int*, double*, double*, int*, double*, int*);
    extern void dcopy_(int*, double*, int*, double*, int*);
    extern double ddot_(int*, double*, int*, double*, int*);
    extern void dgemv_(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
    extern void dger_(int*, int*, double*, double*, int*, double*, int*, double*, int*);
    extern double dnrm2_(int*, double*, int*);
    extern void dscal_(const int*, double*, double*, int*);
    extern void dswap_(int*, double*, int*, double*, int*);
    extern void dtrmm_(char*, char*, char*, char*, int*, int*, double*, double*, int*, double*, int*);

    /* BLAS (single) */

    extern void saxpy_(const int*, float*, float*, int*, float*, int*);
    extern void scopy_(int*, float*, int*, float*, int*);
    extern float sdot_(int*, float*, int*, float*, int*);
    extern void sgemv_(char*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);
    extern void sger_(int*, int*, float*, float*, int*, float*, int*, float*, int*);
    extern float snrm2_(int*, float*, int*);
    extern void sswap_(int*, float*, int*, float*, int*);
    extern void strmm_(char*, char*, char*, char*, int*, int*, float*, float*, int*, float*, int*);
    extern void sscal_(const int*, float*, float*, int*);

    /* BLAS (doublecomplex) */

    extern double dznrm2_(int*, a_dcomplex*, int*);
    extern void zaxpy_(const int*, a_dcomplex*, a_dcomplex*, int*, a_dcomplex*, int*);
    extern void zcopy_(int*, a_dcomplex*, int*, a_dcomplex*, int*);
    extern void zdotc_(a_dcomplex*, int*, a_dcomplex*, int*, a_dcomplex*, int*);
    extern void zgemv_(char*, int*, int*, a_dcomplex*, a_dcomplex*, int*, a_dcomplex*, int*, a_dcomplex*, a_dcomplex*, int*);
    extern void zgeru_(int*, int*, a_dcomplex*, a_dcomplex*, int*, a_dcomplex*, int*, a_dcomplex*, int*);
    extern void zscal_(const int*, a_dcomplex*, a_dcomplex*, int*);
    extern void ztrmm_(char*, char*, char*, char*, int*, int*, a_dcomplex*, a_dcomplex*, int*, a_dcomplex*, int*);

#ifdef __cplusplus
}
#endif
