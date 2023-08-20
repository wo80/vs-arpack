#pragma once

#include "f2c.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /* complex */

    extern double scnrm2_(int *, complex *, int *);
    extern int caxpy_(int *, complex *, complex *, int *, complex *, int *);
    extern int ccopy_(int *, complex *, int *, complex *, int *);
    extern int cgemv_(char *, int *, int *, complex *, complex *, int *, complex *, int *, complex *, complex *, int *);
    extern int cgeqr2_(int *, int *, complex *, int *, complex *, complex *, int *);
    extern int cgeru_(int *, int *, complex *, complex *, int *, complex *, int *, complex *, int *);
    extern int clacpy_(char *, int *, int *, complex *, int *, complex *, int *);
    extern int clahqr_(logical *, logical *, int *, int *, int *, complex *, int *, complex *, int *, int *, complex *, int *, int *);
    extern double clanhs_(char *, int *, complex *, int *, complex *);
    extern int clarnv_(int *, int *, int *, complex *);
    extern int clartg_(complex *, complex *, float *, complex *, complex *);
    extern int clascl_(char *, int *, int *, float *, float *, int *, int *, complex *, int *, int *);
    extern int claset_(char *, int *, int *, complex *, complex *, complex *, int *);
    extern /* Complex */ VOID cdotc_(complex *, int *, complex *, int *, complex *, int *);
    extern int cscal_(int *, complex *, complex *, int *);
    extern int csscal_(int *, float *, complex *, int *);
    extern int ctrevc_(char *, char *, logical *, int *, complex *, int *, complex *, int *, complex *, int *, int *, int *, complex *, float *, int *);
    extern int ctrmm_(char *, char *, char *, char *, int *, int *, complex *, complex *, int *, complex *, int *);
    extern int ctrsen_(char *, char *, logical *, int *, complex *, int *, complex *, int *, complex *, int *, float *, float *, complex *, int *, int *);
    extern int cunm2r_(char *, char *, int *, int *, int *, complex *, int *, complex *, complex *, int *, complex *, int *);

    /* double */

    extern int daxpy_(int *, double *, double *, int *, double *, int *);
    extern int dcopy_(int *, double *, int *, double *, int *);
    extern double ddot_(int *, double *, int *, double *, int *);
    extern int dgemv_(char *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
    extern int dgeqr2_(int *, int *, double *, int *, double *, double *, int *);
    extern int dger_(int *, int *, double *, double *, int *, double *, int *, double *, int *);
    extern int dlabad_(double *, double *);
    extern int dlacpy_(char *, int *, int *, double *, int *, double *, int *);
    extern int dlae2_(double *, double *, double *, double *, double *);
    extern int dlaev2_(double *, double *, double *, double *, double *, double *, double *);
    extern int dlahqr_(logical *, logical *, int *, int *, int *, double *, int *, double *, double *, int *, int *, double *, int *, int *);
    extern double dlamch_(char *);
    extern double dlanhs_(char *, int *, double *, int *, double *);
    extern double dlanst_(char *, int *, double *, double *);
    extern double dlapy2_(double *, double *);
    extern int dlarf_(char *, int *, int *, double *, int *, double *, double *, int *, double *);
    extern int dlarfg_(int *, double *, double *, int *, double *);
    extern int dlarnv_(int *, int *, int *, double *);
    extern int dlartg_(double *, double *, double *, double *, double *);
    extern int dlascl_(char *, int *, int *, double *, double *, int *, int *, double *, int *, int *);
    extern int dlaset_(char *, int *, int *, double *, double *, double *, int *);
    extern int dlasr_(char *, char *, char *, int *, int *, double *, double *, double *, int *);
    extern int dlasrt_(char *, int *, double *, int *);
    extern double dnrm2_(int *, double *, int *);
    extern int dorm2r_(char *, char *, int *, int *, int *, double *, int *, double *, double *, int *, double *, int *);
    extern int dscal_(int *, double *, double *, int *);
    extern int dsteqr_(char *, int *, double *, double *, double *, int *, double *, int *);
    extern int dswap_(int *, double *, int *, double *, int *);
    extern int dtrevc_(char *, char *, logical *, int *, double *, int *, double *, int *, double *, int *, int *, int *, double *, int *);
    extern int dtrmm_(char *, char *, char *, char *, int *, int *, double *, double *, int *, double *, int *);
    extern int dtrsen_(char *, char *, logical *, int *, double *, int *, double *, int *, double *, double *, int *, double *, double *, double *, int *, int *, int *, int *);

    /* single */

    extern int saxpy_(int *, float *, float *, int *, float *, int *);
    extern int scopy_(int *, float *, int *, float *, int *);
    extern double sdot_(int *, float *, int *, float *, int *);
    extern int sgemv_(char *, int *, int *, float *, float *, int *, float *, int *, float *, float *, int *);
    extern int sgeqr2_(int *, int *, float *, int *, float *, float *, int *);
    extern int sger_(int *, int *, float *, float *, int *, float *, int *, float *, int *);
    extern int slabad_(float *, float *);
    extern int slacpy_(char *, int *, int *, float *, int *, float *, int *);
    extern int slae2_(float *, float *, float *, float *, float *);
    extern int slaev2_(float *, float *, float *, float *, float *, float *, float *);
    extern int slahqr_(logical *, logical *, int *, int *, int *, float *, int *, float *, float *, int *, int *, float *, int *, int *);
    extern double slamch_(char *);
    extern double slanhs_(char *, int *, float *, int *, float *);
    extern double slanst_(char *, int *, float *, float *);
    extern double slapy2_(float *, float *);
    extern int slarf_(char *, int *, int *, float *, int *, float *, float *, int *, float *);
    extern int slarfg_(int *, float *, float *, int *, float *);
    extern int slarnv_(int *, int *, int *, float *);
    extern int slartg_(float *, float *, float *, float *, float *);
    extern int slascl_(char *, int *, int *, float *, float *, int *, int *, float *, int *, int *);
    extern int slaset_(char *, int *, int *, float *, float *, float *, int *);
    extern int slasr_(char *, char *, char *, int *, int *, float *, float *, float *, int *);
    extern int slasrt_(char *, int *, float *, int *);
    extern double snrm2_(int *, float *, int *);
    extern int sorm2r_(char *, char *, int *, int *, int *, float *, int *, float *, float *, int *, float *, int *);
    extern int sscal_(int *, float *, float *, int *);
    extern int ssteqr_(char *, int *, float *, float *, float *, int *, float *, int *);
    extern int sswap_(int *, float *, int *, float *, int *);
    extern int strevc_(char *, char *, logical *, int *, float *, int *, float *, int *, float *, int *, int *, int *, float *, int *);
    extern int strmm_(char *, char *, char *, char *, int *, int *, float *, float *, int *, float *, int *);
    extern int strsen_(char *, char *, logical *, int *, float *, int *, float *, int *, float *, float *, int *, float *, float *, float *, int *, int *, int *, int *);

    /* doublecomplex */

    extern double dznrm2_(int *, doublecomplex *, int *);
    extern int zaxpy_(int *, doublecomplex *, doublecomplex *, int *, doublecomplex *, int *);
    extern int zcopy_(int *, doublecomplex *, int *, doublecomplex *, int *);
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, int *, doublecomplex *, int *, doublecomplex *, int *);
    extern int zdscal_(int *, double *, doublecomplex *, int *);
    extern int zgemv_(char *, int *, int *, doublecomplex *, doublecomplex *, int *, doublecomplex *, int *, doublecomplex *, doublecomplex *, int *);
    extern int zgeqr2_(int *, int *, doublecomplex *, int *, doublecomplex *, doublecomplex *, int *);
    extern int zgeru_(int *, int *, doublecomplex *, doublecomplex *, int *, doublecomplex *, int *, doublecomplex *, int *);
    extern int zlacpy_(char *, int *, int *, doublecomplex *, int *, doublecomplex *, int *);
    extern int zlahqr_(logical *, logical *, int *, int *, int *, doublecomplex *, int *, doublecomplex *, int *, int *, doublecomplex *, int *, int *);
    extern double zlanhs_(char *, int *, doublecomplex *, int *, doublecomplex *);
    extern int zlarnv_(int *, int *, int *, doublecomplex *);
    extern int zlartg_(doublecomplex *, doublecomplex *, double *, doublecomplex *, doublecomplex *);
    extern int zlascl_(char *, int *, int *, double *, double *, int *, int *, doublecomplex *, int *, int *);
    extern int zlaset_(char *, int *, int *, doublecomplex *, doublecomplex *, doublecomplex *, int *);
    extern int zscal_(int *, doublecomplex *, doublecomplex *, int *);
    extern int ztrevc_(char *, char *, logical *, int *, doublecomplex *, int *, doublecomplex *, int *, doublecomplex *, int *, int *, int *, doublecomplex *, double *, int *);
    extern int ztrmm_(char *, char *, char *, char *, int *, int *, doublecomplex *, doublecomplex *, int *, doublecomplex *, int *);
    extern int ztrsen_(char *, char *, logical *, int *, doublecomplex *, int *, doublecomplex *, int *, doublecomplex *, int *, double *, double *, doublecomplex *, int *, int *);
    extern int zunm2r_(char *, char *, int *, int *, int *, doublecomplex *, int *, doublecomplex *, doublecomplex *, int *, doublecomplex *, int *);

#ifdef __cplusplus
}
#endif