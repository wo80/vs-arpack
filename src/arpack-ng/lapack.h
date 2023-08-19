#pragma once

#include "f2c.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /* complex */

    extern doublereal scnrm2_(integer *, complex *, integer *);
    extern int caxpy_(integer *, complex *, complex *, integer *, complex *, integer *);
    extern int ccopy_(integer *, complex *, integer *, complex *, integer *);
    extern int cgemv_(char *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *);
    extern int cgeqr2_(integer *, integer *, complex *, integer *, complex *, complex *, integer *);
    extern int cgeru_(integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, integer *);
    extern int clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *);
    extern int clahqr_(logical *, logical *, integer *, integer *, integer *, complex *, integer *, complex *, integer *, integer *, complex *, integer *, integer *);
    extern doublereal clanhs_(char *, integer *, complex *, integer *, complex *);
    extern int clarnv_(integer *, integer *, integer *, complex *);
    extern int clartg_(complex *, complex *, real *, complex *, complex *);
    extern int clascl_(char *, integer *, integer *, real *, real *, integer *, integer *, complex *, integer *, integer *);
    extern int claset_(char *, integer *, integer *, complex *, complex *, complex *, integer *);
    extern /* Complex */ VOID cdotc_(complex *, integer *, complex *, integer *, complex *, integer *);
    extern int cscal_(integer *, complex *, complex *, integer *);
    extern int csscal_(integer *, real *, complex *, integer *);
    extern int ctrevc_(char *, char *, logical *, integer *, complex *, integer *, complex *, integer *, complex *, integer *, integer *, integer *, complex *, real *, integer *);
    extern int ctrmm_(char *, char *, char *, char *, integer *, integer *, complex *, complex *, integer *, complex *, integer *);
    extern int ctrsen_(char *, char *, logical *, integer *, complex *, integer *, complex *, integer *, complex *, integer *, real *, real *, complex *, integer *, integer *);
    extern int cunm2r_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *, complex *, integer *, complex *, integer *);

    /* double */

    extern int daxpy_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, integer *);
    extern int dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
    extern int dgeqr2_(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
    extern int dger_(integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *);
    extern int dlabad_(doublereal *, doublereal *);
    extern int dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *);
    extern int dlae2_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    extern int dlaev2_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    extern int dlahqr_(logical *, logical *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *);
    extern doublereal dlanhs_(char *, integer *, doublereal *, integer *, doublereal *);
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *);
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern int dlarf_(char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *);
    extern int dlarfg_(integer *, doublereal *, doublereal *, integer *, doublereal *);
    extern int dlarnv_(integer *, integer *, integer *, doublereal *);
    extern int dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    extern int dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *, doublereal *, integer *, integer *);
    extern int dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *);
    extern int dlasr_(char *, char *, char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *);
    extern int dlasrt_(char *, integer *, doublereal *, integer *);
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern int dorm2r_(char *, char *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    extern int dsteqr_(char *, integer *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *);
    extern int dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    extern int dtrevc_(char *, char *, logical *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, integer *, integer *, doublereal *, integer *);
    extern int dtrmm_(char *, char *, char *, char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
    extern int dtrsen_(char *, char *, logical *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *, doublereal *, doublereal *, integer *, integer *, integer *, integer *);

    /* single */

    extern int saxpy_(integer *, real *, real *, integer *, real *, integer *);
    extern int scopy_(integer *, real *, integer *, real *, integer *);
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    extern int sgemv_(char *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *);
    extern int sgeqr2_(integer *, integer *, real *, integer *, real *, real *, integer *);
    extern int sger_(integer *, integer *, real *, real *, integer *, real *, integer *, real *, integer *);
    extern int slabad_(real *, real *);
    extern int slacpy_(char *, integer *, integer *, real *, integer *, real *, integer *);
    extern int slae2_(real *, real *, real *, real *, real *);
    extern int slaev2_(real *, real *, real *, real *, real *, real *, real *);
    extern int slahqr_(logical *, logical *, integer *, integer *, integer *, real *, integer *, real *, real *, integer *, integer *, real *, integer *, integer *);
    extern doublereal slamch_(char *);
    extern doublereal slanhs_(char *, integer *, real *, integer *, real *);
    extern doublereal slanst_(char *, integer *, real *, real *);
    extern doublereal slapy2_(real *, real *);
    extern int slarf_(char *, integer *, integer *, real *, integer *, real *, real *, integer *, real *);
    extern int slarfg_(integer *, real *, real *, integer *, real *);
    extern int slarnv_(integer *, integer *, integer *, real *);
    extern int slartg_(real *, real *, real *, real *, real *);
    extern int slascl_(char *, integer *, integer *, real *, real *, integer *, integer *, real *, integer *, integer *);
    extern int slaset_(char *, integer *, integer *, real *, real *, real *, integer *);
    extern int slasr_(char *, char *, char *, integer *, integer *, real *, real *, real *, integer *);
    extern int slasrt_(char *, integer *, real *, integer *);
    extern doublereal snrm2_(integer *, real *, integer *);
    extern int sorm2r_(char *, char *, integer *, integer *, integer *, real *, integer *, real *, real *, integer *, real *, integer *);
    extern int sscal_(integer *, real *, real *, integer *);
    extern int ssteqr_(char *, integer *, real *, real *, real *, integer *, real *, integer *);
    extern int sswap_(integer *, real *, integer *, real *, integer *);
    extern int strevc_(char *, char *, logical *, integer *, real *, integer *, real *, integer *, real *, integer *, integer *, integer *, real *, integer *);
    extern int strmm_(char *, char *, char *, char *, integer *, integer *, real *, real *, integer *, real *, integer *);
    extern int strsen_(char *, char *, logical *, integer *, real *, integer *, real *, integer *, real *, real *, integer *, real *, real *, real *, integer *, integer *, integer *, integer *);

    /* doublecomplex */

    extern doublereal dznrm2_(integer *, doublecomplex *, integer *);
    extern int zaxpy_(integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *);
    extern int zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    extern int zdscal_(integer *, doublereal *, doublecomplex *, integer *);
    extern int zgemv_(char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *);
    extern int zgeqr2_(integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *);
    extern int zgeru_(integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    extern int zlacpy_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    extern int zlahqr_(logical *, logical *, integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, integer *, doublecomplex *, integer *, integer *);
    extern doublereal zlanhs_(char *, integer *, doublecomplex *, integer *, doublecomplex *);
    extern int zlarnv_(integer *, integer *, integer *, doublecomplex *);
    extern int zlartg_(doublecomplex *, doublecomplex *, doublereal *, doublecomplex *, doublecomplex *);
    extern int zlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *, doublecomplex *, integer *, integer *);
    extern int zlaset_(char *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *, integer *);
    extern int zscal_(integer *, doublecomplex *, doublecomplex *, integer *);
    extern int ztrevc_(char *, char *, logical *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, integer *, integer *, doublecomplex *, doublereal *, integer *);
    extern int ztrmm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *);
    extern int ztrsen_(char *, char *, logical *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublereal *, doublereal *, doublecomplex *, integer *, integer *);
    extern int zunm2r_(char *, char *, integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *);

#ifdef __cplusplus
}
#endif