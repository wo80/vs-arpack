#pragma once

#include "f2c.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /* LAPACK (complex) */

    extern void cgeqr2_(int*, int*, a_fcomplex*, int*, a_fcomplex*, a_fcomplex*, int*);
    extern void clacpy_(char*, int*, int*, a_fcomplex*, int*, a_fcomplex*, int*);
    extern void clahqr_(logical*, logical*, int*, int*, int*, a_fcomplex*, int*, a_fcomplex*, int*, int*, a_fcomplex*, int*, int*);
    extern float clanhs_(char*, int*, a_fcomplex*, int*, a_fcomplex*);
    extern void clarnv_(int*, int*, int*, a_fcomplex*);
    extern void clartg_(a_fcomplex*, a_fcomplex*, float*, a_fcomplex*, a_fcomplex*);
    extern void clascl_(char*, int*, int*, float*, float*, int*, int*, a_fcomplex*, int*, int*);
    extern void claset_(char*, int*, int*, a_fcomplex*, a_fcomplex*, a_fcomplex*, int*);
    extern void ctrevc_(char*, char*, logical*, int*, a_fcomplex*, int*, a_fcomplex*, int*, a_fcomplex*, int*, int*, int*, a_fcomplex*, float*, int*);
    extern void ctrsen_(char*, char*, logical*, int*, a_fcomplex*, int*, a_fcomplex*, int*, a_fcomplex*, int*, float*, float*, a_fcomplex*, int*, int*);
    extern void cunm2r_(char*, char*, int*, int*, int*, a_fcomplex*, int*, a_fcomplex*, a_fcomplex*, int*, a_fcomplex*, int*);
    extern void cgttrf_(int *, a_fcomplex *, a_fcomplex *, a_fcomplex *, a_fcomplex *, int *, int *);
    extern void cgttrs_(char *, int *, int *, a_fcomplex *, a_fcomplex *, a_fcomplex *, a_fcomplex *, int *, a_fcomplex *, int *, int *);

    /* LAPACK (double) */

    extern void dorm2r_(char*, char*, int*, int*, int*, double*, int*, double*, double*, int*, double*, int*);
    extern void dgeqr2_(int*, int*, double*, int*, double*, double*, int*);
    extern void dlabad_(double*, double*);
    extern void dlacpy_(char*, int*, int*, double*, int*, double*, int*);
    extern void dlae2_(double*, double*, double*, double*, double*);
    extern void dlaev2_(double*, double*, double*, double*, double*, double*, double*);
    extern void dlahqr_(logical*, logical*, int*, int*, int*, double*, int*, double*, double*, int*, int*, double*, int*, int*);
    extern double dlamch_(char*);
    extern double dlanhs_(char*, int*, double*, int*, double*);
    extern double dlanst_(char*, int*, double*, double*);
    extern double dlapy2_(double*, double*);
    extern void dlarf_(char*, int*, int*, double*, int*, double*, double*, int*, double*);
    extern void dlarfg_(int*, double*, double*, int*, double*);
    extern void dlarnv_(int*, int*, int*, double*);
    extern void dlartg_(double*, double*, double*, double*, double*);
    extern void dlascl_(char*, int*, int*, double*, double*, int*, int*, double*, int*, int*);
    extern void dlaset_(char*, int*, int*, double*, double*, double*, int*);
    extern void dlasr_(char*, char*, char*, int*, int*, double*, double*, double*, int*);
    extern void dlasrt_(char*, int*, double*, int*);
    extern void dsteqr_(char*, int*, double*, double*, double*, int*, double*, int*);
    extern void dtrevc_(char*, char*, logical*, int*, double*, int*, double*, int*, double*, int*, int*, int*, double*, int*);
    extern void dtrsen_(char*, char*, logical*, int*, double*, int*, double*, int*, double*, double*, int*, double*, double*, double*, int*, int*, int*, int*);
    extern void dgttrf_(int *, double *, double *, double *, double *, int *, int *);
    extern void dgttrs_(char *, int *, int *, double *, double *, double *, double *, int *, double *, int *, int *);
    extern void dpttrf_(int *, double *, double *, int *);
    extern void dpttrs_(int *, int *, double *, double *, double *, int *, int *);

    /* LAPACK (single) */

    extern void sgeqr2_(int*, int*, float*, int*, float*, float*, int*);
    extern void slabad_(float*, float*);
    extern void slacpy_(char*, int*, int*, float*, int*, float*, int*);
    extern void slae2_(float*, float*, float*, float*, float*);
    extern void slaev2_(float*, float*, float*, float*, float*, float*, float*);
    extern void slahqr_(logical*, logical*, int*, int*, int*, float*, int*, float*, float*, int*, int*, float*, int*, int*);
    extern float slamch_(char*);
    extern float slanhs_(char*, int*, float*, int*, float*);
    extern float slanst_(char*, int*, float*, float*);
    extern float slapy2_(float*, float*);
    extern void slarf_(char*, int*, int*, float*, int*, float*, float*, int*, float*);
    extern void slarfg_(int*, float*, float*, int*, float*);
    extern void slarnv_(int*, int*, int*, float*);
    extern void slartg_(float*, float*, float*, float*, float*);
    extern void slascl_(char*, int*, int*, float*, float*, int*, int*, float*, int*, int*);
    extern void slaset_(char*, int*, int*, float*, float*, float*, int*);
    extern void slasr_(char*, char*, char*, int*, int*, float*, float*, float*, int*);
    extern void slasrt_(char*, int*, float*, int*);
    extern void sorm2r_(char*, char*, int*, int*, int*, float*, int*, float*, float*, int*, float*, int*);
    extern void ssteqr_(char*, int*, float*, float*, float*, int*, float*, int*);
    extern void strevc_(char*, char*, logical*, int*, float*, int*, float*, int*, float*, int*, int*, int*, float*, int*);
    extern void strsen_(char*, char*, logical*, int*, float*, int*, float*, int*, float*, float*, int*, float*, float*, float*, int*, int*, int*, int*);
    extern void sgttrf_(int *, float *, float *, float *, float *, int *, int *);
    extern void sgttrs_(char *, int *, int *, float *, float *, float *, float *, int *, float *, int *, int *);
    extern void spttrf_(int *, float *, float *, int *);
    extern void spttrs_(int *, int *, float *, float *, float *, int *, int *);

    /* LAPACK (doublecomplex) */

    extern void zdscal_(int*, double*, a_dcomplex*, int*);
    extern void zgeqr2_(int*, int*, a_dcomplex*, int*, a_dcomplex*, a_dcomplex*, int*);
    extern void zlacpy_(char*, int*, int*, a_dcomplex*, int*, a_dcomplex*, int*);
    extern void zlahqr_(logical*, logical*, int*, int*, int*, a_dcomplex*, int*, a_dcomplex*, int*, int*, a_dcomplex*, int*, int*);
    extern double zlanhs_(char*, int*, a_dcomplex*, int*, a_dcomplex*);
    extern void zlarnv_(int*, int*, int*, a_dcomplex*);
    extern void zlartg_(a_dcomplex*, a_dcomplex*, double*, a_dcomplex*, a_dcomplex*);
    extern void zlascl_(char*, int*, int*, double*, double*, int*, int*, a_dcomplex*, int*, int*);
    extern void zlaset_(char*, int*, int*, a_dcomplex*, a_dcomplex*, a_dcomplex*, int*);
    extern void ztrevc_(char*, char*, logical*, int*, a_dcomplex*, int*, a_dcomplex*, int*, a_dcomplex*, int*, int*, int*, a_dcomplex*, double*, int*);
    extern void ztrsen_(char*, char*, logical*, int*, a_dcomplex*, int*, a_dcomplex*, int*, a_dcomplex*, int*, double*, double*, a_dcomplex*, int*, int*);
    extern void zunm2r_(char*, char*, int*, int*, int*, a_dcomplex*, int*, a_dcomplex*, a_dcomplex*, int*, a_dcomplex*, int*);
    extern void zgttrf_(int *, a_dcomplex *,  a_dcomplex *, a_dcomplex *, a_dcomplex *, int *,  int *);
    extern void zgttrs_(char *, int *, int *, a_dcomplex *, a_dcomplex *, a_dcomplex *, a_dcomplex *, int *, a_dcomplex *, int *, int *);

#ifdef __cplusplus
}
#endif
