#pragma once

#include "common.h"
#include "arcomp.h"

#ifdef __cplusplus
extern "C"
{
#endif

    // Method parameters:
    //
    // [in]       which = requested spectrum
    // [in]           k = number of eigenvalues
    // [in]         ncv = number of basis vectors
    // [in]       maxit = maximum number of iterations
    // [in]         tol = tolerance
    // [in]       sigma = sigma for shifted mode
    // [in]           A = csc matrix
    // [in,out]  result = eigenvalues/vectors storage

    // Complex standard problem (Hermitian or not) - regular mode
    EXPORT int ar_zi_ns(char* which, int k, int ncv, int maxit, double tol,
        ar_spmat* A, ar_result* result);

    // Complex standard problem (Hermitian or not) - shift-and-invert mode
    EXPORT int ar_zi_ns_shift(char* which, int k, int ncv, int maxit, double tol, arcomplex<double> sigma,
        ar_spmat* A, ar_result* result);

    // Complex generalized problem (Hermitian or not) - regular mode
    EXPORT int ar_zi_ng(char* which, int k, int ncv, int maxit, double tol,
        ar_spmat* A, ar_spmat* B, ar_result* result);

    // Complex generalized problem (Hermitian or not) - shift-and-invert mode
    EXPORT int ar_zi_ng_shift(char* which, int k, int ncv, int maxit, double tol, arcomplex<double> sigma,
        ar_spmat* A, ar_spmat* B, ar_result* result);

#ifdef __cplusplus
}
#endif