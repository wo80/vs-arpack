#pragma once

#include "common.h"
#include "arcomp.h"

#define EXPORT __declspec(dllexport)

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
	EXPORT int ar_ci_ns(char* which, int k, int ncv, int maxit, double tol,
		spmat *A, ar_result *result);

	// Complex standard problem (Hermitian or not) - shift-and-invert mode
	EXPORT int ar_ci_ns_shift(char* which,  int k, int ncv, int maxit, double tol, arcomplex<float> sigma,
		spmat *A, ar_result *result);

	// Complex generalized problem (Hermitian or not) - regular mode
	EXPORT int ar_ci_ng(char* which, int k, int ncv, int maxit, double tol,
		spmat *A, spmat *B, ar_result *result);

	// Complex generalized problem (Hermitian or not) - shift-and-invert mode
	EXPORT int ar_ci_ng_shift(char* which,  int k, int ncv, int maxit, double tol, arcomplex<float> sigma,
		spmat *A, spmat *B, ar_result *result);

#ifdef __cplusplus
}
#endif