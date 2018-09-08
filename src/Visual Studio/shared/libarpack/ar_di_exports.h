#pragma once

#include "common.h"

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

	// Real symmetric standard problem - regular mode
	EXPORT int ar_di_ss(char* which, int k, int ncv, int maxit, double tol,
		spmat *A, ar_result *result);

	// Real symmetric standard problem - shift-and-invert mode
	EXPORT int ar_di_ss_shift(char* which,  int k, int ncv, int maxit, double tol, double sigma,
		spmat *A, ar_result *result);

	// Real symmetric generalized problem - regular mode
	EXPORT int ar_di_sg(char* which, int k, int ncv, int maxit, double tol,
		spmat *A, spmat *B, ar_result *result);

	// Real symmetric generalized problem - shift-and-invert mode (standard, buckling or Caley)
	EXPORT int ar_di_sg_shift(char* which, char mode, int k, int ncv, int maxit, double tol, double sigma,
		spmat *A, spmat *B, ar_result *result);

	// Real non-symmetric standard problem - regular mode
	EXPORT int ar_di_ns(char* which, int k, int ncv, int maxit, double tol,
		spmat *A, ar_result *result);

	// Real non-symmetric standard problem - shift-and-invert mode
	EXPORT int ar_di_ns_shift(char* which,  int k, int ncv, int maxit, double tol, double sigma,
		spmat *A, ar_result *result);

	// Real non-symmetric generalized problem - regular mode
	EXPORT int ar_di_ng(char* which, int k, int ncv, int maxit, double tol,
		spmat *A, spmat *B, ar_result *result);

	// Real non-symmetric generalized problem - shift-and-invert mode
	EXPORT int ar_di_ng_shift(char* which,  int k, int ncv, int maxit, double tol, double sigma,
		spmat *A, spmat *B, ar_result *result);

	// Real non-symmetric generalized problem - complex shift-and-invert mode
	EXPORT int ar_di_ng_shift_cx(char* which,  int k, int ncv, int maxit, double tol, char part, double sigma_r, double sigma_i,
		spmat *A, spmat *B, ar_result *result);

#ifdef __cplusplus
}
#endif