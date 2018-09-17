#pragma warning (disable : 4250)

#include "ar_zi_exports.h"

#include "arlsmat.h"
#include "arlnsmat.h"
#include "arlssym.h"
#include "arlgsym.h"
#include "arlsnsym.h"
#include "arlgnsym.h"
#include "arlscomp.h"
#include "arlgcomp.h"

#define AR_TYPE arcomplex<double>

int ar_zi_ns(char* which, int k, int ncv, int maxit, double tol,
		   ar_spmat *A, ar_result *result)
{
	int nconv = 0; // Number of converged eigenvalues.
	
	AR_TYPE *eigvec = (AR_TYPE *)result->eigvec;
	AR_TYPE *eigval = (AR_TYPE *)result->eigvalr;

	try
	{
		ARluNonSymMatrix<AR_TYPE, double> matrix(A->n, A->nnz, (AR_TYPE *)A->x, A->i, A->p);
		ARluCompStdEig<double> prob(k, matrix, which, ncv, tol, maxit);

		if (eigvec == NULL)
		{
			nconv = prob.Eigenvalues(eigval);
		}
		else
		{
			nconv = prob.EigenValVectors(eigvec, eigval);
		}

		result->iterations = prob.GetIter();
		result->info = 0;
	}
	catch(const ArpackError& e)
	{
		result->info = e.Status();
	}

	return nconv;
}

int ar_zi_ns_shift(char* which,  int k, int ncv, int maxit, double tol, arcomplex<double> sigma,
				 ar_spmat *A, ar_result *result)
{
	int nconv = 0; // Number of converged eigenvalues.
	
	AR_TYPE *eigvec = (AR_TYPE *)result->eigvec;
	AR_TYPE *eigval = (AR_TYPE *)result->eigvalr;

	try
	{
		ARluNonSymMatrix<AR_TYPE, double> matrix(A->n, A->nnz, (AR_TYPE *)A->x, A->i, A->p);
		ARluCompStdEig<double> prob(k, matrix, sigma, which, ncv, tol, maxit);

		if (eigvec == NULL)
		{
			nconv = prob.Eigenvalues(eigval);
		}
		else
		{
			nconv = prob.EigenValVectors(eigvec, eigval);
		}

		result->iterations = prob.GetIter();
		result->info = 0;
	}
	catch(const ArpackError& e)
	{
		result->info = e.Status();
	}

	return nconv;
}

int ar_zi_ng(char* which, int k, int ncv, int maxit, double tol,
		   ar_spmat *A, ar_spmat *B, ar_result *result)
{
	int nconv = 0; // Number of converged eigenvalues.
	
	AR_TYPE *eigvec = (AR_TYPE *)result->eigvec;
	AR_TYPE *eigval = (AR_TYPE *)result->eigvalr;

	try
	{
		ARluNonSymMatrix<AR_TYPE, double> matrixA(A->n, A->nnz, (AR_TYPE *)A->x, A->i, A->p);
		ARluNonSymMatrix<AR_TYPE, double> matrixB(B->n, B->nnz, (AR_TYPE *)B->x, B->i, B->p);
		ARluCompGenEig<double> prob(k, matrixA, matrixB, which, ncv, tol, maxit);

		if (eigvec == NULL)
		{
			nconv = prob.Eigenvalues(eigval);
		}
		else
		{
			nconv = prob.EigenValVectors(eigvec, eigval);
		}

		result->iterations = prob.GetIter();
		result->info = 0;
	}
	catch(const ArpackError& e)
	{
		result->info = e.Status();
	}

	return nconv;
}

int ar_zi_ng_shift(char* which,  int k, int ncv, int maxit, double tol, arcomplex<double> sigma,
				 ar_spmat *A, ar_spmat *B, ar_result *result)
{
	int nconv = 0; // Number of converged eigenvalues.
	
	AR_TYPE *eigvec = (AR_TYPE *)result->eigvec;
	AR_TYPE *eigval = (AR_TYPE *)result->eigvalr;

	try
	{
		ARluNonSymMatrix<AR_TYPE, double> matrixA(A->n, A->nnz, (AR_TYPE *)A->x, A->i, A->p);
		ARluNonSymMatrix<AR_TYPE, double> matrixB(B->n, B->nnz, (AR_TYPE *)B->x, B->i, B->p);
		ARluCompGenEig<double> prob(k, matrixA, matrixB, sigma, which, ncv, tol, maxit);

		if (eigvec == NULL)
		{
			nconv = prob.Eigenvalues(eigval);
		}
		else
		{
			nconv = prob.EigenValVectors(eigvec, eigval);
		}

		result->iterations = prob.GetIter();
		result->info = 0;
	}
	catch(const ArpackError& e)
	{
		result->info = e.Status();
	}

	return nconv;
}

#undef AR_TYPE
