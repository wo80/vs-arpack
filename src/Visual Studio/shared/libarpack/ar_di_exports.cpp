#pragma warning (disable : 4250)

#include "ar_di_exports.h"

#include "arlsmat.h"
#include "arlnsmat.h"
#include "arlssym.h"
#include "arlgsym.h"
#include "arlsnsym.h"
#include "arlgnsym.h"
#include "arlscomp.h"
#include "arlgcomp.h"

#define AR_TYPE double

int ar_di_ss(char* which, int k, int ncv, int maxit, double tol,
		   ar_spmat *A, ar_result *result)
{
	char uplo = 'U';
	int nconv = 0; // Number of converged eigenvalues.
	
	AR_TYPE *eigvec = (AR_TYPE *)result->eigvec;
	AR_TYPE *eigvalr = (AR_TYPE *)result->eigvalr;

	try
	{
		// Defining the eigenvalue problem.
		ARluSymMatrix<AR_TYPE> matrix(A->n, A->nnz, (AR_TYPE *)A->x, A->i, A->p, uplo);
		ARluSymStdEig<AR_TYPE> prob(k, matrix, which, ncv, tol, maxit);

		// Finding eigenvalues.
		if (eigvec == NULL)
		{
			nconv = prob.Eigenvalues(eigvalr);
		}
		else
		{
			nconv = prob.EigenValVectors(eigvec, eigvalr, 0);
		}

		result->iterations = prob.GetIter();
		result->ncv = prob.GetNcv();
		result->info = 0;
	}
	catch(const ArpackError& e)
	{
		result->info = e.Status();
	}

	return nconv;
}

int ar_di_ss_shift(char* which,  int k, int ncv, int maxit, double tol, double sigma,
				 ar_spmat *A, ar_result *result)
{
	char uplo = 'U';
	int nconv = 0; // Number of converged eigenvalues.

	AR_TYPE *eigvec = (AR_TYPE *)result->eigvec;
	AR_TYPE *eigvalr = (AR_TYPE *)result->eigvalr;

	try
	{
		// Defining the eigenvalue problem.
		ARluSymMatrix<AR_TYPE> matrix(A->n, A->nnz, (AR_TYPE *)A->x, A->i, A->p, uplo);
		ARluSymStdEig<AR_TYPE> prob(k, matrix, sigma, which, ncv, tol, maxit);

		// Finding eigenvalues.
		if (eigvec == NULL)
		{
			nconv = prob.Eigenvalues(eigvalr);
		}
		else
		{
			nconv = prob.EigenValVectors(eigvec, eigvalr, 0);
		}

		result->iterations = prob.GetIter();
		result->ncv = prob.GetNcv();
		result->info = 0;
	}
	catch(const ArpackError& e)
	{
		result->info = e.Status();
	}

	return nconv;
}

int ar_di_sg(char* which, int k, int ncv, int maxit, double tol,
	ar_spmat *A, ar_spmat *B, ar_result *result)
{
	char uplo = 'U';
	int nconv = 0; // Number of converged eigenvalues.

	AR_TYPE *eigvec = (AR_TYPE *)result->eigvec;
	AR_TYPE *eigvalr = (AR_TYPE *)result->eigvalr;

	try
	{
		// Defining the eigenvalue problem.
		ARluSymMatrix<AR_TYPE> matrixA(A->n, A->nnz, (AR_TYPE *)A->x, A->i, A->p, uplo);
		ARluSymMatrix<AR_TYPE> matrixB(B->n, B->nnz, (AR_TYPE *)B->x, B->i, B->p, uplo);
		ARluSymGenEig<AR_TYPE> prob(k, matrixA, matrixB, which, ncv, tol, maxit);

		// Finding eigenvalues.
		if (eigvec == NULL)
		{
			nconv = prob.Eigenvalues(eigvalr);
		}
		else
		{
			nconv = prob.EigenValVectors(eigvec, eigvalr, 0);
		}

		result->iterations = prob.GetIter();
		result->ncv = prob.GetNcv();
		result->info = 0;
	}
	catch (const ArpackError& e)
	{
		result->info = e.Status();
	}

	return nconv;
}

int ar_di_sg_shift(char* which, char mode, int k, int ncv, int maxit, double tol, double sigma,
	ar_spmat *A, ar_spmat *B, ar_result *result)
{
	char uplo = 'U';
	int nconv = 0; // Number of converged eigenvalues.

	AR_TYPE *eigvec = (AR_TYPE *)result->eigvec;
	AR_TYPE *eigvalr = (AR_TYPE *)result->eigvalr;

	try
	{
		// Defining the eigenvalue problem.
		ARluSymMatrix<AR_TYPE> matrixA(A->n, A->nnz, (AR_TYPE *)A->x, A->i, A->p, uplo);
		ARluSymMatrix<AR_TYPE> matrixB(B->n, B->nnz, (AR_TYPE *)B->x, B->i, B->p, uplo);
		ARluSymGenEig<AR_TYPE> prob(mode, k, matrixA, matrixB, sigma, which, ncv, tol, maxit);

		// Finding eigenvalues.
		if (eigvec == NULL)
		{
			nconv = prob.Eigenvalues(eigvalr);
		}
		else
		{
			nconv = prob.EigenValVectors(eigvec, eigvalr, 0);
		}

		result->iterations = prob.GetIter();
		result->ncv = prob.GetNcv();
		result->info = 0;
	}
	catch (const ArpackError& e)
	{
		result->info = e.Status();
	}

	return nconv;
}

int ar_di_ns(char* which, int k, int ncv, int maxit, double tol,
		   ar_spmat *A, ar_result *result)
{
	int nconv = 0; // Number of converged eigenvalues.
	
	AR_TYPE *eigvec = (AR_TYPE *)result->eigvec;
	AR_TYPE *eigvalr = (AR_TYPE *)result->eigvalr;
	AR_TYPE *eigvali = (AR_TYPE *)result->eigvali;

	try
	{
		ARluNonSymMatrix<AR_TYPE, double> matrix(A->n, A->nnz, (AR_TYPE *)A->x, A->i, A->p);
		ARluNonSymStdEig<AR_TYPE> prob(k, matrix, which, ncv, tol, maxit);

		if (eigvec == NULL)
		{
			nconv = prob.Eigenvalues(eigvalr, eigvali);
		}
		else
		{
			nconv = prob.EigenValVectors(eigvec, eigvalr, eigvali, 0);
		}

		result->iterations = prob.GetIter();
		result->ncv = prob.GetNcv();
		result->info = 0;
	}
	catch(const ArpackError& e)
	{
		result->info = e.Status();
	}

	return nconv;
}

int ar_di_ns_shift(char* which,  int k, int ncv, int maxit, double tol, double sigma,
				 ar_spmat *A, ar_result *result)
{
	int nconv = 0; // Number of converged eigenvalues.
	
	AR_TYPE *eigvec = (AR_TYPE *)result->eigvec;
	AR_TYPE *eigvalr = (AR_TYPE *)result->eigvalr;
	AR_TYPE *eigvali = (AR_TYPE *)result->eigvali;

	try
	{
		ARluNonSymMatrix<AR_TYPE, double> matrix(A->n, A->nnz, (AR_TYPE *)A->x, A->i, A->p);
		ARluNonSymStdEig<AR_TYPE> prob(k, matrix, sigma, which, ncv, tol, maxit);

		if (eigvec == NULL)
		{
			nconv = prob.Eigenvalues(eigvalr, eigvali);
		}
		else
		{
			nconv = prob.EigenValVectors(eigvec, eigvalr, eigvali, 0);
		}

		result->iterations = prob.GetIter();
		result->ncv = prob.GetNcv();
		result->info = 0;
	}
	catch(const ArpackError& e)
	{
		result->info = e.Status();
	}

	return nconv;
}

int ar_di_ng(char* which, int k, int ncv, int maxit, double tol,
		   ar_spmat *A, ar_spmat *B, ar_result *result)
{
	int nconv = 0; // Number of converged eigenvalues.
	
	AR_TYPE *eigvec = (AR_TYPE *)result->eigvec;
	AR_TYPE *eigvalr = (AR_TYPE *)result->eigvalr;
	AR_TYPE *eigvali = (AR_TYPE *)result->eigvali;

	try
	{
		ARluNonSymMatrix<AR_TYPE, double> matrixA(A->n, A->nnz, (AR_TYPE *)A->x, A->i, A->p);
		ARluNonSymMatrix<AR_TYPE, double> matrixB(B->n, B->nnz, (AR_TYPE *)B->x, B->i, B->p);
		ARluNonSymGenEig<AR_TYPE> prob(k, matrixA, matrixB, which, ncv, tol, maxit);

		if (eigvec == NULL)
		{
			nconv = prob.Eigenvalues(eigvalr, eigvali);
		}
		else
		{
			nconv = prob.EigenValVectors(eigvec, eigvalr, eigvali, 0);
		}

		result->iterations = prob.GetIter();
		result->ncv = prob.GetNcv();
		result->info = 0;
	}
	catch(const ArpackError& e)
	{
		result->info = e.Status();
	}

	return nconv;
}

int ar_di_ng_shift(char* which,  int k, int ncv, int maxit, double tol, double sigma,
				 ar_spmat *A, ar_spmat *B, ar_result *result)
{
	int nconv = 0; // Number of converged eigenvalues.
	
	AR_TYPE *eigvec = (AR_TYPE *)result->eigvec;
	AR_TYPE *eigvalr = (AR_TYPE *)result->eigvalr;
	AR_TYPE *eigvali = (AR_TYPE *)result->eigvali;

	try
	{
		ARluNonSymMatrix<AR_TYPE, double> matrixA(A->n, A->nnz, (AR_TYPE *)A->x, A->i, A->p);
		ARluNonSymMatrix<AR_TYPE, double> matrixB(B->n, B->nnz, (AR_TYPE *)B->x, B->i, B->p);
		ARluNonSymGenEig<AR_TYPE> prob(k, matrixA, matrixB, sigma, which, ncv, tol, maxit);

		if (eigvec == NULL)
		{
			nconv = prob.Eigenvalues(eigvalr, eigvali);
		}
		else
		{
			nconv = prob.EigenValVectors(eigvec, eigvalr, eigvali, 0);
		}

		result->iterations = prob.GetIter();
		result->ncv = prob.GetNcv();
		result->info = 0;
	}
	catch(const ArpackError& e)
	{
		result->info = e.Status();
	}

	return nconv;
}

int ar_di_ng_shift_cx(char* which,  int k, int ncv, int maxit, double tol, char part, double sigma_r, double sigma_i,
				 ar_spmat *A, ar_spmat *B, ar_result *result)
{
	int nconv = 0; // Number of converged eigenvalues.
	
	AR_TYPE *eigvec = (AR_TYPE *)result->eigvec;
	AR_TYPE *eigvalr = (AR_TYPE *)result->eigvalr;
	AR_TYPE *eigvali = (AR_TYPE *)result->eigvali;

	try
	{
		ARluNonSymMatrix<AR_TYPE, double> matrixA(A->n, A->nnz, (AR_TYPE *)A->x, A->i, A->p);
		ARluNonSymMatrix<AR_TYPE, double> matrixB(B->n, B->nnz, (AR_TYPE *)B->x, B->i, B->p);
		ARluNonSymGenEig<AR_TYPE> prob(k, matrixA, matrixB, part, sigma_r, sigma_i, which, ncv, tol, maxit);

		if (eigvec == NULL)
		{
			nconv = prob.Eigenvalues(eigvalr, eigvali);
		}
		else
		{
			nconv = prob.EigenValVectors(eigvec, eigvalr, eigvali, 0);
		}

		result->iterations = prob.GetIter();
		result->ncv = prob.GetNcv();
		result->info = 0;
	}
	catch(const ArpackError& e)
	{
		result->info = e.Status();
	}

	return nconv;
}

#undef AR_TYPE
