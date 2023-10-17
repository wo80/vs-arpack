#pragma once

#ifdef _WIN32
# ifdef WIN_IMPORT
#   define EXPORT __declspec( dllimport )
# else
#   define EXPORT __declspec( dllexport )
# endif
#else
# define EXPORT
#endif

/* Sparse matrix in column compressed format. */
typedef struct ar_spmat_t {
	/* Number of rows. */
    int  m;
	/* Number of columns. */
    int  n;
	/* Array of nonzero values. */
    void *x;
	/* Array of column indices. */
    int  *i;
	/* Array of row pointers. */
    int  *p;
	/* Number of nonzeros in the matrix. */
    int  nnz;
} ar_spmat;

typedef struct ar_result_t {
	/* Eigenvalues. */
    void *eigvalr;
	/* Eigenvalues (imaginary part) - only used for real problem with complex eigenvalues. */
    void *eigvali;
	/* Eigenvectors. */
    void *eigvec;
	/* Number of iterations taken. */
    int  iterations;
	/* Number of Arnoldi vectors. */
    int  ncv;
	/* Error info. */
    int  info;
} ar_result;
