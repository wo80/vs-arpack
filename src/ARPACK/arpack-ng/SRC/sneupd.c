/* D:\Projekte\ARPACK\arpack-ng\SRC\sneupd.f -- translated by f2c (version 20100827). */

#include "arpack.h"


/**
 * \BeginDoc
 *
 * \Name: sneupd
 *
 * \Description:
 *
 *  This subroutine returns the converged approximations to eigenvalues
 *  of A*z = lambda*B*z and (optionally):
 *
 *      (1) The corresponding approximate eigenvectors;
 *
 *      (2) An orthonormal basis for the associated approximate
 *          invariant subspace;
 *
 *      (3) Both.
 *
 *  There is negligible additional cost to obtain eigenvectors.  An orthonormal
 *  basis is always computed.  There is an additional storage cost of n*nev
 *  if both are requested (in this case a separate array Z must be supplied).
 *
 *  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
 *  are derived from approximate eigenvalues and eigenvectors of
 *  of the linear operator OP prescribed by the MODE selection in the
 *  call to SNAUPD.  SNAUPD must be called before this routine is called.
 *  These approximate eigenvalues and vectors are commonly called Ritz
 *  values and Ritz vectors respectively.  They are referred to as such
 *  in the comments that follow.  The computed orthonormal basis for the
 *  invariant subspace corresponding to these Ritz values is referred to as a
 *  Schur basis.
 *
 *  See documentation in the header of the subroutine SNAUPD for
 *  definition of OP as well as other terms and the relation of computed
 *  Ritz values and Ritz vectors of OP with respect to the given problem
 *  A*z = lambda*B*z.  For a brief description, see definitions of
 *  IPARAM(7), MODE and WHICH in the documentation of SNAUPD.
 *
 * \Usage:
 *  call sneupd
 *     ( RVEC, HOWMNY, SELECT, DR, DI, Z, LDZ, SIGMAR, SIGMAI, WORKEV, BMAT,
 *       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL,
 *       LWORKL, INFO )
 *
 * \Arguments:
 *  RVEC    LOGICAL  (INPUT)
 *          Specifies whether a basis for the invariant subspace corresponding
 *          to the converged Ritz value approximations for the eigenproblem
 *          A*z = lambda*B*z is computed.
 *
 *             RVEC = .FALSE.     Compute Ritz values only.
 *
 *             RVEC = .TRUE.      Compute the Ritz vectors or Schur vectors.
 *                                See Remarks below.
 *
 *  HOWMNY  Character*1  (INPUT)
 *          Specifies the form of the basis for the invariant subspace
 *          corresponding to the converged Ritz values that is to be computed.
 *
 *          = 'A': Compute NEV Ritz vectors;
 *          = 'P': Compute NEV Schur vectors;
 *          = 'S': compute some of the Ritz vectors, specified
 *                 by the logical array SELECT.
 *
 *  SELECT  Logical array of dimension NCV.  (INPUT)
 *          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
 *          computed. To select the Ritz vector corresponding to a
 *          Ritz value (DR(j), DI(j)), SELECT(j) must be set to .TRUE..
 *          If HOWMNY = 'A' or 'P', SELECT is used as internal workspace.
 *
 *  DR      Real  array of dimension NEV+1.  (OUTPUT)
 *          If IPARAM(7) = 1,2 or 3 and SIGMAI=0.0  then on exit: DR contains
 *          the real part of the Ritz  approximations to the eigenvalues of
 *          A*z = lambda*B*z.
 *          If IPARAM(7) = 3, 4 and SIGMAI is not equal to zero, then on exit:
 *          DR contains the real part of the Ritz values of OP computed by
 *          SNAUPD. A further computation must be performed by the user
 *          to transform the Ritz values computed for OP by SNAUPD to those
 *          of the original system A*z = lambda*B*z. See remark 3 below.
 *
 *  DI      Real  array of dimension NEV+1.  (OUTPUT)
 *          On exit, DI contains the imaginary part of the Ritz value
 *          approximations to the eigenvalues of A*z = lambda*B*z associated
 *          with DR.
 *
 *          NOTE: When Ritz values are complex, they will come in complex
 *                conjugate pairs.  If eigenvectors are requested, the
 *                corresponding Ritz vectors will also come in conjugate
 *                pairs and the real and imaginary parts of these are
 *                represented in two consecutive columns of the array Z
 *                (see below).
 *
 *  Z       Real  N by NEV+1 array if RVEC = .TRUE. and HOWMNY = 'A'. (OUTPUT)
 *          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of
 *          Z represent approximate eigenvectors (Ritz vectors) corresponding
 *          to the NCONV=IPARAM(5) Ritz values for eigensystem
 *          A*z = lambda*B*z.
 *
 *          The complex Ritz vector associated with the Ritz value
 *          with positive imaginary part is stored in two consecutive
 *          columns.  The first column holds the real part of the Ritz
 *          vector and the second column holds the imaginary part.  The
 *          Ritz vector associated with the Ritz value with negative
 *          imaginary part is simply the complex conjugate of the Ritz vector
 *          associated with the positive imaginary part.
 *
 *          If  RVEC = .FALSE. or HOWMNY = 'P', then Z is not referenced.
 *
 *          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
 *          the array Z may be set equal to first NEV+1 columns of the Arnoldi
 *          basis array V computed by SNAUPD.  In this case the Arnoldi basis
 *          will be destroyed and overwritten with the eigenvector basis.
 *
 *  LDZ     Integer.  (INPUT)
 *          The leading dimension of the array Z.  If Ritz vectors are
 *          desired, then  LDZ >= max( 1, N ).  In any case,  LDZ >= 1.
 *
 *  SIGMAR  Real   (INPUT)
 *          If IPARAM(7) = 3 or 4, represents the real part of the shift.
 *          Not referenced if IPARAM(7) = 1 or 2.
 *
 *  SIGMAI  Real   (INPUT)
 *          If IPARAM(7) = 3 or 4, represents the imaginary part of the shift.
 *          Not referenced if IPARAM(7) = 1 or 2. See remark 3 below.
 *
 *  WORKEV  Real  work array of dimension 3*NCV.  (WORKSPACE)
 *
 *  **** The remaining arguments MUST be the same as for the   ****
 *  **** call to SNAUPD that was just completed.               ****
 *
 *  NOTE: The remaining arguments
 *
 *           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
 *           WORKD, WORKL, LWORKL, INFO
 *
 *         must be passed directly to SNEUPD following the last call
 *         to SNAUPD.  These arguments MUST NOT BE MODIFIED between
 *         the the last call to SNAUPD and the call to SNEUPD.
 *
 *  Three of these parameters (V, WORKL, INFO) are also output parameters:
 *
 *  V       Real  N by NCV array.  (INPUT/OUTPUT)
 *
 *          Upon INPUT: the NCV columns of V contain the Arnoldi basis
 *                      vectors for OP as constructed by SNAUPD .
 *
 *          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
 *                       contain approximate Schur vectors that span the
 *                       desired invariant subspace.  See Remark 2 below.
 *
 *          NOTE: If the array Z has been set equal to first NEV+1 columns
 *          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
 *          Arnoldi basis held by V has been overwritten by the desired
 *          Ritz vectors.  If a separate array Z has been passed then
 *          the first NCONV=IPARAM(5) columns of V will contain approximate
 *          Schur vectors that span the desired invariant subspace.
 *
 *  WORKL   Real  work array of length LWORKL.  (OUTPUT/WORKSPACE)
 *          WORKL(1:ncv*ncv+3*ncv) contains information obtained in
 *          snaupd.  They are not changed by sneupd.
 *          WORKL(ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) holds the
 *          real and imaginary part of the untransformed Ritz values,
 *          the upper quasi-triangular matrix for H, and the
 *          associated matrix representation of the invariant subspace for H.
 *
 *          Note: IPNTR(9:13) contains the pointer into WORKL for addresses
 *          of the above information computed by sneupd.
 *          -------------------------------------------------------------
 *          IPNTR(9):  pointer to the real part of the NCV RITZ values of the
 *                     original system.
 *          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of
 *                     the original system.
 *          IPNTR(11): pointer to the NCV corresponding error bounds.
 *          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
 *                     Schur matrix for H.
 *          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
 *                     of the upper Hessenberg matrix H. Only referenced by
 *                     sneupd if RVEC = .TRUE. See Remark 2 below.
 *          -------------------------------------------------------------
 *
 *  INFO    Integer.  (OUTPUT)
 *          Error flag on output.
 *
 *          =  0: Normal exit.
 *
 *          =  1: The Schur form computed by LAPACK routine slahqr
 *                could not be reordered by LAPACK routine strsen.
 *                Re-enter subroutine sneupd with IPARAM(5)=NCV and
 *                increase the size of the arrays DR and DI to have
 *                dimension at least dimension NCV and allocate at least NCV
 *                columns for Z. NOTE: Not necessary if Z and V share
 *                the same space. Please notify the authors if this error
 *                occurs.
 *
 *          = -1: N must be positive.
 *          = -2: NEV must be positive.
 *          = -3: NCV-NEV >= 2 and less than or equal to N.
 *          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
 *          = -6: BMAT must be one of 'I' or 'G'.
 *          = -7: Length of private work WORKL array is not sufficient.
 *          = -8: Error return from calculation of a real Schur form.
 *                Informational error from LAPACK routine slahqr.
 *          = -9: Error return from calculation of eigenvectors.
 *                Informational error from LAPACK routine strevc.
 *          = -10: IPARAM(7) must be 1,2,3,4.
 *          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
 *          = -12: HOWMNY = 'S' not yet implemented
 *          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
 *          = -14: SNAUPD did not find any eigenvalues to sufficient
 *                 accuracy.
 *          = -15: DNEUPD got a different count of the number of converged
 *                 Ritz values than DNAUPD got.  This indicates the user
 *                 probably made an error in passing data from DNAUPD to
 *                 DNEUPD or that the data was modified before entering
 *                 DNEUPD
 * \EndDoc */
/** \BeginLib
 *
 * \References:
 *  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
 *     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
 *     pp 357-385.
 *  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
 *     Restarted Arnoldi Iteration", Rice University Technical Report
 *     TR95-13, Department of Computational and Applied Mathematics.
 *  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
 *     Real Matrices", Linear Algebra and its Applications, vol 88/89,
 *     pp 575-595, (1987).
 *
 * \Routines called:
 *     ivout   ARPACK utility routine that prints integers.
 *     smout   ARPACK utility routine that prints matrices
 *     svout   ARPACK utility routine that prints vectors.
 *     sgeqr2  LAPACK routine that computes the QR factorization of
 *             a matrix.
 *     slacpy  LAPACK matrix copy routine.
 *     slahqr  LAPACK routine to compute the real Schur form of an
 *             upper Hessenberg matrix.
 *     slamch  LAPACK routine that determines machine constants.
 *     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     slaset  LAPACK matrix initialization routine.
 *     sorm2r  LAPACK routine that applies an orthogonal matrix in
 *             factored form.
 *     strevc  LAPACK routine to compute the eigenvectors of a matrix
 *             in upper quasi-triangular form.
 *     strsen  LAPACK routine that re-orders the Schur form.
 *     strmm   Level 3 BLAS matrix times an upper triangular matrix.
 *     sger    Level 2 BLAS rank one update to a matrix.
 *     scopy   Level 1 BLAS that copies one vector to another .
 *     sdot    Level 1 BLAS that computes the scalar product of two vectors.
 *     snrm2   Level 1 BLAS that computes the norm of a vector.
 *     sscal   Level 1 BLAS that scales a vector.
 *
 * \Remarks
 *
 *  1. Currently only HOWMNY = 'A' and 'P' are implemented.
 *
 *     Let trans(X) denote the transpose of X.
 *
 *  2. Schur vectors are an orthogonal representation for the basis of
 *     Ritz vectors. Thus, their numerical properties are often superior.
 *     If RVEC = .TRUE. then the relationship
 *             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
 *     trans(V(:,1:IPARAM(5))) * V(:,1:IPARAM(5)) = I are approximately
 *     satisfied. Here T is the leading submatrix of order IPARAM(5) of the
 *     real upper quasi-triangular matrix stored workl(ipntr(12)). That is,
 *     T is block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
 *     each 2-by-2 diagonal block has its diagonal elements equal and its
 *     off-diagonal elements of opposite sign.  Corresponding to each 2-by-2
 *     diagonal block is a complex conjugate pair of Ritz values. The real
 *     Ritz values are stored on the diagonal of T.
 *
 *  3. If IPARAM(7) = 3 or 4 and SIGMAI is not equal zero, then the user must
 *     form the IPARAM(5) Rayleigh quotients in order to transform the Ritz
 *     values computed by SNAUPD for OP to those of A*z = lambda*B*z.
 *     Set RVEC = .true. and HOWMNY = 'A', and
 *     compute
 *           trans(Z(:,I)) * A * Z(:,I) if DI(I) = 0.
 *     If DI(I) is not equal to zero and DI(I+1) = - D(I),
 *     then the desired real and imaginary parts of the Ritz value are
 *           trans(Z(:,I)) * A * Z(:,I) +  trans(Z(:,I+1)) * A * Z(:,I+1),
 *           trans(Z(:,I)) * A * Z(:,I+1) -  trans(Z(:,I+1)) * A * Z(:,I),
 *     respectively.
 *     Another possibility is to set RVEC = .true. and HOWMNY = 'P' and
 *     compute trans(V(:,1:IPARAM(5))) * A * V(:,1:IPARAM(5)) and then an upper
 *     quasi-triangular matrix of order IPARAM(5) is computed. See remark
 *     2 above.
 *
 * \Authors
 *     Danny Sorensen               Phuong Vu
 *     Richard Lehoucq              CRPC / Rice University
 *     Chao Yang                    Houston, Texas
 *     Dept. of Computational &
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * \SCCS Information: @(#)
 * FILE: neupd.F   SID: 2.7   DATE OF SID: 09/20/00   RELEASE: 2
 *
 * \EndLib
 */

/* Subroutine */ int sneupd_(logical *rvec, char *howmny, logical *select, 
	real *dr, real *di, real *z__, integer *ldz, real *sigmar, real *
	sigmai, real *workev, char *bmat, integer *n, char *which, integer *
	nev, real *tol, real *resid, integer *ncv, real *v, integer *ldv, 
	integer *iparam, integer *ipntr, real *workd, real *workl, integer *
	lworkl, integer *info)
{
    /* System generated locals */
    integer v_dim1, v_offset, z_dim1, z_offset, i__1;
    real r__1, r__2;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    integer j, k, ih, jj, np;
    real vl[1]	/* was [1][1] */;
    integer ibd, ldh, ldq, iri;
    real sep;
    integer irr, wri, wrr, mode;
    real eps23;
    extern /* Subroutine */ int sger_(integer *, integer *, real *, real *, 
	    integer *, real *, integer *, real *, integer *);
    integer ierr;
    real temp;
    integer iwev;
    char type__[6];
    real temp1;
    extern doublereal snrm2_(integer *, real *, integer *);
    integer ihbds, iconj;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    real conds;
    logical reord;
    extern /* Subroutine */ int sgemv_(char *, integer *, integer *, real *, 
	    real *, integer *, real *, integer *, real *, real *, integer *);
    integer nconv, iwork[1];
    real rnorm;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *);
    integer ritzi;
    extern /* Subroutine */ int strmm_(char *, char *, char *, char *, 
	    integer *, integer *, real *, real *, integer *, real *, integer *),
	    ivout_(integer *, integer *, 
	    integer *, integer *, char *, ftnlen), smout_(integer *, integer *
	    , integer *, real *, integer *, integer *, char *, ftnlen);
    integer ritzr;
    extern /* Subroutine */ int svout_(integer *, integer *, real *, integer *
	    , char *, ftnlen), sgeqr2_(integer *, integer *, real *, integer *
	    , real *, real *, integer *);
    integer nconv2;
    extern doublereal slapy2_(real *, real *);
    extern /* Subroutine */ int sorm2r_(char *, char *, integer *, integer *, 
	    integer *, real *, integer *, real *, real *, integer *, real *, 
	    integer *);
    integer iheigi, iheigr, bounds, invsub, iuptri, msglvl, outncv, ishift, 
	    numcnv;
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, real *, 
	    integer *, real *, integer *), slahqr_(logical *, logical 
	    *, integer *, integer *, integer *, real *, integer *, real *, 
	    real *, integer *, integer *, real *, integer *, integer *), 
	    slaset_(char *, integer *, integer *, real *, real *, real *, 
	    integer *), strevc_(char *, char *, logical *, integer *, 
	    real *, integer *, real *, integer *, real *, integer *, integer *
	    , integer *, real *, integer *), strsen_(char *, 
	    char *, logical *, integer *, real *, integer *, real *, integer *
	    , real *, real *, integer *, real *, real *, real *, integer *, 
	    integer *, integer *, integer *);
    extern doublereal slamch_(char *);
    extern /* Subroutine */ int sngets_(integer *, char *, integer *, integer 
	    *, real *, real *, real *, real *, real *);


/*     %----------------------------------------------------% */
/*     | Include files for debugging and timing information | */
/*     %----------------------------------------------------% */


/* \SCCS Information: @(#) */
/* FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 */

/*     %---------------------------------% */
/*     | See debug.doc for documentation | */
/*     %---------------------------------% */

/*     %------------------% */
/*     | Scalar Arguments | */
/*     %------------------% */

/*     %--------------------------------% */
/*     | See stat.doc for documentation | */
/*     %--------------------------------% */

/* \SCCS Information: @(#) */
/* FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 */



/*     %-----------------% */
/*     | Array Arguments | */
/*     %-----------------% */


/*     %------------% */
/*     | Parameters | */
/*     %------------% */


/*     %---------------% */
/*     | Local Scalars | */
/*     %---------------% */


/*     %----------------------% */
/*     | External Subroutines | */
/*     %----------------------% */


/*     %--------------------% */
/*     | External Functions | */
/*     %--------------------% */


/*     %---------------------% */
/*     | Intrinsic Functions | */
/*     %---------------------% */


/*     %-----------------------% */
/*     | Executable Statements | */
/*     %-----------------------% */

/*     %------------------------% */
/*     | Set default parameters | */
/*     %------------------------% */

    /* Parameter adjustments */
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --workd;
    --resid;
    --di;
    --dr;
    --workev;
    --select;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --iparam;
    --ipntr;
    --workl;

    /* Function Body */
    msglvl = debug_1.mneupd;
    mode = iparam[7];
    nconv = iparam[5];
    *info = 0;

/*     %---------------------------------% */
/*     | Get machine dependent constant. | */
/*     %---------------------------------% */

    eps23 = slamch_("Epsilon-Machine");
    d__1 = (doublereal) eps23;
    eps23 = pow_dd(&d__1, &d_23);

/*     %--------------% */
/*     | Quick return | */
/*     %--------------% */

    ierr = 0;

    if (nconv <= 0) {
	ierr = -14;
    } else if (*n <= 0) {
	ierr = -1;
    } else if (*nev <= 0) {
	ierr = -2;
    } else if (*ncv <= *nev + 1 || *ncv > *n) {
	ierr = -3;
    } else if (s_cmp(which, "LM", (ftnlen)2, (ftnlen)2) != 0 && s_cmp(which, 
	    "SM", (ftnlen)2, (ftnlen)2) != 0 && s_cmp(which, "LR", (ftnlen)2, 
	    (ftnlen)2) != 0 && s_cmp(which, "SR", (ftnlen)2, (ftnlen)2) != 0 
	    && s_cmp(which, "LI", (ftnlen)2, (ftnlen)2) != 0 && s_cmp(which, 
	    "SI", (ftnlen)2, (ftnlen)2) != 0) {
	ierr = -5;
    } else if (*(unsigned char *)bmat != 'I' && *(unsigned char *)bmat != 'G')
	     {
	ierr = -6;
    } else /* if(complicated condition) */ {
/* Computing 2nd power */
	i__1 = *ncv;
	if (*lworkl < i__1 * i__1 * 3 + *ncv * 6) {
	    ierr = -7;
	} else if (*(unsigned char *)howmny != 'A' && *(unsigned char *)
		howmny != 'P' && *(unsigned char *)howmny != 'S' && *rvec) {
	    ierr = -13;
	} else if (*(unsigned char *)howmny == 'S') {
	    ierr = -12;
	}
    }

    if (mode == 1 || mode == 2) {
	s_copy(type__, "REGULR", (ftnlen)6, (ftnlen)6);
    } else if (mode == 3 && *sigmai == 0.f) {
	s_copy(type__, "SHIFTI", (ftnlen)6, (ftnlen)6);
    } else if (mode == 3) {
	s_copy(type__, "REALPT", (ftnlen)6, (ftnlen)6);
    } else if (mode == 4) {
	s_copy(type__, "IMAGPT", (ftnlen)6, (ftnlen)6);
    } else {
	ierr = -10;
    }
    if (mode == 1 && *(unsigned char *)bmat == 'G') {
	ierr = -11;
    }

/*     %------------% */
/*     | Error Exit | */
/*     %------------% */

    if (ierr != 0) {
	*info = ierr;
	goto L9000;
    }

/*     %--------------------------------------------------------% */
/*     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q   | */
/*     | etc... and the remaining workspace.                    | */
/*     | Also update pointer to be used on output.              | */
/*     | Memory is laid out as follows:                         | */
/*     | workl(1:ncv*ncv) := generated Hessenberg matrix        | */
/*     | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary   | */
/*     |                                   parts of ritz values | */
/*     | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds   | */
/*     %--------------------------------------------------------% */

/*     %-----------------------------------------------------------% */
/*     | The following is used and set by SNEUPD.                  | */
/*     | workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed | */
/*     |                             real part of the Ritz values. | */
/*     | workl(ncv*ncv+4*ncv+1:ncv*ncv+5*ncv) := The untransformed | */
/*     |                        imaginary part of the Ritz values. | */
/*     | workl(ncv*ncv+5*ncv+1:ncv*ncv+6*ncv) := The untransformed | */
/*     |                           error bounds of the Ritz values | */
/*     | workl(ncv*ncv+6*ncv+1:2*ncv*ncv+6*ncv) := Holds the upper | */
/*     |                             quasi-triangular matrix for H | */
/*     | workl(2*ncv*ncv+6*ncv+1: 3*ncv*ncv+6*ncv) := Holds the    | */
/*     |       associated matrix representation of the invariant   | */
/*     |       subspace for H.                                     | */
/*     | GRAND total of NCV * ( 3 * NCV + 6 ) locations.           | */
/*     %-----------------------------------------------------------% */

    ih = ipntr[5];
    ritzr = ipntr[6];
    ritzi = ipntr[7];
    bounds = ipntr[8];
    ldh = *ncv;
    ldq = *ncv;
    iheigr = bounds + ldh;
    iheigi = iheigr + ldh;
    ihbds = iheigi + ldh;
    iuptri = ihbds + ldh;
    invsub = iuptri + ldh * *ncv;
    ipntr[9] = iheigr;
    ipntr[10] = iheigi;
    ipntr[11] = ihbds;
    ipntr[12] = iuptri;
    ipntr[13] = invsub;
    wrr = 1;
    wri = *ncv + 1;
    iwev = wri + *ncv;

/*     %-----------------------------------------% */
/*     | irr points to the REAL part of the Ritz | */
/*     |     values computed by _neigh before    | */
/*     |     exiting _naup2.                     | */
/*     | iri points to the IMAGINARY part of the | */
/*     |     Ritz values computed by _neigh      | */
/*     |     before exiting _naup2.              | */
/*     | ibd points to the Ritz estimates        | */
/*     |     computed by _neigh before exiting   | */
/*     |     _naup2.                             | */
/*     %-----------------------------------------% */

    irr = ipntr[14] + *ncv * *ncv;
    iri = irr + *ncv;
    ibd = iri + *ncv;

/*     %------------------------------------% */
/*     | RNORM is B-norm of the RESID(1:N). | */
/*     %------------------------------------% */

    rnorm = workl[ih + 2];
    workl[ih + 2] = 0.f;

    if (msglvl > 2) {
	svout_(&debug_1.logfil, ncv, &workl[irr], &debug_1.ndigit, "_neupd: "
		"Real part of Ritz values passed in from _NAUPD.", (ftnlen)55);
	svout_(&debug_1.logfil, ncv, &workl[iri], &debug_1.ndigit, "_neupd: "
		"Imag part of Ritz values passed in from _NAUPD.", (ftnlen)55);
	svout_(&debug_1.logfil, ncv, &workl[ibd], &debug_1.ndigit, "_neupd: "
		"Ritz estimates passed in from _NAUPD.", (ftnlen)45);
    }

    if (*rvec) {

	reord = FALSE_;

/*        %---------------------------------------------------% */
/*        | Use the temporary bounds array to store indices   | */
/*        | These will be used to mark the select array later | */
/*        %---------------------------------------------------% */

	i__1 = *ncv;
	for (j = 1; j <= i__1; ++j) {
	    workl[bounds + j - 1] = (real) j;
	    select[j] = FALSE_;
/* L10: */
	}

/*        %-------------------------------------% */
/*        | Select the wanted Ritz values.      | */
/*        | Sort the Ritz values so that the    | */
/*        | wanted ones appear at the tailing   | */
/*        | NEV positions of workl(irr) and     | */
/*        | workl(iri).  Move the corresponding | */
/*        | error estimates in workl(bound)     | */
/*        | accordingly.                        | */
/*        %-------------------------------------% */

	np = *ncv - *nev;
	ishift = 0;
	sngets_(&ishift, which, nev, &np, &workl[irr], &workl[iri], &workl[
		bounds], &workl[1], &workl[np + 1]);

	if (msglvl > 2) {
	    svout_(&debug_1.logfil, ncv, &workl[irr], &debug_1.ndigit, "_neu"
		    "pd: Real part of Ritz values after calling _NGETS.", (
		    ftnlen)54);
	    svout_(&debug_1.logfil, ncv, &workl[iri], &debug_1.ndigit, "_neu"
		    "pd: Imag part of Ritz values after calling _NGETS.", (
		    ftnlen)54);
	    svout_(&debug_1.logfil, ncv, &workl[bounds], &debug_1.ndigit, 
		    "_neupd: Ritz value indices after calling _NGETS.", (
		    ftnlen)48);
	}

/*        %-----------------------------------------------------% */
/*        | Record indices of the converged wanted Ritz values  | */
/*        | Mark the select array for possible reordering       | */
/*        %-----------------------------------------------------% */

	numcnv = 0;
	i__1 = *ncv;
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	    r__1 = eps23, r__2 = slapy2_(&workl[irr + *ncv - j], &workl[iri + 
		    *ncv - j]);
	    temp1 = dmax(r__1,r__2);
	    jj = workl[bounds + *ncv - j];
	    if (numcnv < nconv && workl[ibd + jj - 1] <= *tol * temp1) {
		select[jj] = TRUE_;
		++numcnv;
		if (jj > nconv) {
		    reord = TRUE_;
		}
	    }
/* L11: */
	}

/*        %-----------------------------------------------------------% */
/*        | Check the count (numcnv) of converged Ritz values with    | */
/*        | the number (nconv) reported by dnaupd.  If these two      | */
/*        | are different then there has probably been an error       | */
/*        | caused by incorrect passing of the dnaupd data.           | */
/*        %-----------------------------------------------------------% */

	if (msglvl > 2) {
	    ivout_(&debug_1.logfil, &c__1, &numcnv, &debug_1.ndigit, "_neupd"
		    ": Number of specified eigenvalues", (ftnlen)39);
	    ivout_(&debug_1.logfil, &c__1, &nconv, &debug_1.ndigit, "_neupd:"
		    " Number of \"converged\" eigenvalues", (ftnlen)41);
	}

	if (numcnv != nconv) {
	    *info = -15;
	    goto L9000;
	}

/*        %-----------------------------------------------------------% */
/*        | Call LAPACK routine slahqr to compute the real Schur form | */
/*        | of the upper Hessenberg matrix returned by SNAUPD.        | */
/*        | Make a copy of the upper Hessenberg matrix.               | */
/*        | Initialize the Schur vector matrix Q to the identity.     | */
/*        %-----------------------------------------------------------% */

	i__1 = ldh * *ncv;
	scopy_(&i__1, &workl[ih], &c__1, &workl[iuptri], &c__1);
	slaset_("All", ncv, ncv, &s_zero, &s_one, &workl[invsub], &ldq);
	slahqr_(&c_true, &c_true, ncv, &c__1, ncv, &workl[iuptri], &ldh, &
		workl[iheigr], &workl[iheigi], &c__1, ncv, &workl[invsub], &
		ldq, &ierr);
	scopy_(ncv, &workl[invsub + *ncv - 1], &ldq, &workl[ihbds], &c__1);

	if (ierr != 0) {
	    *info = -8;
	    goto L9000;
	}

	if (msglvl > 1) {
	    svout_(&debug_1.logfil, ncv, &workl[iheigr], &debug_1.ndigit, 
		    "_neupd: Real part of the eigenvalues of H", (ftnlen)41);
	    svout_(&debug_1.logfil, ncv, &workl[iheigi], &debug_1.ndigit, 
		    "_neupd: Imaginary part of the Eigenvalues of H", (ftnlen)
		    46);
	    svout_(&debug_1.logfil, ncv, &workl[ihbds], &debug_1.ndigit, 
		    "_neupd: Last row of the Schur vector matrix", (ftnlen)43)
		    ;
	    if (msglvl > 3) {
		smout_(&debug_1.logfil, ncv, ncv, &workl[iuptri], &ldh, &
			debug_1.ndigit, "_neupd: The upper quasi-triangular "
			"matrix ", (ftnlen)42);
	    }
	}

	if (reord) {

/*           %-----------------------------------------------------% */
/*           | Reorder the computed upper quasi-triangular matrix. | */
/*           %-----------------------------------------------------% */

	    strsen_("None", "V", &select[1], ncv, &workl[iuptri], &ldh, &
		    workl[invsub], &ldq, &workl[iheigr], &workl[iheigi], &
		    nconv2, &conds, &sep, &workl[ihbds], ncv, iwork, &c__1, &
		    ierr);

	    if (nconv2 < nconv) {
		nconv = nconv2;
	    }
	    if (ierr == 1) {
		*info = 1;
		goto L9000;
	    }

	    if (msglvl > 2) {
		svout_(&debug_1.logfil, ncv, &workl[iheigr], &debug_1.ndigit, 
			"_neupd: Real part of the eigenvalues of H--reordered"
			, (ftnlen)52);
		svout_(&debug_1.logfil, ncv, &workl[iheigi], &debug_1.ndigit, 
			"_neupd: Imag part of the eigenvalues of H--reordered"
			, (ftnlen)52);
		if (msglvl > 3) {
		    smout_(&debug_1.logfil, ncv, ncv, &workl[iuptri], &ldq, &
			    debug_1.ndigit, "_neupd: Quasi-triangular matrix"
			    " after re-ordering", (ftnlen)49);
		}
	    }

	}

/*        %---------------------------------------% */
/*        | Copy the last row of the Schur vector | */
/*        | into workl(ihbds).  This will be used | */
/*        | to compute the Ritz estimates of      | */
/*        | converged Ritz values.                | */
/*        %---------------------------------------% */

	scopy_(ncv, &workl[invsub + *ncv - 1], &ldq, &workl[ihbds], &c__1);

/*        %----------------------------------------------------% */
/*        | Place the computed eigenvalues of H into DR and DI | */
/*        | if a spectral transformation was not used.         | */
/*        %----------------------------------------------------% */

	if (s_cmp(type__, "REGULR", (ftnlen)6, (ftnlen)6) == 0) {
	    scopy_(&nconv, &workl[iheigr], &c__1, &dr[1], &c__1);
	    scopy_(&nconv, &workl[iheigi], &c__1, &di[1], &c__1);
	}

/*        %----------------------------------------------------------% */
/*        | Compute the QR factorization of the matrix representing  | */
/*        | the wanted invariant subspace located in the first NCONV | */
/*        | columns of workl(invsub,ldq).                            | */
/*        %----------------------------------------------------------% */

	sgeqr2_(ncv, &nconv, &workl[invsub], &ldq, &workev[1], &workev[*ncv + 
		1], &ierr);

/*        %---------------------------------------------------------% */
/*        | * Postmultiply V by Q using sorm2r.                     | */
/*        | * Copy the first NCONV columns of VQ into Z.            | */
/*        | * Postmultiply Z by R.                                  | */
/*        | The N by NCONV matrix Z is now a matrix representation  | */
/*        | of the approximate invariant subspace associated with   | */
/*        | the Ritz values in workl(iheigr) and workl(iheigi)      | */
/*        | The first NCONV columns of V are now approximate Schur  | */
/*        | vectors associated with the real upper quasi-triangular | */
/*        | matrix of order NCONV in workl(iuptri)                  | */
/*        %---------------------------------------------------------% */

	sorm2r_("Right", "Notranspose", n, ncv, &nconv, &workl[invsub], &ldq, 
		&workev[1], &v[v_offset], ldv, &workd[*n + 1], &ierr);
	slacpy_("All", n, &nconv, &v[v_offset], ldv, &z__[z_offset], ldz);

	i__1 = nconv;
	for (j = 1; j <= i__1; ++j) {

/*           %---------------------------------------------------% */
/*           | Perform both a column and row scaling if the      | */
/*           | diagonal element of workl(invsub,ldq) is negative | */
/*           | I'm lazy and don't take advantage of the upper    | */
/*           | quasi-triangular form of workl(iuptri,ldq)        | */
/*           | Note that since Q is orthogonal, R is a diagonal  | */
/*           | matrix consisting of plus or minus ones           | */
/*           %---------------------------------------------------% */

	    if (workl[invsub + (j - 1) * ldq + j - 1] < 0.f) {
		sscal_(&nconv, &s_m1, &workl[iuptri + j - 1], &ldq);
		sscal_(&nconv, &s_m1, &workl[iuptri + (j - 1) * ldq], &c__1);
	    }

/* L20: */
	}

	if (*(unsigned char *)howmny == 'A') {

/*           %--------------------------------------------% */
/*           | Compute the NCONV wanted eigenvectors of T | */
/*           | located in workl(iuptri,ldq).              | */
/*           %--------------------------------------------% */

	    i__1 = *ncv;
	    for (j = 1; j <= i__1; ++j) {
		if (j <= nconv) {
		    select[j] = TRUE_;
		} else {
		    select[j] = FALSE_;
		}
/* L30: */
	    }

	    strevc_("Right", "Select", &select[1], ncv, &workl[iuptri], &ldq, 
		    vl, &c__1, &workl[invsub], &ldq, ncv, &outncv, &workev[1],
		     &ierr);

	    if (ierr != 0) {
		*info = -9;
		goto L9000;
	    }

/*           %------------------------------------------------% */
/*           | Scale the returning eigenvectors so that their | */
/*           | Euclidean norms are all one. LAPACK subroutine | */
/*           | strevc returns each eigenvector normalized so  | */
/*           | that the element of largest magnitude has      | */
/*           | magnitude 1;                                   | */
/*           %------------------------------------------------% */

	    iconj = 0;
	    i__1 = nconv;
	    for (j = 1; j <= i__1; ++j) {

		if (workl[iheigi + j - 1] == 0.f) {

/*                 %----------------------% */
/*                 | real eigenvalue case | */
/*                 %----------------------% */

		    temp = snrm2_(ncv, &workl[invsub + (j - 1) * ldq], &c__1);
		    r__1 = 1.f / temp;
		    sscal_(ncv, &r__1, &workl[invsub + (j - 1) * ldq], &c__1);

		} else {

/*                 %-------------------------------------------% */
/*                 | Complex conjugate pair case. Note that    | */
/*                 | since the real and imaginary part of      | */
/*                 | the eigenvector are stored in consecutive | */
/*                 | columns, we further normalize by the      | */
/*                 | square root of two.                       | */
/*                 %-------------------------------------------% */

		    if (iconj == 0) {
			r__1 = snrm2_(ncv, &workl[invsub + (j - 1) * ldq], &
				c__1);
			r__2 = snrm2_(ncv, &workl[invsub + j * ldq], &c__1);
			temp = slapy2_(&r__1, &r__2);
			r__1 = 1.f / temp;
			sscal_(ncv, &r__1, &workl[invsub + (j - 1) * ldq], &
				c__1);
			r__1 = 1.f / temp;
			sscal_(ncv, &r__1, &workl[invsub + j * ldq], &c__1);
			iconj = 1;
		    } else {
			iconj = 0;
		    }

		}

/* L40: */
	    }

	    sgemv_("T", ncv, &nconv, &s_one, &workl[invsub], &ldq, &workl[
		    ihbds], &c__1, &s_zero, &workev[1], &c__1);

	    iconj = 0;
	    i__1 = nconv;
	    for (j = 1; j <= i__1; ++j) {
		if (workl[iheigi + j - 1] != 0.f) {

/*                 %-------------------------------------------% */
/*                 | Complex conjugate pair case. Note that    | */
/*                 | since the real and imaginary part of      | */
/*                 | the eigenvector are stored in consecutive | */
/*                 %-------------------------------------------% */

		    if (iconj == 0) {
			workev[j] = slapy2_(&workev[j], &workev[j + 1]);
			workev[j + 1] = workev[j];
			iconj = 1;
		    } else {
			iconj = 0;
		    }
		}
/* L45: */
	    }

	    if (msglvl > 2) {
		scopy_(ncv, &workl[invsub + *ncv - 1], &ldq, &workl[ihbds], &
			c__1);
		svout_(&debug_1.logfil, ncv, &workl[ihbds], &debug_1.ndigit, 
			"_neupd: Last row of the eigenvector matrix for T", (
			ftnlen)48);
		if (msglvl > 3) {
		    smout_(&debug_1.logfil, ncv, ncv, &workl[invsub], &ldq, &
			    debug_1.ndigit, "_neupd: The eigenvector matrix "
			    "for T", (ftnlen)36);
		}
	    }

/*           %---------------------------------------% */
/*           | Copy Ritz estimates into workl(ihbds) | */
/*           %---------------------------------------% */

	    scopy_(&nconv, &workev[1], &c__1, &workl[ihbds], &c__1);

/*           %---------------------------------------------------------% */
/*           | Compute the QR factorization of the eigenvector matrix  | */
/*           | associated with leading portion of T in the first NCONV | */
/*           | columns of workl(invsub,ldq).                           | */
/*           %---------------------------------------------------------% */

	    sgeqr2_(ncv, &nconv, &workl[invsub], &ldq, &workev[1], &workev[*
		    ncv + 1], &ierr);

/*           %----------------------------------------------% */
/*           | * Postmultiply Z by Q.                       | */
/*           | * Postmultiply Z by R.                       | */
/*           | The N by NCONV matrix Z is now contains the  | */
/*           | Ritz vectors associated with the Ritz values | */
/*           | in workl(iheigr) and workl(iheigi).          | */
/*           %----------------------------------------------% */

	    sorm2r_("Right", "Notranspose", n, ncv, &nconv, &workl[invsub], &
		    ldq, &workev[1], &z__[z_offset], ldz, &workd[*n + 1], &
		    ierr, (ftnlen)5, (ftnlen)11);

	    strmm_("Right", "Upper", "No transpose", "Non-unit", n, &nconv, &
		    s_one, &workl[invsub], &ldq, &z__[z_offset], ldz, (ftnlen)
		    5, (ftnlen)5, (ftnlen)12, (ftnlen)8);

	}

    } else {

/*        %------------------------------------------------------% */
/*        | An approximate invariant subspace is not needed.     | */
/*        | Place the Ritz values computed SNAUPD into DR and DI | */
/*        %------------------------------------------------------% */

	scopy_(&nconv, &workl[ritzr], &c__1, &dr[1], &c__1);
	scopy_(&nconv, &workl[ritzi], &c__1, &di[1], &c__1);
	scopy_(&nconv, &workl[ritzr], &c__1, &workl[iheigr], &c__1);
	scopy_(&nconv, &workl[ritzi], &c__1, &workl[iheigi], &c__1);
	scopy_(&nconv, &workl[bounds], &c__1, &workl[ihbds], &c__1);
    }

/*     %------------------------------------------------% */
/*     | Transform the Ritz values and possibly vectors | */
/*     | and corresponding error bounds of OP to those  | */
/*     | of A*x = lambda*B*x.                           | */
/*     %------------------------------------------------% */

    if (s_cmp(type__, "REGULR", (ftnlen)6, (ftnlen)6) == 0) {

	if (*rvec) {
	    sscal_(ncv, &rnorm, &workl[ihbds], &c__1);
	}

    } else {

/*        %---------------------------------------% */
/*        |   A spectral transformation was used. | */
/*        | * Determine the Ritz estimates of the | */
/*        |   Ritz values in the original system. | */
/*        %---------------------------------------% */

	if (s_cmp(type__, "SHIFTI", (ftnlen)6, (ftnlen)6) == 0) {

	    if (*rvec) {
		sscal_(ncv, &rnorm, &workl[ihbds], &c__1);
	    }

	    i__1 = *ncv;
	    for (k = 1; k <= i__1; ++k) {
		temp = slapy2_(&workl[iheigr + k - 1], &workl[iheigi + k - 1])
			;
		workl[ihbds + k - 1] = (r__1 = workl[ihbds + k - 1], dabs(
			r__1)) / temp / temp;
/* L50: */
	    }

	} else if (s_cmp(type__, "REALPT", (ftnlen)6, (ftnlen)6) == 0) {

	    i__1 = *ncv;
	    for (k = 1; k <= i__1; ++k) {
/* L60: */
	    }

	} else if (s_cmp(type__, "IMAGPT", (ftnlen)6, (ftnlen)6) == 0) {

	    i__1 = *ncv;
	    for (k = 1; k <= i__1; ++k) {
/* L70: */
	    }

	}

/*        %-----------------------------------------------------------% */
/*        | *  Transform the Ritz values back to the original system. | */
/*        |    For TYPE = 'SHIFTI' the transformation is              | */
/*        |             lambda = 1/theta + sigma                      | */
/*        |    For TYPE = 'REALPT' or 'IMAGPT' the user must from     | */
/*        |    Rayleigh quotients or a projection. See remark 3 above.| */
/*        | NOTES:                                                    | */
/*        | *The Ritz vectors are not affected by the transformation. | */
/*        %-----------------------------------------------------------% */

	if (s_cmp(type__, "SHIFTI", (ftnlen)6, (ftnlen)6) == 0) {

	    i__1 = *ncv;
	    for (k = 1; k <= i__1; ++k) {
		temp = slapy2_(&workl[iheigr + k - 1], &workl[iheigi + k - 1])
			;
		workl[iheigr + k - 1] = workl[iheigr + k - 1] / temp / temp + 
			*sigmar;
		workl[iheigi + k - 1] = -workl[iheigi + k - 1] / temp / temp 
			+ *sigmai;
/* L80: */
	    }

	    scopy_(&nconv, &workl[iheigr], &c__1, &dr[1], &c__1);
	    scopy_(&nconv, &workl[iheigi], &c__1, &di[1], &c__1);

	} else if (s_cmp(type__, "REALPT", (ftnlen)6, (ftnlen)6) == 0 || 
		s_cmp(type__, "IMAGPT", (ftnlen)6, (ftnlen)6) == 0) {

	    scopy_(&nconv, &workl[iheigr], &c__1, &dr[1], &c__1);
	    scopy_(&nconv, &workl[iheigi], &c__1, &di[1], &c__1);

	}

    }

    if (s_cmp(type__, "SHIFTI", (ftnlen)6, (ftnlen)6) == 0 && msglvl > 1) {
	svout_(&debug_1.logfil, &nconv, &dr[1], &debug_1.ndigit, "_neupd: Un"
		"transformed real part of the Ritz valuess.", (ftnlen)52);
	svout_(&debug_1.logfil, &nconv, &di[1], &debug_1.ndigit, "_neupd: Un"
		"transformed imag part of the Ritz valuess.", (ftnlen)52);
	svout_(&debug_1.logfil, &nconv, &workl[ihbds], &debug_1.ndigit, "_ne"
		"upd: Ritz estimates of untransformed Ritz values.", (ftnlen)
		52);
    } else if (s_cmp(type__, "REGULR", (ftnlen)6, (ftnlen)6) == 0 && msglvl > 
	    1) {
	svout_(&debug_1.logfil, &nconv, &dr[1], &debug_1.ndigit, "_neupd: Re"
		"al parts of converged Ritz values.", (ftnlen)44);
	svout_(&debug_1.logfil, &nconv, &di[1], &debug_1.ndigit, "_neupd: Im"
		"ag parts of converged Ritz values.", (ftnlen)44);
	svout_(&debug_1.logfil, &nconv, &workl[ihbds], &debug_1.ndigit, "_ne"
		"upd: Associated Ritz estimates.", (ftnlen)34);
    }

/*     %-------------------------------------------------% */
/*     | Eigenvector Purification step. Formally perform | */
/*     | one of inverse subspace iteration. Only used    | */
/*     | for MODE = 2.                                   | */
/*     %-------------------------------------------------% */

    if (*rvec && *(unsigned char *)howmny == 'A' && s_cmp(type__, "SHIFTI", (
	    ftnlen)6, (ftnlen)6) == 0) {

/*        %------------------------------------------------% */
/*        | Purify the computed Ritz vectors by adding a   | */
/*        | little bit of the residual vector:             | */
/*        |                      T                         | */
/*        |          resid(:)*( e    s ) / theta           | */
/*        |                      NCV                       | */
/*        | where H s = s theta. Remember that when theta  | */
/*        | has nonzero imaginary part, the corresponding  | */
/*        | Ritz vector is stored across two columns of Z. | */
/*        %------------------------------------------------% */

	iconj = 0;
	i__1 = nconv;
	for (j = 1; j <= i__1; ++j) {
	    if (workl[iheigi + j - 1] == 0.f && workl[iheigr + j - 1] != 0.f) 
		    {
		workev[j] = workl[invsub + (j - 1) * ldq + *ncv - 1] / workl[
			iheigr + j - 1];
	    } else if (iconj == 0) {
		temp = slapy2_(&workl[iheigr + j - 1], &workl[iheigi + j - 1])
			;
		if (temp != 0.f) {
		    workev[j] = (workl[invsub + (j - 1) * ldq + *ncv - 1] * 
			    workl[iheigr + j - 1] + workl[invsub + j * ldq + *
			    ncv - 1] * workl[iheigi + j - 1]) / temp / temp;
		    workev[j + 1] = (workl[invsub + j * ldq + *ncv - 1] * 
			    workl[iheigr + j - 1] - workl[invsub + (j - 1) * 
			    ldq + *ncv - 1] * workl[iheigi + j - 1]) / temp / 
			    temp;
		}
		iconj = 1;
	    } else {
		iconj = 0;
	    }
/* L110: */
	}

/*        %---------------------------------------% */
/*        | Perform a rank one update to Z and    | */
/*        | purify all the Ritz vectors together. | */
/*        %---------------------------------------% */

	sger_(n, &nconv, &s_one, &resid[1], &c__1, &workev[1], &c__1, &z__[
		z_offset], ldz);

    }

L9000:

    return 0;

/*     %---------------% */
/*     | End of SNEUPD | */
/*     %---------------% */

} /* sneupd_ */

