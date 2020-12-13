/* arpack-ng\SRC\zgetv0.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \BeginDoc
 *
 * \Name: zgetv0
 *
 * \Description:
 *  Generate a random initial residual vector for the Arnoldi process.
 *  Force the residual vector to be in the range of the operator OP.
 *
 * \Usage:
 *  call zgetv0
 *     ( IDO, BMAT, ITRY, INITV, N, J, V, LDV, RESID, RNORM,
 *       IPNTR, WORKD, IERR )
 *
 * \Arguments
 *  IDO     Integer.  (INPUT/OUTPUT)
 *          Reverse communication flag.  IDO must be zero on the first
 *          call to zgetv0.
 *          -------------------------------------------------------------
 *          IDO =  0: first call to the reverse communication interface
 *          IDO = -1: compute  Y = OP * X  where
 *                    IPNTR(1) is the pointer into WORKD for X,
 *                    IPNTR(2) is the pointer into WORKD for Y.
 *                    This is for the initialization phase to force the
 *                    starting vector into the range of OP.
 *          IDO =  2: compute  Y = B * X  where
 *                    IPNTR(1) is the pointer into WORKD for X,
 *                    IPNTR(2) is the pointer into WORKD for Y.
 *          IDO = 99: done
 *          -------------------------------------------------------------
 *
 *  BMAT    Character*1.  (INPUT)
 *          BMAT specifies the type of the matrix B in the (generalized)
 *          eigenvalue problem A*x = lambda*B*x.
 *          B = 'I' -> standard eigenvalue problem A*x = lambda*x
 *          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
 *
 *  ITRY    Integer.  (INPUT)
 *          ITRY counts the number of times that zgetv0 is called.
 *          It should be set to 1 on the initial call to zgetv0.
 *
 *  INITV   Logical variable.  (INPUT)
 *          .TRUE.  => the initial residual vector is given in RESID.
 *          .FALSE. => generate a random initial residual vector.
 *
 *  N       Integer.  (INPUT)
 *          Dimension of the problem.
 *
 *  J       Integer.  (INPUT)
 *          Index of the residual vector to be generated, with respect to
 *          the Arnoldi process.  J > 1 in case of a "restart".
 *
 *  V       Complex*16 N by J array.  (INPUT)
 *          The first J-1 columns of V contain the current Arnoldi basis
 *          if this is a "restart".
 *
 *  LDV     Integer.  (INPUT)
 *          Leading dimension of V exactly as declared in the calling
 *          program.
 *
 *  RESID   Complex*16 array of length N.  (INPUT/OUTPUT)
 *          Initial residual vector to be generated.  If RESID is
 *          provided, force RESID into the range of the operator OP.
 *
 *  RNORM   Double precision scalar.  (OUTPUT)
 *          B-norm of the generated residual.
 *
 *  IPNTR   Integer array of length 3.  (OUTPUT)
 *
 *  WORKD   Complex*16 work array of length 2*N.  (REVERSE COMMUNICATION).
 *          On exit, WORK(1:N) = B*RESID to be used in SSAITR.
 *
 *  IERR    Integer.  (OUTPUT)
 *          =  0: Normal exit.
 *          = -1: Cannot generate a nontrivial restarted residual vector
 *                in the range of the operator OP.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Local variables:
 *     xxxxxx  Complex*16
 *
 * \References:
 *  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
 *     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
 *     pp 357-385.
 *
 * \Routines called:
 *     arscnd  ARPACK utility routine for timing.
 *     zvout   ARPACK utility routine that prints vectors.
 *     zlarnv  LAPACK routine for generating a random vector.
 *     zgemv   Level 2 BLAS routine for matrix vector multiplication.
 *     zcopy   Level 1 BLAS that copies one vector to another.
 *     zdotc   Level 1 BLAS that computes the scalar product of two vectors.
 *     dznrm2  Level 1 BLAS that computes the norm of a vector.
 *
 * \Author
 *     Danny Sorensen               Phuong Vu
 *     Richard Lehoucq              CRPC / Rice University
 *     Dept. of Computational &     Houston, Texas
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * \SCCS Information: @(#)
 * FILE: getv0.F   SID: 2.3   DATE OF SID: 08/27/96   RELEASE: 2
 *
 * \EndLib
 */
int zgetv0_(int *ido, char *bmat, int *itry, bool *initv, int *n, int *j,
            zomplex *v, int *ldv, zomplex *resid, double *rnorm, int *ipntr, zomplex *workd,
            int *ierr)
{
    /* Initialized data */

    static bool inits = true;

    /* System generated locals */
    int v_dim1, v_offset, i__1, i__2;
    double d__1, d__2;
    zomplex z__1;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static float t0, t1, t2, t3;
    int jj;
    static int iter;
    static bool orth;
    static int iseed[4];
    int idist;
    zomplex cnorm;
    static bool first;
    static double rnorm0;
    static int msglvl;

    /* Parameter adjustments */
    --workd;
    --resid;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --ipntr;

    /* Function Body */


    /* --------------------------------- */
    /* Initialize the seed of the LAPACK */
    /* random number generator           */
    /* --------------------------------- */

    if (inits)
    {
        iseed[0] = 1;
        iseed[1] = 3;
        iseed[2] = 5;
        iseed[3] = 7;
        inits = false;
    }

    if (*ido == 0)
    {
        /* ----------------------------- */
        /* Initialize timing statistics  */
        /* & message level for debugging */
        /* ----------------------------- */

#ifndef NO_TIMER
        arscnd_(&t0);
#endif

        msglvl = debug_1.mgetv0;

        *ierr = 0;
        iter = 0;
        first = false;
        orth = false;

        /* --------------------------------------------------- */
        /* Possibly generate a random starting vector in RESID */
        /* Use a LAPACK random number generator used by the    */
        /* matrix generation routines.                         */
        /*    idist = 1: uniform (0,1)  distribution;          */
        /*    idist = 2: uniform (-1,1) distribution;          */
        /*    idist = 3: normal  (0,1)  distribution;          */
        /* --------------------------------------------------- */

        if (! (*initv))
        {
            idist = 2;
            zlarnv_(&idist, iseed, n, &resid[1]);
        }

        /* -------------------------------------------------------- */
        /* Force the starting vector into the range of OP to handle */
        /* the generalized problem when B is possibly (singular).   */
        /* -------------------------------------------------------- */

#ifndef NO_TIMER
        arscnd_(&t2);
#endif

        if (*itry == 1)
        {
            ++timing_1.nopx;
            ipntr[1] = 1;
            ipntr[2] = *n + 1;
            zcopy_(n, &resid[1], &c__1, &workd[1], &c__1);
            *ido = -1;
            goto L9000;
        }
        else if (*itry > 1 && *bmat == 'G')
        {
            zcopy_(n, &resid[1], &c__1, &workd[*n + 1], &c__1);
        }
    }

    /* --------------------------------------- */
    /* Back from computing OP*(initial-vector) */
    /* --------------------------------------- */

    if (first)
    {
        goto L20;
    }

    /* ---------------------------------------------- */
    /* Back from computing OP*(orthogonalized-vector) */
    /* ---------------------------------------------- */

    if (orth)
    {
        goto L40;
    }

#ifndef NO_TIMER
    arscnd_(&t3);
    timing_1.tmvopx += t3 - t2;
#endif

    /* ---------------------------------------------------- */
    /* Starting vector is now in the range of OP; r = OP*r; */
    /* Compute B-norm of starting vector.                   */
    /* ---------------------------------------------------- */

#ifndef NO_TIMER
    arscnd_(&t2);
#endif

    first = true;
    if (*itry == 1)
    {
        zcopy_(n, &workd[*n + 1], &c__1, &resid[1], &c__1);
    }
    if (*bmat == 'G')
    {
        ++timing_1.nbx;
        ipntr[1] = *n + 1;
        ipntr[2] = 1;
        *ido = 2;
        goto L9000;
    }
    else if (*bmat == 'I')
    {
        zcopy_(n, &resid[1], &c__1, &workd[1], &c__1);
    }

L20:

#ifndef NO_TIMER
    if (*bmat == 'G')
    {
        arscnd_(&t3);
        timing_1.tmvbx += t3 - t2;
    }
#endif

    first = false;
    if (*bmat == 'G')
    {
        zdotc_(&z__1, n, &resid[1], &c__1, &workd[1], &c__1);
        cnorm.r = z__1.r, cnorm.i = z__1.i;
        d__1 = cnorm.r;
        d__2 = cnorm.i;
        rnorm0 = sqrt(dlapy2_(&d__1, &d__2));
    }
    else if (*bmat == 'I')
    {
        rnorm0 = dznrm2_(n, &resid[1], &c__1);
    }
    *rnorm = rnorm0;

    /* ------------------------------------------- */
    /* Exit if this is the very first Arnoldi step */
    /* ------------------------------------------- */

    if (*j == 1)
    {
        goto L50;
    }

    /* ------------------------------------------------------------- */
    /* Otherwise need to B-orthogonalize the starting vector against */
    /* the current Arnoldi basis using Gram-Schmidt with iter. ref.  */
    /* This is the case where an invariant subspace is encountered   */
    /* in the middle of the Arnoldi factorization.                   */
    /*                                                               */
    /*       s = V^{T}*B*r;   r = r - V*s;                           */
    /*                                                               */
    /* Stopping criteria used for iter. ref. is discussed in         */
    /* Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   */
    /* ------------------------------------------------------------- */

    orth = true;
L30:

    i__1 = *j - 1;
    zgemv_("C", n, &i__1, &z_one, &v[v_offset], ldv, &workd[1], &c__1, &z_zero, &workd[*n + 1], &c__1);
    i__1 = *j - 1;
    z__1.r = -1., z__1.i = -0.0;
    zgemv_("N", n, &i__1, &z__1, &v[v_offset], ldv, &workd[*n + 1], &c__1, &z_one, &resid[1], &c__1);

    /* -------------------------------------------------------- */
    /* Compute the B-norm of the orthogonalized starting vector */
    /* -------------------------------------------------------- */

#ifndef NO_TIMER
    arscnd_(&t2);
#endif

    if (*bmat == 'G')
    {
        ++timing_1.nbx;
        zcopy_(n, &resid[1], &c__1, &workd[*n + 1], &c__1);
        ipntr[1] = *n + 1;
        ipntr[2] = 1;
        *ido = 2;
        goto L9000;
    }
    else if (*bmat == 'I')
    {
        zcopy_(n, &resid[1], &c__1, &workd[1], &c__1);
    }

L40:

#ifndef NO_TIMER
    if (*bmat == 'G')
    {
        arscnd_(&t3);
        timing_1.tmvbx += t3 - t2;
    }
#endif

    if (*bmat == 'G')
    {
        zdotc_(&z__1, n, &resid[1], &c__1, &workd[1], &c__1);
        cnorm.r = z__1.r, cnorm.i = z__1.i;
        d__1 = cnorm.r;
        d__2 = cnorm.i;
        *rnorm = sqrt(dlapy2_(&d__1, &d__2));
    }
    else if (*bmat == 'I')
    {
        *rnorm = dznrm2_(n, &resid[1], &c__1);
    }

    /* ------------------------------------ */
    /* Check for further orthogonalization. */
    /* ------------------------------------ */

#ifndef NO_TRACE
    if (msglvl > 2)
    {
        dvout_(&c__1, &rnorm0, &debug_1.ndigit, "_getv0: re-orthonalization ; rnorm0 is");
        dvout_(&c__1, rnorm, &debug_1.ndigit, "_getv0: re-orthonalization ; rnorm is");
    }
#endif

    if (*rnorm > rnorm0 * 0.717)
    {
        goto L50;
    }

    ++iter;
    if (iter <= 1)
    {
        /* --------------------------------- */
        /* Perform iterative refinement step */
        /* --------------------------------- */

        rnorm0 = *rnorm;
        goto L30;
    }
    else
    {
        /* ---------------------------------- */
        /* Iterative refinement step "failed" */
        /* ---------------------------------- */

        i__1 = *n;
        for (jj = 1; jj <= i__1; ++jj)
        {
            i__2 = jj;
            resid[i__2].r = 0.0, resid[i__2].i = 0.0;

        }
        *rnorm = 0.0;
        *ierr = -1;
    }

L50:

#ifndef NO_TRACE
    if (msglvl > 0)
    {
        dvout_(&c__1, rnorm, &debug_1.ndigit, "_getv0: B-norm of initial / restarted starting vector");
    }

    if (msglvl > 2)
    {
        zvout_(n, &resid[1], &debug_1.ndigit, "_getv0: initial / restarted starting vector");
    }
#endif

    *ido = 99;

#ifndef NO_TIMER
    arscnd_(&t1);
    timing_1.tgetv0 += t1 - t0;
#endif

L9000:
    return 0;

    /* ------------- */
    /* End of zgetv0 */
    /* ------------- */

} /* zgetv0_ */

