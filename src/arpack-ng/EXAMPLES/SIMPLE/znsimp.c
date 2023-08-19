/* EXAMPLES\SIMPLE\znsimp.f -- translated by f2c (version 20100827). */

#include <stdlib.h>
#include "arpack_internal.h"

int znsimp_av_(const int nx, zomplex* v, zomplex* w);
int znsimp_tv_(const int nx, zomplex* x, zomplex* y);

/**
 * \BeginDoc
 *
/**
 * \BeginDoc
 *
 *     This example program is intended to illustrate the
 *     simplest case of using ARPACK in considerable detail.
 *     This code may be used to understand basic usage of ARPACK
 *     and as a template for creating an interface to ARPACK.
 *
 *     This code shows how to use ARPACK to find a few eigenvalues
 *     (lambda) and corresponding eigenvectors (x) for the standard
 *     eigenvalue problem:
 *
 *                        A*x = lambda*x
 *
 *     where A is a general n by n complex matrix.
 *
 *     The main points illustrated here are
 *
 *        1) How to declare sufficient memory to find NEV
 *           eigenvalues of largest magnitude.  Other options
 *           are available.
 *
 *        2) Illustration of the reverse communication interface
 *           needed to utilize the top level ARPACK routine ZNAUPD
 *           that computes the quantities needed to construct
 *           the desired eigenvalues and eigenvectors(if requested).
 *
 *        3) How to extract the desired eigenvalues and eigenvectors
 *           using the ARPACK routine ZNEUPD .
 *
 *     The only thing that must be supplied in order to use this
 *     routine on your problem is to change the array dimensions
 *     appropriately, to specify WHICH eigenvalues you want to compute
 *     and to supply a matrix-vector product
 *
 *                         w <-  Av
 *
 *     in place of the call to AV( )  below.
 *
 *     Once usage of this routine is understood, you may wish to explore
 *     the other available options to improve convergence, to solve generalized
 *     problems, etc.  Look at the file ex-complex.doc in DOCUMENTS directory.
 *     This codes implements
 *
 * \Example-1
 *     ... Suppose we want to solve A*x = lambda*x in regular mode,
 *     ... OP = A  and  B = I.
 *     ... Assume "call av (nx,x,y)" computes y = A*x
 *     ... Use mode 1 of ZNAUPD .
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Routines called
 *     znaupd   ARPACK reverse communication interface routine.
 *     zneupd   ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     dznrm2   Level 1 BLAS that computes the norm of a complex vector.
 *     zaxpy    Level 1 BLAS that computes y <- alpha*x+y.
 *     av      Matrix vector multiplication routine that computes A*x.
 *     tv      Matrix vector multiplication routine that computes T*x,
 *             where T is a tridiagonal matrix.  It is used in routine
 *             av.
 *
 * \EndLib
 */
int main()
{
    /* System generated locals */
    int i__1, i__2;
    zomplex z__1;

    /* Local variables */
    double rd[90] /* (3 * MAXNCV) */;
    double rwork[30];

    zomplex d[30];
    zomplex workev[60];

    zomplex sigma;

    int j;
    int ierr, nconv;
    int ishfts, maxitr, mode;

    int ipntr[14];
    int iparam[11];
    bool select[30];
    bool rvec;

    /* ---------------------------------------------------- */
    /* Storage Declarations:                                */
    /*                                                      */
    /* The maximum dimensions for all arrays are            */
    /* set here to accommodate a problem size of            */
    /* N .le. MAXN                                          */
    /*                                                      */
    /* NEV is the number of eigenvalues requested.          */
    /*     See specifications for ARPACK usage below.       */
    /*                                                      */
    /* NCV is the largest number of basis vectors that will */
    /*     be used in the Implicitly Restarted Arnoldi      */
    /*     Process.  Work per major iteration is            */
    /*     proportional to N*NCV*NCV.                       */
    /*                                                      */
    /* You must set:                                        */
    /*                                                      */
    /* MAXN:   Maximum dimension of the A allowed.          */
    /* MAXNEV: Maximum NEV allowed.                         */
    /* MAXNCV: Maximum NCV allowed.                         */
    /* ---------------------------------------------------- */

    /* ----------------------------------------------- */
    /* The following include statement and assignments */
    /* initiate trace output from the internal         */
    /* actions of ARPACK.  See debug.doc in the        */
    /* DOCUMENTS directory for usage.  Initially, the  */
    /* most useful information will be a breakdown of  */
    /* time spent in the various stages of computation */
    /* given by setting mcaupd = 1                     */
    /* ----------------------------------------------- */

    /* ------------------------------- */
    /* See debug.doc for documentation */
    /* ------------------------------- */
    debug_1.ndigit = -3;
    debug_1.logfil = 6;
    debug_1.mcaitr = 0;
    debug_1.mcapps = 0;
    debug_1.mcaupd = 1;
    debug_1.mcaup2 = 0;
    debug_1.mceigh = 0;
    debug_1.mceupd = 0;

    /* ----------------------------------------------- */
    /* The following sets dimensions for this problem. */
    /* ----------------------------------------------- */

    int nx = 10;
    int n = nx * nx;

    /* --------------------------------------------- */
    /*                                               */
    /* Specifications for ARPACK usage are set       */
    /* below:                                        */
    /*                                               */
    /*    1) NEV = 4  asks for 4 eigenvalues to be   */
    /*       computed.                               */
    /*                                               */
    /*    2) NCV = 20 sets the length of the Arnoldi */
    /*       factorization                           */
    /*                                               */
    /*    3) This is a standard problem              */
    /*         (indicated by bmat  = 'I')            */
    /*                                               */
    /*    4) Ask for the NEV eigenvalues of          */
    /*       largest magnitude                       */
    /*         (indicated by which = 'LM')           */
    /*       See documentation in ZNAUPD  for the     */
    /*       other options SM, LR, SR, LI, SI.       */
    /*                                               */
    /* Note: NEV and NCV must satisfy the following  */
    /* conditions:                                   */
    /*              NEV <= MAXNEV                    */
    /*          NEV + 2 <= NCV <= MAXNCV             */
    /*                                               */
    /* --------------------------------------------- */

    int nev = 4;
    int ncv = 20;
    char* bmat = "I";
    char* which = "LM";

    if (n > 256)
    {
        printf(" ERROR with _NSIMP: N is greater than MAXN \n");
        return 0;
    }
    else if (nev > 12)
    {
        printf(" ERROR with _NSIMP: NEV is greater than MAXNEV \n");
        return 0;
    }
    else if (ncv > 30)
    {
        printf(" ERROR with _NSIMP: NCV is greater than MAXNCV \n");
        return 0;
    }

    /* --------------------------------------------------- */
    /*                                                     */
    /* Specification of stopping rules and initial         */
    /* conditions before calling ZNAUPD                     */
    /*                                                     */
    /* TOL  determines the stopping criterion.             */
    /*                                                     */
    /*      Expect                                         */
    /*           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) */
    /*               computed   true                       */
    /*                                                     */
    /*      If TOL .le. 0,  then TOL <- macheps            */
    /*           (machine precision) is used.              */
    /*                                                     */
    /* IDO  is the REVERSE COMMUNICATION parameter         */
    /*      used to specify actions to be taken on return  */
    /*      from ZNAUPD . (see usage below)                 */
    /*                                                     */
    /*      It MUST initially be set to 0 before the first */
    /*      call to ZNAUPD .                                */
    /*                                                     */
    /* INFO on entry specifies starting vector information */
    /*      and on return indicates error codes            */
    /*                                                     */
    /*      Initially, setting INFO=0 indicates that a     */
    /*      random starting vector is requested to         */
    /*      start the ARNOLDI iteration.  Setting INFO to  */
    /*      a nonzero value on the initial call is used    */
    /*      if you want to specify your own starting       */
    /*      vector (This vector must be placed in RESID).  */
    /*                                                     */
    /* The work array WORKL is used in ZNAUPD  as           */
    /* workspace.  Its dimension LWORKL is set as          */
    /* illustrated below.                                  */
    /*                                                     */
    /* --------------------------------------------------- */

    int lworkl = ncv * ncv * 3 + ncv * 5;
    double tol = 0.0;
    int ido = 0;
    int info = 0;

    zomplex* ax = (zomplex*)malloc(n * sizeof(zomplex));
    zomplex* resid = (zomplex*)malloc(n * sizeof(zomplex));
    zomplex* v = (zomplex*)malloc(n * ncv * sizeof(zomplex));
    zomplex* workl = (zomplex*)malloc(lworkl * sizeof(zomplex));
    zomplex* workd = (zomplex*)malloc(3 * n * sizeof(zomplex));

    /* ------------------------------------------------- */
    /* Specification of Algorithm Mode:                  */
    /*                                                   */
    /* This program uses the exact shift strategy        */
    /* (indicated by setting IPARAM(1) = 1).             */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 1 of ZNAUPD  is used     */
    /* (IPARAM(7) = 1). All these options can be changed */
    /* by the user. For details see the documentation in */
    /* ZNAUPD .                                           */
    /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode = 1;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

    /* ---------------------------------------------- */
    /* M A I N   L O O P (Reverse Communication Loop) */
    /* ---------------------------------------------- */

L10:

    /* ------------------------------------------- */
    /* Repeatedly call the routine ZNAUPD  and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */
    znaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, rwork, &info);

    if (ido == -1 || ido == 1)
    {
        /* ----------------------------------------- */
        /* Perform matrix vector multiplication      */
        /*                                           */
        /*                y <--- A*x                 */
        /*                                           */
        /* The user should supply his/her own        */
        /* matrix vector multiplication routine here */
        /* that takes workd(ipntr(1)) as the input   */
        /* vector x , and returns the resulting      */
        /* matrix-vector product y = A*x in the      */
        /* array workd(ipntr(2)).                    */
        /* ----------------------------------------- */

        znsimp_av_(nx, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call ZNAUPD  again. */
        /* --------------------------------------- */

        goto L10;
    }

    /* -------------------------------------- */
    /* Either we have convergence or there is */
    /* an error.                              */
    /* -------------------------------------- */

    if (info < 0)
    {
        /* ------------------------ */
        /* Error message, check the */
        /* documentation in ZNAUPD   */
        /* ------------------------ */

        printf(" \n");
        printf(" Error with _naupd info = %d\n", info);
        printf(" Check the documentation of _naupd\n");
        printf(" \n");

        ierr = info;
        goto EXIT;
    }

    /* ----------------------------------------- */
    /* No fatal errors occurred.                 */
    /* Post-Process using ZNEUPD .                */
    /*                                           */
    /* Computed eigenvalues may be extracted.    */
    /*                                           */
    /* Eigenvectors may be also computed now if  */
    /* desired.  (indicated by rvec = .true.)    */
    /*                                           */
    /* The routine ZNEUPD  now called to do this  */
    /* post processing (Other modes may require  */
    /* more complicated post processing than     */
    /* mode.0)                                   */
    /*                                           */
    /* ----------------------------------------- */

    rvec = true;

    zneupd_(&rvec, "A", select, d, v, &n, &sigma, workev, bmat, &n,which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);

    /* --------------------------------------------- */
    /* Eigenvalues are returned in the one           */
    /* dimensional array D and the corresponding     */
    /* eigenvectors are returned in the first        */
    /* NCONV (=IPARAM(5)) columns of the two         */
    /* dimensional array V if requested.  Otherwise, */
    /* an orthogonal basis for the invariant         */
    /* subspace corresponding to the eigenvalues in  */
    /* D is returned in V.                           */
    /* --------------------------------------------- */

    if (ierr != 0)
    {
        /* ---------------------------------- */
        /* Error condition:                   */
        /* Check the documentation of ZNEUPD . */
        /* ---------------------------------- */

        printf(" \n");
        printf(" Error with _neupd info = %d\n", ierr);
        printf(" Check the documentation of _neupd. \n");
        printf(" \n");

        goto EXIT;
    }

    nconv = iparam[4];
    for (j = 1; j <= nconv; ++j)
    {
        int k = (j - 1) * n;

        /* ------------------------- */
        /* Compute the residual norm */
        /*                           */
        /*   ||  A*x - lambda*x ||   */
        /*                           */
        /* for the NCONV accurately  */
        /* computed eigenvalues and  */
        /* eigenvectors.  (iparam(5) */
        /* indicates how many are    */
        /* accurate to the requested */
        /* tolerance)                */
        /* ------------------------- */

        znsimp_av_(nx, &v[k], ax);
        i__2 = j - 1;
        z__1.r = -d[i__2].r, z__1.i = -d[i__2].i;
        zaxpy_(&n, &z__1, &v[k], &c__1, ax, &c__1);
        i__2 = j - 1;
        rd[j - 1] = d[i__2].r;
        rd[j + 29] = d[i__2].i;
        rd[j + 59] = dznrm2_(&n, ax, &c__1);
        rd[j + 59] /= dlapy2_(&rd[j - 1], &rd[j + 29]);

    }

    /* --------------------------- */
    /* Display computed residuals. */
    /* --------------------------- */

    dmout_(&nconv, &c__3, rd, &c__30, &c_n6, "Ritz values (Real, Imag) and relative residuals");

    /* ----------------------------------------- */
    /* Print additional convergence information. */
    /* ----------------------------------------- */

    if (info == 1)
    {
        printf(" \n");
        printf(" Maximum number of iterations reached.\n");
        printf(" \n");
    }
    else if (info == 3)
    {
        printf(" \n");
        printf(" No shifts could be applied during implicit\n");
        printf(" Arnoldi update try increasing NCV.\n");
        printf(" \n");
    }

    printf(" \n");
    printf("_NSIMP \n");
    printf("====== \n");
    printf(" \n");
    printf(" Size of the matrix is %d\n", n);
    printf(" The number of Ritz values requested is %d\n", nev);
    printf(" The number of Arnoldi vectors generated (NCV) is %d\n", ncv);
    printf(" What portion of the spectrum: %s\n", which);
    printf(" The number of converged Ritz values is %d\n", nconv);
    printf(" The number of Implicit Arnoldi update iterations taken is %d\n", iparam[2]);
    printf(" The number of OP*x is %d\n", iparam[8]);
    printf(" The convergence criterion is %e\n", tol);
    printf(" \n");

EXIT:

    free(ax);
    free(resid);
    free(v);
    free(workl);
    free(workd);

    /* ------------------------- */
    /* Done with program znsimp . */
    /* ------------------------- */

    return ierr;
}

/** Matrix vector subroutine.
 *
 * The matrix used is the convection-diffusion operator
 * discretized using centered difference.
 *
 * Computes w <--- OP*v, where OP is the nx*nx by nx*nx block
 * tridiagonal matrix
 *
 *                  | T -I          |
 *                  |-I  T -I       |
 *             OP = |   -I  T       |
 *                  |        ...  -I|
 *                  |           -I T|
 *
 * derived from the standard central difference discretization
 * of the 2-dimensional convection-diffusion operator
 *              (Laplacian u) + rho*(du/dx)
 * on the unit square with zero boundary condition.
 *
 * The subroutine TV is called to computed y<---T*x.
 */
int znsimp_av_(const int nx, zomplex *v, zomplex *w)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    int j;
    int lo;
    double h2;
    zomplex z;

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    h2 = 1.0 / (double)((nx + 1) * (nx + 1));

    z.r = -1.0 / h2;
    z.i = -0.0;

    znsimp_tv_(nx, &v[1], &w[1]);
    zaxpy_(&nx, &z, &v[nx + 1], &c__1, &w[1], &c__1);

    i__1 = nx - 1;
    for (j = 2; j <= i__1; ++j)
    {
        lo = (j - 1) * nx;
        znsimp_tv_(nx, &v[lo + 1], &w[lo + 1]);
        zaxpy_(&nx, &z, &v[lo - nx + 1], &c__1, &w[lo + 1], &c__1);
        zaxpy_(&nx, &z, &v[lo + nx + 1], &c__1, &w[lo + 1], &c__1);
    }

    lo = (nx - 1) * nx;
    znsimp_tv_(nx, &v[lo + 1], &w[lo + 1]);
    zaxpy_(&nx, &z, &v[lo - nx + 1], &c__1, &w[lo + 1], &c__1);

    return 0;
} /* av_ */

/**
 * Compute the matrix vector multiplication y<---T*x
 * where T is a nx by nx tridiagonal matrix with DD on the
 * diagonal, DL on the subdiagonal, and DU on the superdiagonal
 */
int znsimp_tv_(const int nx, zomplex *x, zomplex *y)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    int j;
    double h, h2, dd, dl, du;

    /* Function Body */
    h = 1.0 / (double)(nx + 1);
    h2 = h * h;
    dd = 4.0 / h2;
    dl = -1.0 / h2 - 50.0 / h;
    du = -1.0 / h2 + 50.0 / h;

    y[0].r = dd * x[0].r + du * x[1].r;
    y[0].i = dd * x[0].i + du * x[1].i;

    i__1 = nx - 1;
    for (j = 1; j < i__1; ++j)
    {
        i__2 = j - 1;
        i__3 = j + 1;
        y[j].r = (dl * x[i__2].r + dd * x[j].r) + du * x[i__3].r;
        y[j].i = (dl * x[i__2].i + dd * x[j].i) + du * x[i__3].i;
    }
    i__1 = nx - 1;
    i__2 = nx - 2;
    y[i__1].r = dl * x[i__2].r + dd * x[i__1].r;
    y[i__1].i = dl * x[i__2].i + dd * x[i__1].i;
    return 0;
} /* tv_ */
