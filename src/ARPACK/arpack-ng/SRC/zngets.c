/* arpack-ng\SRC\zngets.f -- translated by f2c (version 20100827). */

#include "arpack_internal.h"

/**
 * \BeginDoc
 *
 * \Name: zngets
 *
 * \Description:
 *  Given the eigenvalues of the upper Hessenberg matrix H,
 *  computes the NP shifts AMU that are zeros of the polynomial of
 *  degree NP which filters out components of the unwanted eigenvectors
 *  corresponding to the AMU's based on some given criteria.
 *
 *  NOTE: call this even in the case of user specified shifts in order
 *  to sort the eigenvalues, and error bounds of H for later use.
 *
 * \Usage:
 *  call zngets
 *      ( ISHIFT, WHICH, KEV, NP, RITZ, BOUNDS )
 *
 * \Arguments
 *  ISHIFT  Integer.  (INPUT)
 *          Method for selecting the implicit shifts at each iteration.
 *          ISHIFT = 0: user specified shifts
 *          ISHIFT = 1: exact shift with respect to the matrix H.
 *
 *  WHICH   Character*2.  (INPUT)
 *          Shift selection criteria.
 *          'LM' -> want the KEV eigenvalues of largest magnitude.
 *          'SM' -> want the KEV eigenvalues of smallest magnitude.
 *          'LR' -> want the KEV eigenvalues of largest REAL part.
 *          'SR' -> want the KEV eigenvalues of smallest REAL part.
 *          'LI' -> want the KEV eigenvalues of largest imaginary part.
 *          'SI' -> want the KEV eigenvalues of smallest imaginary part.
 *
 *  KEV     Integer.  (INPUT)
 *          The number of desired eigenvalues.
 *
 *  NP      Integer.  (INPUT)
 *          The number of shifts to compute.
 *
 *  RITZ    Complex*16 array of length KEV+NP.  (INPUT/OUTPUT)
 *          On INPUT, RITZ contains the the eigenvalues of H.
 *          On OUTPUT, RITZ are sorted so that the unwanted
 *          eigenvalues are in the first NP locations and the wanted
 *          portion is in the last KEV locations.  When exact shifts are
 *          selected, the unwanted part corresponds to the shifts to
 *          be applied. Also, if ISHIFT .eq. 1, the unwanted eigenvalues
 *          are further sorted so that the ones with largest Ritz values
 *          are first.
 *
 *  BOUNDS  Complex*16 array of length KEV+NP.  (INPUT/OUTPUT)
 *          Error bounds corresponding to the ordering in RITZ.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Local variables:
 *     xxxxxx  Complex*16
 *
 * \Routines called:
 *     zsortc  ARPACK sorting routine.
 *     ivout   ARPACK utility routine that prints integers.
 *     arscnd  ARPACK utility routine for timing.
 *     zvout   ARPACK utility routine that prints vectors.
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
 * FILE: ngets.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2
 *
 * \Remarks
 *     1. This routine does not keep complex conjugate pairs of
 *        eigenvalues together.
 *
 * \EndLib
 */
int zngets_(int *ishift, char *which, int *kev, int *np, zomplex *ritz,
            zomplex *bounds)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static float t0, t1;
    int msglvl;


    /* ----------------------------- */
    /* Initialize timing statistics  */
    /* & message level for debugging */
    /* ----------------------------- */

    /* Function Body */
#ifndef NO_TIMER
    arscnd_(&t0);
#endif

    msglvl = debug_1.mcgets;

    i__1 = *kev + *np;
    zsortc_(which, &c_true, &i__1, ritz, bounds);

    if (*ishift == 1)
    {
        /* ----------------------------------------------------- */
        /* Sort the unwanted Ritz values used as shifts so that  */
        /* the ones with largest Ritz estimates are first        */
        /* This will tend to minimize the effects of the         */
        /* forward instability of the iteration when the shifts  */
        /* are applied in subroutine znapps.                     */
        /* Be careful and use 'SM' since we want to sort BOUNDS! */
        /* ----------------------------------------------------- */

        zsortc_("SM", &c_true, np, bounds, ritz);
    }

#ifndef NO_TIMER
    arscnd_(&t1);
    timing_1.tcgets += t1 - t0;
#endif

#ifndef NO_TRACE
    if (msglvl > 0)
    {
        ivout_(&c__1, kev, &debug_1.ndigit, "_ngets: KEV is");
        ivout_(&c__1, np, &debug_1.ndigit, "_ngets: NP is");
        zvout_(&i__1, ritz, &debug_1.ndigit, "_ngets: Eigenvalues of current H matrix ");
        zvout_(&i__1, bounds, &debug_1.ndigit, "_ngets: Ritz estimates of the current KEV+NP Ritz values");
    }
#endif

    return 0;

    /* ------------- */
    /* End of zngets */
    /* ------------- */

} /* zngets_ */

