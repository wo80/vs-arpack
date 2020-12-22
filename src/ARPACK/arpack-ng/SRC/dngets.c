/* arpack-ng\SRC\dngets.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \BeginDoc
 *
 * \Name: dngets
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
 *  call dngets
 *     ( ISHIFT, WHICH, KEV, NP, RITZR, RITZI, BOUNDS, SHIFTR, SHIFTI )
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
 *          'LR' -> want the KEV eigenvalues of largest real part.
 *          'SR' -> want the KEV eigenvalues of smallest real part.
 *          'LI' -> want the KEV eigenvalues of largest imaginary part.
 *          'SI' -> want the KEV eigenvalues of smallest imaginary part.
 *
 *  KEV      Integer.  (INPUT/OUTPUT)
 *           INPUT: KEV+NP is the size of the matrix H.
 *           OUTPUT: Possibly increases KEV by one to keep complex conjugate
 *           pairs together.
 *
 *  NP       Integer.  (INPUT/OUTPUT)
 *           Number of implicit shifts to be computed.
 *           OUTPUT: Possibly decreases NP by one to keep complex conjugate
 *           pairs together.
 *
 *  RITZR,  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
 *  RITZI   On INPUT, RITZR and RITZI contain the real and imaginary
 *          parts of the eigenvalues of H.
 *          On OUTPUT, RITZR and RITZI are sorted so that the unwanted
 *          eigenvalues are in the first NP locations and the wanted
 *          portion is in the last KEV locations.  When exact shifts are
 *          selected, the unwanted part corresponds to the shifts to
 *          be applied. Also, if ISHIFT .eq. 1, the unwanted eigenvalues
 *          are further sorted so that the ones with largest Ritz values
 *          are first.
 *
 *  BOUNDS  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
 *          Error bounds corresponding to the ordering in RITZ.
 *
 *  SHIFTR, SHIFTI  *** USE deprecated as of version 2.1. ***
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Local variables:
 *     xxxxxx  real
 *
 * \Routines called:
 *     dsortc  ARPACK sorting routine.
 *     dcopy   Level 1 BLAS that copies one vector to another .
 *
 * \Author
 *     Danny Sorensen               Phuong Vu
 *     Richard Lehoucq              CRPC / Rice University
 *     Dept. of Computational &     Houston, Texas
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * \Revision history:
 *     xx/xx/92: Version ' 2.1'
 *
 * \SCCS Information: @(#)
 * FILE: ngets.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
 *
 * \Remarks
 *     1. xxxx
 *
 * \EndLib
 */
int dngets_(int *ishift, char *which, int *kev, int *np, double *ritzr,
            double *ritzi, double *bounds,double *shiftr, double *shifti)
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

    msglvl = debug_1.mngets;

    /* -------------------------------------------------- */
    /* LM, SM, LR, SR, LI, SI case.                       */
    /* Sort the eigenvalues of H into the desired order   */
    /* and apply the resulting order to BOUNDS.           */
    /* The eigenvalues are sorted so that the wanted part */
    /* are always in the last KEV locations.              */
    /* We first do a pre-processing sort in order to keep */
    /* complex conjugate pairs together                   */
    /* -------------------------------------------------- */

    i__1 = *kev + *np;

    if (strcmp(which, "LM") == 0)
    {
        dsortc_("LR", &c_true, &i__1, ritzr, ritzi, bounds);
    }
    else if (strcmp(which, "SM") == 0)
    {
        dsortc_("SR", &c_true, &i__1, ritzr, ritzi, bounds);
    }
    else if (strcmp(which, "LR") == 0)
    {
        dsortc_("LM", &c_true, &i__1, ritzr, ritzi, bounds);
    }
    else if (strcmp(which, "SR") == 0)
    {
        dsortc_("SM", &c_true, &i__1, ritzr, ritzi, bounds);
    }
    else if (strcmp(which, "LI") == 0)
    {
        dsortc_("LM", &c_true, &i__1, ritzr, ritzi, bounds);
    }
    else if (strcmp(which, "SI") == 0)
    {
        dsortc_("SM", &c_true, &i__1, ritzr, ritzi, bounds);
    }

    dsortc_(which, &c_true, &i__1, ritzr, ritzi, bounds);

    /* ----------------------------------------------------- */
    /* Increase KEV by one if the ( ritzr(np),ritzi(np) )    */
    /* = ( ritzr(np+1),-ritzi(np+1) ) and ritz(np) .ne. zero */
    /* Accordingly decrease NP by one. In other words keep   */
    /* complex conjugate pairs together.                     */
    /* ----------------------------------------------------- */

    if (ritzr[*np] - ritzr[*np - 1] == 0.0 && ritzi[*np] + ritzi[*np - 1] == 0.0)
    {
        --(*np);
        ++(*kev);
    }

    if (*ishift == 1)
    {
        /* ----------------------------------------------------- */
        /* Sort the unwanted Ritz values used as shifts so that  */
        /* the ones with largest Ritz estimates are first        */
        /* This will tend to minimize the effects of the         */
        /* forward instability of the iteration when they shifts */
        /* are applied in subroutine dnapps.                     */
        /* Be careful and use 'SR' since we want to sort BOUNDS! */
        /* ----------------------------------------------------- */

        dsortc_("SR", &c_true, np, bounds, ritzr, ritzi);
    }

#ifndef NO_TIMER
    arscnd_(&t1);
    timing_1.tngets += t1 - t0;
#endif

#ifndef NO_TRACE
    if (msglvl > 0)
    {
        i__1 = *kev + *np;
        ivout_(&c__1, kev, &debug_1.ndigit, "_ngets: KEV is");
        ivout_(&c__1, np, &debug_1.ndigit, "_ngets: NP is");
        dvout_(&i__1, ritzr, &debug_1.ndigit, "_ngets: Eigenvalues of current H matrix -- float part");
        dvout_(&i__1, ritzi, &debug_1.ndigit, "_ngets: Eigenvalues of current H matrix -- imag part");
        dvout_(&i__1, bounds, &debug_1.ndigit, "_ngets: Ritz estimates of the current KEV+NP Ritz values");
    }
#endif

    return 0;

    /* ------------- */
    /* End of dngets */
    /* ------------- */

} /* dngets_ */

