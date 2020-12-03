/* D:\Projekte\ARPACK\arpack-ng\SRC\ssgets.f -- translated by f2c (version 20100827). */

#include "arpack.h"

/**
 * \BeginDoc
 *
 * \Name: ssgets
 *
 * \Description:
 *  Given the eigenvalues of the symmetric tridiagonal matrix H,
 *  computes the NP shifts AMU that are zeros of the polynomial of
 *  degree NP which filters out components of the unwanted eigenvectors
 *  corresponding to the AMU's based on some given criteria.
 *
 *  NOTE: This is called even in the case of user specified shifts in
 *  order to sort the eigenvalues, and error bounds of H for later use.
 *
 * \Usage:
 *  call ssgets
 *     ( ISHIFT, WHICH, KEV, NP, RITZ, BOUNDS, SHIFTS )
 *
 * \Arguments
 *  ISHIFT  Integer.  (INPUT)
 *          Method for selecting the implicit shifts at each iteration.
 *          ISHIFT = 0: user specified shifts
 *          ISHIFT = 1: exact shift with respect to the matrix H.
 *
 *  WHICH   Character*2.  (INPUT)
 *          Shift selection criteria.
 *          'LM' -> KEV eigenvalues of largest magnitude are retained.
 *          'SM' -> KEV eigenvalues of smallest magnitude are retained.
 *          'LA' -> KEV eigenvalues of largest value are retained.
 *          'SA' -> KEV eigenvalues of smallest value are retained.
 *          'BE' -> KEV eigenvalues, half from each end of the spectrum.
 *                  If KEV is odd, compute one more from the high end.
 *
 *  KEV      Integer.  (INPUT)
 *          KEV+NP is the size of the matrix H.
 *
 *  NP      Integer.  (INPUT)
 *          Number of implicit shifts to be computed.
 *
 *  RITZ    Real array of length KEV+NP.  (INPUT/OUTPUT)
 *          On INPUT, RITZ contains the eigenvalues of H.
 *          On OUTPUT, RITZ are sorted so that the unwanted eigenvalues
 *          are in the first NP locations and the wanted part is in
 *          the last KEV locations.  When exact shifts are selected, the
 *          unwanted part corresponds to the shifts to be applied.
 *
 *  BOUNDS  Real array of length KEV+NP.  (INPUT/OUTPUT)
 *          Error bounds corresponding to the ordering in RITZ.
 *
 *  SHIFTS  Real array of length NP.  (INPUT/OUTPUT)
 *          On INPUT:  contains the user specified shifts if ISHIFT = 0.
 *          On OUTPUT: contains the shifts sorted into decreasing order
 *          of magnitude with respect to the Ritz estimates contained in
 *          BOUNDS. If ISHIFT = 0, SHIFTS is not modified on exit.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * \Local variables:
 *     xxxxxx  float
 *
 * \Routines called:
 *     ssortr  ARPACK utility sorting routine.
 *     ivout   ARPACK utility routine that prints int32_ts.
 *     arscnd  ARPACK utility routine for timing.
 *     svout   ARPACK utility routine that prints vectors.
 *     scopy   Level 1 BLAS that copies one vector to another.
 *     sswap   Level 1 BLAS that swaps the contents of two vectors.
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
 *     xx/xx/93: Version ' 2.1'
 *
 * \SCCS Information: @(#)
 * FILE: sgets.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2
 *
 * \Remarks
 *
 * \EndLib
 */

int ssgets_(int32_t *ishift, char *which, int32_t *kev, 
	int32_t *np, float *ritz, float *bounds, float *shifts)
{
    /* System generated locals */
    int32_t i__1;

    /* Builtin functions */
    int32_t s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static float t0, t1;
    int32_t kevd2;
    int32_t msglvl;

     /* --------------------- */
     /* Executable Statements */
     /* --------------------- */

     /* ----------------------------- */
     /* Initialize timing statistics  */
     /* & message level for debugging */
     /* ----------------------------- */

    /* Parameter adjustments */
    --shifts;
    --bounds;
    --ritz;

    /* Function Body */
    arscnd_(&t0);
    msglvl = debug_1.msgets;

    if (s_cmp(which, "BE", (ftnlen)2, (ftnlen)2) == 0) {

        /* --------------------------------------------------- */
        /* Both ends of the spectrum are requested.            */
        /* Sort the eigenvalues into algebraically increasing  */
        /* order first then swap high end of the spectrum next */
        /* to low end in appropriate locations.                */
        /* NOTE: when np < floor(kev/2) be careful not to swap */
        /* overlapping locations.                              */
        /* --------------------------------------------------- */

	i__1 = *kev + *np;
	ssortr_("LA", &c_true, &i__1, &ritz[1], &bounds[1]);
	kevd2 = *kev / 2;
	if (*kev > 1) {
	    i__1 = min(kevd2,*np);
	    sswap_(&i__1, &ritz[1], &c__1, &ritz[max(kevd2,*np) + 1], &c__1);
	    i__1 = min(kevd2,*np);
	    sswap_(&i__1, &bounds[1], &c__1, &bounds[max(kevd2,*np) + 1], &
		    c__1);
	}

    } else {

        /* -------------------------------------------------- */
        /* LM, SM, LA, SA case.                               */
        /* Sort the eigenvalues of H into the desired order   */
        /* and apply the resulting order to BOUNDS.           */
        /* The eigenvalues are sorted so that the wanted part */
        /* are always in the last KEV locations.              */
        /* -------------------------------------------------- */

	i__1 = *kev + *np;
	ssortr_(which, &c_true, &i__1, &ritz[1], &bounds[1]);
    }

    if (*ishift == 1 && *np > 0) {

        /* ----------------------------------------------------- */
        /* Sort the unwanted Ritz values used as shifts so that  */
        /* the ones with largest Ritz estimates are first.       */
        /* This will tend to minimize the effects of the         */
        /* forward instability of the iteration when the shifts  */
        /* are applied in subroutine ssapps.                     */
        /* ----------------------------------------------------- */

	ssortr_("SM", &c_true, np, &bounds[1], &ritz[1]);
	scopy_(np, &ritz[1], &c__1, &shifts[1], &c__1);
    }

    arscnd_(&t1);
    timing_1.tsgets += t1 - t0;

    if (msglvl > 0) {
	ivout_(&debug_1.logfil, &c__1, kev, &debug_1.ndigit, "_sgets: KEV is",
		 (ftnlen)14);
	ivout_(&debug_1.logfil, &c__1, np, &debug_1.ndigit, "_sgets: NP is", (
		ftnlen)13);
	i__1 = *kev + *np;
	svout_(&debug_1.logfil, &i__1, &ritz[1], &debug_1.ndigit, "_sgets: E"
		"igenvalues of current H matrix", (ftnlen)39);
	i__1 = *kev + *np;
	svout_(&debug_1.logfil, &i__1, &bounds[1], &debug_1.ndigit, "_sgets:"
		" Associated Ritz estimates", (ftnlen)33);
    }

    return 0;

     /* ------------- */
     /* End of ssgets */
     /* ------------- */

} /* ssgets_ */

