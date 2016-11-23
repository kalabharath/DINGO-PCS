#include <stdio.h>
#include <math.h>
#include "mydef.h"
#include "util.h"

MINT eigen(MFLOAT *a, MFLOAT *r__, MINT *n, MINT *mv)
{
    MINT    i__1, i__2, i__3;
    MFLOAT  d__1;

    MFLOAT  cosx, sinx, cosx2, sinx2;
    MINT    i__, j, k, l, m;
    MFLOAT  x, y, range, anorm, sincs, anrmx;
    MINT    ia, ij, il, im, ll, lm, iq, mm, jq, lq, mq, ind, ilq, imq, ilr, imr;
    MFLOAT  thr;


/* CCCCC W.F. VAN GUNSTEREN, CAMBRIDGE, JUNE 1979 CCCCCCCCCCCCCCCCCCCCCCCC
                                                                      C 
     SUBROUTINE EIGEN (A,R,N,MV)                                      C 
                                                                      C 
         EIGEN COMPUTES EIGENVALUES AND EIGENVECTORS OF THE REAL      C
     SYMMETRIC N*N MATRIX A, USING THE DIAGONALIZATION METHOD         C 
     DESCRIBED IN "MATHEMATICAL METHODS FOR DIGITAL COMPUTERS", EDS.  C 
     A.RALSTON AND H.S.WILF, WILEY, NEW YORK, 1962, CHAPTER 7.        C 
     IT HAS BEEN COPIED FROM THE IBM SCIENTIFIC SUBROUTINE PACKAGE.   C 
                                                                      C 
     A(1..N*(N+1)/2) = MATRIX TO BE DIAGONALIZED, STORED IN SYMMETRIC C 
                       STORAGE MODE, VIZ. THE I,J-TH ELEMENT (I.GE.J) C 
                       IS STORED AT THE LOCATION K=I*(I-1)/2+J IN A;  C 
                       THE EIGENVALUES ARE DELIVERED IN DESCENDING    C 
                       ORDER ON THE DIAGONAL, VIZ. AT THE LOCATIONS   C 
                       K=I*(I+1)/2                                    C 
     R(1..N,1..N) = DELIVERED WITH THE CORRESPONDING EIGENVECTORS     C 
                    STORED COLUMNWISE                                 C 
     N = ORDER OF MATRICES A AND R                                    C 
     MV = 0 : EIGENVALUES AND EIGENVECTORS ARE COMPUTED               C 
        = 1 : ONLY EIGENVALUES ARE COMPUTED                           C 
                                                                      C 
 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*/

    /* Parameter adjustments */
    --r__;
    --a;



    range = 1e-12;
    if (*mv != 1) {
        iq = -(*n);
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            iq += *n;
            i__2 = *n;
            for (i__ = 1; i__ <= i__2; ++i__) {
                ij = iq + i__;
                r__[ij] = 0.;
                if (i__ - j == 0)
                    r__[ij] = 1.;
            }
        }
    }

/* *****COMPUTE INITIAL AND FINAL NORMS (ANORM AND ANRMX) */
    anorm = 0.0;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = *n;
	for (j = i__; j <= i__1; ++j) {
	    if (i__ - j != 0) {
                ia = i__ + (j * j - j) / 2;
                anorm += a[ia] * a[ia];
            }
	}
    }
    if (anorm > 0.0) {
        anorm = sqrt(anorm) * 1.414;
        anrmx = anorm * range / (MFLOAT) (*n);

/* *****INITIALIZE INDICATORS AND COMPUTE THRESHOLD, THR */
        thr = anorm;

        do {  /* ----- while (thr - anrmx > 0.0) ------ */
            thr /= (MFLOAT) (*n);

            ind = TRUE;
            while (ind) {
                ind = FALSE;

                l = 1;
                for (l = 1; l <= *n-1; l++) {
                    for (m = l+1; m <= *n; m++) {

/* *****COMPUT SIN AND COS */

                        mq = (m * m - m) / 2;
                        lq = (l * l - l) / 2;
                        lm = l + mq;
                        d__1 = a[lm];
                        if (ABS(d__1) - thr >= 0.0) {
                            ind = TRUE;
                            ll = l + lq;
                            mm = m + mq;
                            x = (a[ll] - a[mm]) * .5;
                            y = -a[lm] / sqrt(a[lm] * a[lm] + x * x);
                            if (x < 0.0)
                                y = -y;
                            sinx = y / sqrt((sqrt(1.0 - y * y) + 1.0) * 2.0);
                            sinx2 = sinx * sinx;
                            cosx = sqrt(1.0 - sinx2);
                            cosx2 = cosx * cosx;
                            sincs = sinx * cosx;

/* *****ROTATE L AND M COLUMNS */
                            ilq = *n * (l - 1);
                            imq = *n * (m - 1);
                            i__1 = *n;
                            for (i__ = 1; i__ <= i__1; ++i__) {
                                iq = (i__ * i__ - i__) / 2;
                                if (i__ - l != 0) {

                                    i__2 = i__ - m;
                                    if (i__2 != 0) {
                                        if (i__2 < 0)
                                            im = i__ + mq;
                                        else
                                            im = m + iq;

                                        if (i__ - l >= 0)
                                            il = l + iq;
                                        else
                                            il = i__ + lq;
                                        x = a[il] * cosx - a[im] * sinx;
                                        a[im] = a[il] * sinx + a[im] * cosx;
                                        a[il] = x;
                                    } /* ---- (i__2 != 0) ---- */
                                } /* ------ if (i__ - l != 0) ---- */

                                if (*mv != 1) {
                                    ilr = ilq + i__;
                                    imr = imq + i__;
                                    x = r__[ilr] * cosx - r__[imr] * sinx;
                                    r__[imr] = r__[ilr] * sinx + r__[imr] * cosx;
                                    r__[ilr] = x;
                                }
                            }
                            x = a[lm] * 2. * sincs;
                            y = a[ll] * cosx2 + a[mm] * sinx2 - x;
                            x = a[ll] * sinx2 + a[mm] * cosx2 + x;
                            a[lm] = (a[ll] - a[mm]) * sincs + a[lm] * (cosx2 - sinx2);
                            a[ll] = y;
                            a[mm] = x;

                        } /* --- if ((d__1 = a[lm], ABS(d__1)) - thr >= 0.0) -----
*/
/* *****TESTS FOR COMPLETION */

/* *****TEST FOR M = LAST COLUMN */
                    } /* ---- for (m = l+1; m <= *n; m++)  ---- */

/* *****TEST FOR L = SECOND FROM LAST COLUMN */
                } /* ----- for (l = 1; l < *n-1; l++) ------ */

            } /* --- while (ind) --- */
/* *****COMPARE THRESHOLD WITH FINAL NORM */

        } while (thr > anrmx);


    } /* ---- if (anorm > 0) ------ */

/* *****SORT EIGENVALUES AND EIGENVECTORS */

    iq = -(*n);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iq += *n;
	ll = i__ + (i__ * i__ - i__) / 2;
	jq = *n * (i__ - 2);
	i__2 = *n;
	for (j = i__; j <= i__2; ++j) {
	    jq += *n;
	    mm = j + (j * j - j) / 2;
	    if (a[ll] - a[mm] < 0.0) {
                x = a[ll];
                a[ll] = a[mm];
                a[mm] = x;
                if (*mv != 1) {
                    i__3 = *n;
                    for (k = 1; k <= i__3; ++k) {
                        ilr = iq + k;
                        imr = jq + k;
                        x = r__[ilr];
                        r__[ilr] = r__[imr];
                        r__[imr] = x;
                    }
                }
            } /* ---- if (a[ll] - a[mm] < 0.0) ---- */
	}
    }

    return TRUE;
} /* eigen */

