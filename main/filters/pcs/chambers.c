#include <stdio.h>
#include <math.h>
#include "mydef.h"
#include "util.h"
#include "chambers.h"

#define ONE_EPS   1.0e-17   /* 1.0 + ONE_EPS > 1.0 */


MINT Givens_Add(MFLOAT **R, MINT q, MFLOAT *W, MFLOAT *y, MFLOAT *x)
{
    MINT     i, j;
    MFLOAT   d, cosx, sinx, c1, c2;
    MFLOAT   Rii, xi;

    for (i = 0; i < q; i++) {
        if (ABS(x[i]) > ONE_EPS * R[i][i]) {
            Rii = R[i][i];
            xi = x[i];
            d = sqrt(Rii*Rii + xi*xi);
            cosx = Rii / d;
            sinx = xi / d;
            R[i][i] = d;
            for (j = i+1; j < q; j++) {
                c1 =  R[i][j]*cosx + x[j]*sinx;
                c2 = -R[i][j]*sinx + x[j]*cosx;
                R[i][j] = c1;
                x[j] = c2;
            }
            c1 =  W[i] * cosx + *y * sinx;
            c2 = -W[i] * sinx + *y * cosx;
            W[i] = c1;
            *y = c2;
        }
    }
    return TRUE;
}
#undef ONE_EPS


MINT Back_Subst(MFLOAT **R, MINT q, MFLOAT *W)
{
    MINT    i, j;

    for (i = q-1; i > 0; i--)  {

        /* Scale the row so the first element is 1 */
        /* This row should have been reduced to all zeros except two */
        if (R[i][i] == 0.0) {
            fprintf(stderr, "Bull shi., singular matrix !!!!!");
            return FALSE;
        }
        W[i] /= R[i][i];
        R[i][i] = (double) 1.0;

        for (j = 0; j < i; j++) {
            /* subtract the row from every one above */
            W[j] += -R[j][i] * W[i];
            R[j][i] = (double) 0.0;
        }
    }

    /* Scale the top row */
    W[0] /= R[0][0];
    R[0][0] = (double) 1.0;

    return TRUE;
}


#ifdef TEST_ME

FILE  *_outF = stdout;   /* stdout file */
FILE  *_errF = stderr;   /* stderr file */

MINT main ()
{
    MINT    i, j, n, q;
    MFLOAT  y, x[100], W[100], **R;

    fscanf(stdin, "%i%i", &n, &q);
    R = DMatrix(q, q, "matrix");
    for (i = 0; i < q; i++) {
        for (j = 0; j < q; j++)
            R[i][j] = 0.0;
        W[i] = 0.0;
    }

    for (i = 0; i < n; i++) {
        fscanf(stdin, "%lf", &y);
        for (j = 0; j < q; j++)
            fscanf(stdin, "%lf", &x[j]);
        Givens_Add(R, q, W, &y, x);
    }
    for (i = 0; i < q; i++) {
        for (j = 0; j < q; j++) {
            if (j < i)
                R[i][j] = R[j][i];
            fprintf(stdout, "%10.5lf", R[i][j]);
        }
        fprintf(stdout, "%10.5lf\n", W[i]);
    }

    Back_Subst(R, q, W);
    fprintf(stdout, "solution:\n");
    for (i = 0; i < q; i++)
        fprintf(stdout, "%10.5lf\n", W[i]);


    GR_FREE(R);
    GR_FREE(*R);
    return TRUE;
}
#endif
