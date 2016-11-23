#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mydef.h"
#include "util.h"
#include "chambers.h"

MINT eigen(MFLOAT *a, MFLOAT *r__, MINT *n, MINT *mv);


#define FACT_USI_PRECALC_FOR_A_3 795.77471545947673803312
void fill_A_lines(MINT n, MFLOAT **A_line, MFLOAT xM, MFLOAT yM, MFLOAT zM, MFLOAT **xyz)
{
    MINT   i;
    MFLOAT value_1_4_PI_r5;
    MFLOAT r2;
    MFLOAT d;
    MFLOAT r5;
    MFLOAT x2, y2, z2;
    MFLOAT v_x, v_y, v_z;

    for (i = 0; i < n; i++) {
        v_x = xyz[i][0] - xM;
        v_y = xyz[i][1] - yM;
        v_z = xyz[i][2] - zM;
        x2 = v_x * v_x;
        y2 = v_y * v_y;
        z2 = v_z * v_z;

        r2 = x2 + y2 + z2;
        d = sqrt(r2);
        r5 = r2 * r2 * d;

        value_1_4_PI_r5 = FACT_USI_PRECALC_FOR_A_3 / r5;

        A_line[i][0] = value_1_4_PI_r5 * (x2 - z2);
        A_line[i][1] = value_1_4_PI_r5 * (2.0 * v_x * v_y);
        A_line[i][2] = value_1_4_PI_r5 * (2.0 * v_x * v_z);
        A_line[i][3] = value_1_4_PI_r5 * (y2 - z2);
        A_line[i][4] = value_1_4_PI_r5 * (2.0 * v_y * v_z);
    }
}

double **MakeDMatrix(int i, int j)
{
    return DMatrix(i, j, "DArray");
}
void FreeDMatrix(double **a)
{
    if ((*a)) free(*a);
    if ((a)) free(a);
}

void FreeDArray(double *a)
{
    if ((a)) free(a);
}
double GetDArray(int i, int j, double **a)
{
    return a[i][j];
}
void SetDArray(int i, int j, double **a, double v)
{
    a[i][j] = v;
}
double GetDvector(int i, double *a)
{
    return a[i];
}
void SetDvector(int i, double *a, double v)
{
    a[i] = v;
}
double *MakeDvector(int n)
{
    return Dvector(n, "r vector");
}
int *MakeIvector(int n)
{
    return Ivector(n, "r vector");
}
int GetIvector(int i, int *a)
{
    return a[i];
}
void SetIvector(int i, int *a, int v)
{
    a[i] = v;
}



void Saupe2X(double *s, double *X)
{
    int    n = 3, mode = 1;
    double m[6], dummy, ev[3];
    m[0] = s[0];
    m[1] = s[1];
    m[2] = s[2];
    m[3] = s[3];
    m[4] = s[4];
    m[5] = -s[0]-s[3];
    eigen(m, &dummy, &n, &mode);
    ev[0] = ABS(m[0]);
    ev[1] = ABS(m[2]);
    ev[2] = ABS(m[5]);

    /* sort ev  (ugly but fast) */
    if (ev[0] > ev[1]) {
        if (ev[1] < ev[2]) {
            if (ev[0] < ev[2]) {
                ev[0] = m[5];
                ev[1] = m[0];
                ev[2] = m[2];
            }
            else {
                ev[0] = m[0];
                ev[1] = m[5];
                ev[2] = m[2];
            }
        }
        else {
            ev[0] = m[0];
            ev[1] = m[2];
            ev[2] = m[5];
        }
    }
    else {
        if (ev[1] < ev[2]) {
            ev[0] = m[5];
            ev[1] = m[2];
            ev[2] = m[0];
        }
        else {
            if (ev[0] < ev[2]) {
                ev[0] = m[2];
                ev[1] = m[5];
                ev[2] = m[0];
            }
            else {
                ev[0] = m[2];
                ev[1] = m[0];
                ev[2] = m[5];
            }
        }
    }
    X[0] = ev[0] - 0.5*(ev[1] + ev[2]);
    X[1] = ev[2] - ev[1];
}


double rfastT1FM_multi(int nM, double *rMx, double *rMy, double *rMz, int nsets, int npcs, double **xyz, double **pcs, double **tensor, double **Xaxrh_range)
{
    MINT    i, ii, j, jj, q = 5;
    MFLOAT  **W, **R, Xaxrh[2];

    MFLOAT  **A_line, c, cbest, dx;
    MFLOAT  Atmp[5];


    R = DMatrix(q, q, "R matrix");
    W = DMatrix(nsets, q, "W matrix");
    A_line = DMatrix(npcs, q, "A matrix");


    cbest = 1.0e+30;
    for (ii = 0; ii < nM; ii++) {

        c = 0.0;
        for (jj = 0; jj < nsets; jj++) {

            for (i = 0; i < q; i++) {
                for (j = 0; j < q; j++)
                    R[i][j] = 0.0;
                W[jj][i] = 0.0;
            }

            fill_A_lines(npcs, A_line, rMx[ii], rMy[ii], rMz[ii], xyz);
            for (i = 0; i < npcs; i++) {
                if (pcs[jj][i] != 999.999) {
                    Atmp[0] = A_line[i][0];
                    Atmp[1] = A_line[i][1];
                    Atmp[2] = A_line[i][2];
                    Atmp[3] = A_line[i][3];
                    Atmp[4] = A_line[i][4];
                    dx = pcs[jj][i];
                    Givens_Add(R, q, W[jj], &dx, Atmp);
                }
            }

            for (i = 0; i < q; i++)
                for (j = 0; j < i; j++)
                    R[i][j] = R[j][i];

            Back_Subst(R, q, W[jj]);

            for (i = 0; i < npcs; i++) {
                if (pcs[jj][i] != 999.999) {
                    dx = W[jj][0]*A_line[i][0]+W[jj][1]*A_line[i][1]+W[jj][2]*A_line[i][2]+W[jj][3]*A_line[i][3]+W[jj][4]*A_line[i][4] - pcs[jj][i];
                    c += dx*dx;
                }
            }

            Saupe2X(W[jj], Xaxrh);

            if ((Xaxrh[0] < Xaxrh_range[jj][0]) || (Xaxrh[0] > Xaxrh_range[jj][1]) || (Xaxrh[1] < Xaxrh_range[jj][2]) || (Xaxrh[1] > Xaxrh_range[jj][3])) {
                c = cbest + 1.0;
                jj = nsets;             /* finish nsets loop */
            }
        }


        if (c < cbest) {
            for (jj = 0; jj < nsets; jj++) {
                tensor[jj][0] = rMx[ii];
                tensor[jj][1] = rMy[ii];
                tensor[jj][2] = rMz[ii];
                for (j = 0; j < 5; j++)
                    tensor[jj][j+3] = W[jj][j];
                cbest = c;
            }
        }
    }

    GR_FREE(*R);
    GR_FREE(R);
    GR_FREE(*W);
    GR_FREE(W);
    GR_FREE(*A_line);
    GR_FREE(A_line);
    return cbest;
}
