/* fastT1FM.i */
 %module fastT1FM
 %{
 /* Put header files here or function declarations like below */
     extern double **MakeDMatrix(int i, int j);
     extern void FreeDMatrix(double **a);
     extern void FreeDArray(double *a);
     extern double GetDArray(int i, int j, double **a);
     extern void SetDArray(int i, int j, double **a, double v);
     extern double GetDvector(int i, double *a);
     extern void SetDvector(int i, double *a, double v);
     extern double *MakeDvector(int n);
     extern int *MakeIvector(int n);
     extern int GetIvector(int i, int *a);
     extern void SetIvector(int i, int *a, int v);
     extern void Saupe2X(double *s, double *X);
     extern double rfastT1FM_multi(int nM, double *rMx, double *rMy, double *rMz, int nsets, int npcs, double **xyz, double **pcs, double **tensor, double **Xaxrh_range);
 %}
 
extern double **MakeDMatrix(int i, int j);
extern void FreeDMatrix(double **a);
extern void FreeDArray(double *a);
extern double GetDArray(int i, int j, double **a);
extern void SetDArray(int i, int j, double **a, double v);
extern double GetDvector(int i, double *a);
extern void SetDvector(int i, double *a, double v);
extern double *MakeDvector(int n);
extern int *MakeIvector(int n);
extern int GetIvector(int i, int *a);
extern void SetIvector(int i, int *a, int v);
extern void Saupe2X(double *s, double *X);
extern double rfastT1FM_multi(int nM, double *rMx, double *rMy, double *rMz, int nsets, int npcs, double **xyz, double **pcs, double **tensor, double **Xaxrh_range);
