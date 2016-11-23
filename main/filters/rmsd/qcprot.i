/* qcprot.i */
 %module qcprot
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
extern double CalcRMSDRotationalMatrix(double **coords1, double **coords2, const int len, double *rot_matrix);
extern int FastCalcRMSDAndRotation(double *rot, double *A, double *rmsd, double E0, int len, double minScore);
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
extern double CalcRMSDRotationalMatrix(double **coords1, double **coords2, const int len, double *rot_matrix);
extern int FastCalcRMSDAndRotation(double *rot, double *A, double *rmsd, double E0, int len, double minScore);
