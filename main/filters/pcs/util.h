/* "$Id: util.h,v 1.5 1995/08/07 03:37:34 huber Exp huber $" */

#ifndef UTIL_H 
#define UTIL_H


#include <time.h>
#include <sys/time.h>

#ifndef MALLOC_H
#include <malloc.h>
#endif


#define FALSE 0
#define TRUE  !FALSE

#ifndef PI
#define PI 3.14159265358979323846
#endif

#define EOL   10
#define TAB    9
#define BLANK 32

#define MAX_LINE_LEN  200

#define COMMENT_CHAR    '#'
#define IS_COMMENT(a) ((a) == COMMENT_CHAR)
/* !strncmp('COMMENT_CHAR', (a), 1) */


/* ----- some macros for functions ----------------------- */

/* ----- absolute value ---------------------------------- */
#define	ABS(x) ((x) < 0 ? -(x) : (x))

/* ----- maximum value of a and b ------------------------ */
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ----- minimum value of a and b ------------------------ */
#define MIN(a,b) ((a) > (b) ? (b) : (a))

/* ---- fortran SIGN function, transfers sign of 2. argument to 1. --- */
#define SIGN(x,y) ((y) > 0.0 ? ABS(x) : -ABS(x))

/* ---- bit operations ----------------------------------- */

#define SETBIT(a, b) ((a) + (1 << (b)))

#define IS_BITSET(a, b) (((a) & (1 << (b))) >> b)


/* ---- save macro for free ----------------------------------- */

#define GR_FREE(a) if ((a)) free(a)

/* ----------- proc declarations ------------------------------------------ */

void error_exit(char *message);

void Get_Line(FILE *f, char *line, MINT *nr_line);

void *Gr_malloc(size_t size, char *message);

/* void Gr_free(void *old_mem); */

MINT *Ivector(MINT length, char *message);

MFLOAT *Dvector(MINT length, char *message);

float *Fvector(MINT length, char *message);

RPoint *Rvector(MINT length, char *message);

MINT **IMatrix(MINT n_rows, MINT n_cols, char *message);

MFLOAT **DMatrix(MINT n_rows, MINT n_cols, char *message);

float **FMatrix(MINT n_rows, MINT n_cols, char *message);

RPoint **RMatrix(MINT n_rows, MINT n_cols, char *message);

char **CMatrix(MINT n_rows, MINT n_cols, char *message);

void Gr_start_Timer(time_t *start_time, clock_t *start_clock);

void Gr_stop_Timer(FILE *f, time_t *start_time, clock_t *start_clock);

MFLOAT Frand0();

MFLOAT Frand1();

#endif
