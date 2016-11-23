#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <malloc.h>
#include "mydef.h"
#include "util.h"


/* ------ global variable for checking memory clashes ---- */
#ifdef MEM_CHECK
#define MAX_MEMORY_CHECK  100
MINT            n_memory_check, nr_memory_check = 0;
size_t sizeof_memory[MAX_MEMORY_CHECK];
unsigned MLONG  memory_address[MAX_MEMORY_CHECK];

#endif


/*----------------------------------------------------------------------*/
void error_exit(char *message)
/*----------------------------------------------------------------------*/
{
#ifdef SHOW_RCSID
    char     RCSID[] = "$Id: util.c,v 1.5 1995/08/07 03:37:21 huber Exp huber $";
#endif

    fprintf(stderr, "\n\n\n PROGRAM ENDED WITH ERROR !!!!\n\n");
    fprintf(stderr, "%s\n", message);

    exit(255);
}


/* ----------------------------------------------------------- */
void Get_Line(FILE *f, char *line, MINT *nr_line)
/* ----------------------------------------------------------- */
{
#ifdef SHOW_RCSID
    char     RCSID[] = "$Id: util.c,v 1.5 1995/08/07 03:37:21 huber Exp huber $";
#endif
    char   *check;

    check = fgets(line, MAX_LINE_LEN, f);
    *nr_line += 1;
    while (IS_COMMENT(*line) && (check != 0)) {
        check = fgets(line, MAX_LINE_LEN, f);
        *nr_line += 1;
    }
    if (check == 0) {
        fprintf(stderr, "ERROR READING FILE\n");
        fprintf(stderr, "line number = %i\n", *nr_line);
        exit(255);
    }
}


/*----------------------------------------------------------------------*/
void *Gr_malloc(size_t size, char *message)
/*----------------------------------------------------------------------*/
{
#ifdef SHOW_RCSID
    char     RCSID[] = "$Id: util.c,v 1.5 1995/08/07 03:37:21 huber Exp huber $";
#endif

    void    *new_mem;

    if (size == 0)
        return NULL;

    new_mem = malloc(size);
    if (new_mem == NULL) {
        //fprintf(stderr, "can't allocate enough memory: %d bytes\n", size);
        error_exit(message);
    }
#ifdef MEM_CHECK
    memory_address[nr_memory_check] = (unsigned MLONG) new_mem;
    sizeof_memory[nr_memory_check] =  size;
    fprintf(stderr, "%s : addr = %p, size = %ld\n", message,
           new_mem, sizeof_memory[nr_memory_check]);
    for (n_memory_check = 0; n_memory_check < nr_memory_check; n_memory_check++) {
        if (((memory_address[n_memory_check] < memory_address[nr_memory_check]) &&
             (memory_address[n_memory_check]+sizeof_memory[n_memory_check] >
              memory_address[nr_memory_check])) ||
            ((memory_address[nr_memory_check] < memory_address[n_memory_check]) &&
             (memory_address[nr_memory_check]+sizeof_memory[nr_memory_check] >
              memory_address[n_memory_check])))
            error_exit("Memory clash !!!!!!");
    }
    ++nr_memory_check;
#endif
    return new_mem;
}


/*---------------------NOW A MACRO !! ----------------------------------
void Gr_free(void *old_mem)
----------------------------------------------------------------------
{
#ifdef SHOW_RCSID
    char     RCSID[] = "$Id: util.c,v 1.5 1995/08/07 03:37:21 huber Exp huber $";
#endif

    if (!old_mem)
        free((char *) old_mem);
}
*/

/*--------------------------------------------------------------*/
MINT *Ivector(MINT length, char *message)
/*--------------------------------------------------------------*/
{
#ifdef SHOW_RCSID
    char    RCSID[] = "$Id: util.c,v 1.5 1995/08/07 03:37:21 huber Exp huber $";
#endif
    MINT  *r;

    r = (MINT *) Gr_malloc(length * sizeof(MINT), message);
    return r;
}

/*--------------------------------------------------------------*/
MFLOAT *Dvector(MINT length, char *message)
/*--------------------------------------------------------------*/
{
#ifdef SHOW_RCSID
    char    RCSID[] = "$Id: util.c,v 1.5 1995/08/07 03:37:21 huber Exp huber $";
#endif
    MFLOAT  *r;

    r = (MFLOAT *) Gr_malloc(length * sizeof(MFLOAT), message);
    return r;
}


/*--------------------------------------------------------------*/
float *Fvector(MINT length, char *message)
/*--------------------------------------------------------------*/
{
#ifdef SHOW_RCSID
    char    RCSID[] = "$Id: util.c,v 1.5 1995/08/07 03:37:21 huber Exp huber $";
#endif
    float  *r;

    r = (float *) Gr_malloc(length * sizeof(float), message);
    return r;
}

/*--------------------------------------------------------------*/
RPoint *Rvector(MINT length, char *message)
/*--------------------------------------------------------------*/
{
#ifdef SHOW_RCSID
    char    RCSID[] = "$Id: util.c,v 1.5 1995/08/07 03:37:21 huber Exp huber $";
#endif
    RPoint  *r;

    r = (RPoint *) Gr_malloc(length * sizeof(RPoint), message);
    return r;
}

/*--------------------------------------------------------------*/
MINT **IMatrix(MINT n_rows, MINT n_cols, char *message)
/*--------------------------------------------------------------*/
{
#ifdef SHOW_RCSID
    char    RCSID[] = "$Id: util.c,v 1.5 1995/08/07 03:37:21 huber Exp huber $";
#endif
    MINT   i;
    MINT   **matrix;

    matrix = (MINT **) Gr_malloc(n_rows*sizeof(MINT *), message);
    matrix[0] = (MINT *) Gr_malloc(n_rows*n_cols*sizeof(MINT), message);

    for (i = 1; i < n_rows; i++)
        matrix[i] = matrix[i-1] + n_cols;
    return matrix;
}


/*--------------------------------------------------------------*/
MFLOAT **DMatrix(MINT n_rows, MINT n_cols, char *message)
/*--------------------------------------------------------------*/
{
#ifdef SHOW_RCSID
    char    RCSID[] = "$Id: util.c,v 1.5 1995/08/07 03:37:21 huber Exp huber $";
#endif
    MINT   i;
    MFLOAT **matrix;

    matrix = (MFLOAT **) Gr_malloc(n_rows*sizeof(MFLOAT *), message);
    matrix[0] = (MFLOAT *) Gr_malloc(n_rows*n_cols*sizeof(MFLOAT), message);

    for (i = 1; i < n_rows; i++)
        matrix[i] = matrix[i-1] + n_cols;
    return matrix;
}

/*--------------------------------------------------------------*/
float **FMatrix(MINT n_rows, MINT n_cols, char *message)
/*--------------------------------------------------------------*/
{
#ifdef SHOW_RCSID
    char    RCSID[] = "$Id: util.c,v 1.5 1995/08/07 03:37:21 huber Exp huber $";
#endif
    MINT   i;
    float **matrix;

    matrix = (float **) Gr_malloc(n_rows*sizeof(float *), message);
    matrix[0] = (float *) Gr_malloc(n_rows*n_cols*sizeof(float), message);

    for (i = 1; i < n_rows; i++)
        matrix[i] = matrix[i-1] + n_cols;
    return matrix;
}



/*--------------------------------------------------------------*/
RPoint **RMatrix(MINT n_rows, MINT n_cols, char *message)
/*--------------------------------------------------------------*/
{
#ifdef SHOW_RCSID
    char    RCSID[] = "$Id: util.c,v 1.5 1995/08/07 03:37:21 huber Exp huber $";
#endif
    MINT   i;
    RPoint **matrix;

    matrix = (RPoint **) Gr_malloc(n_rows*sizeof(RPoint *), message);
    matrix[0] = (RPoint *) Gr_malloc(n_rows*n_cols*sizeof(RPoint), message);

    for (i = 1; i < n_rows; i++)
        matrix[i] = matrix[i-1] + n_cols;
    return matrix;
}


/*--------------------------------------------------------------*/
char **CMatrix(MINT n_rows, MINT n_cols, char *message)
/*--------------------------------------------------------------*/
{
#ifdef SHOW_RCSID
    char    RCSID[] = "$Id: util.c,v 1.5 1995/08/07 03:37:21 huber Exp huber $";
#endif
    MINT   i;
    char   **matrix;

    matrix = (char **) Gr_malloc(n_rows*sizeof(char *), message);
    matrix[0] = (char *) Gr_malloc(n_rows*n_cols*sizeof(char), message);

    for (i = 1; i < n_rows; i++)
        matrix[i] = matrix[i-1] + n_cols;
    return matrix;
}


/*--------------------------------------------------------------*/
void Gr_start_Timer(time_t *start_time, clock_t *start_clock)
/*--------------------------------------------------------------*/
{
#ifdef SHOW_RCSID
    char    RCSID[] = "$Id: util.c,v 1.5 1995/08/07 03:37:21 huber Exp huber $";
#endif

    *start_time = time(NULL);
    *start_clock = clock();
}


/*--------------------------------------------------------------*/
void Gr_stop_Timer(FILE *f, time_t *start_time, clock_t *start_clock)
/*--------------------------------------------------------------*/
{
#ifdef SHOW_RCSID
    char    RCSID[] = "$Id: util.c,v 1.5 1995/08/07 03:37:21 huber Exp huber $";
#endif
    time_t      end_time;
    clock_t     end_clock;
    MFLOAT      used_time, used_clock;
    MINT        days, hours, minutes, seconds;



    end_clock = clock();
    end_time = time(NULL);
    fprintf(f, "\n\n---------------- time measurement ---------------\n");
    fprintf(f, "start of simulation : %s",ctime(start_time));
    fprintf(f, "end of simulation : %s",ctime(&end_time));
#ifdef CLOCKS_PER_SEC
    used_time = difftime(end_time, *start_time);
    seconds = (MINT) used_time;
    days = used_time / 86400;
    seconds -= days*86400;
    hours = seconds / 3600;
    seconds -= hours * 3600;
    minutes = seconds / 60;
    seconds -= minutes * 60;
    
    fprintf(f, "used real time:\n");
    fprintf(f, "%i days %i hours %i minutes %i seconds\n", days, hours, minutes, seconds);

    used_clock = difftime(end_clock, *start_clock) / CLOCKS_PER_SEC;
    seconds = (MINT) used_clock;
    used_clock -= seconds;
    days = seconds / 86400;
    seconds -= days*86400;
    hours = seconds / 3600;
    seconds -= hours * 3600;
    minutes = seconds / 60;
    seconds -= minutes * 60;
    used_clock += seconds;
    fprintf(f, "used cpu time:\n");
    fprintf(f, STR(%i days %i hours %i minutes %6.2f seconds\n), days, hours, minutes, used_clock);
#endif
}

/*--------------------------------------------------------------*/
MFLOAT Frand0()
/*--------------------------------------------------------------*/
{
    MFLOAT r;
    r = rand();
    r /= (MFLOAT) RAND_MAX;
    return r;
}

/*--------------------------------------------------------------*/
MFLOAT Frand1()
/*--------------------------------------------------------------*/
{
    MFLOAT r;
    r = rand();
    r /= (MFLOAT) RAND_MAX;
    r -= 0.5;
    r *= 2;
    return r;
}
