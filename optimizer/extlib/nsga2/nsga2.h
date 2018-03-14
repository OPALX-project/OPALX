/*========================================================================
  PISA  (www.tik.ee.ethz.ch/pisa/)
  ========================================================================
  Computer Engineering (TIK)
  ETH Zurich
  ========================================================================
  NSGA2

  Implementation in C for the selector side.
  
  Header file.
  
  file: nsga2.h
  author: Marco Laumanns, laumanns@tik.ee.ethz.ch

  revision by: Stefan Bleuler, bleuler@tik.ee.ethz.ch
  last change: $date$
  ========================================================================
*/

#ifndef NSGA2_H
#define NSGA2_H

/*-----------------------| specify Operating System |------------------*/
/* necessary for wait() */

/* #define PISA_WIN */
#define PISA_UNIX

/*----------------------------| macro |----------------------------------*/

#define PISA_ERROR(x) fprintf(stderr, "\nError: " x "\n"), fflush(stderr), exit(EXIT_FAILURE)

/*---------------------------| constants |-------------------------------*/
#define FILE_NAME_LENGTH 128 /* maximal length of filenames */
#define CFG_ENTRY_LENGTH 128 /* maximal length of entries in cfg file */
#define PISA_MAXDOUBLE 1E99  /* Internal maximal value for double */

/*----------------------------| structs |--------------------------------*/

typedef struct ind_st  /* an individual */
{
    int index;
    double *f; /* objective vector */
    double fitness;
} ind;

typedef struct pop_st  /* a population */
{
    int size;
    int maxsize;
    ind **ind_array;
} pop;


/*-------------| functions for control flow (in nsga2.c) |------------*/

void write_flag(char *filename, int flag);
int read_flag(char *filename);
void nsga_wait(double sec);

/*---------| initialization function (in nsga2_functions.c) |---------*/

void initialize(char *paramfile, char *filenamebase);

/*--------| memory allocation functions (in nsga2_functions.c) |------*/

void* chk_malloc(size_t size);
pop* create_pop(int size, int dim);
ind* create_ind(int dim);

void free_memory(void);
void free_pop(pop *pp);
void complete_free_pop(pop *pp);
void free_ind(ind *p_ind);

/*-----| functions implementing the selection (nsga2_functions.c) |---*/

void selection();
void mergeOffspring();
void calcFitnesses();
void calcDistances();
int getNN(int index, int k);
double getNNd(int index, int k);
void environmentalSelection();
void truncate_nondominated();
void truncate_dominated();
void matingSelection();

void select_initial();
void select_normal();
int dominates(ind *p_ind_a, ind *p_ind_b);
int is_equal(ind *p_ind_a, ind *p_ind_b);
double calcDistance(ind *p_ind_a, ind *p_ind_b);
int irand(int range);

/*--------------------| data exchange functions |------------------------*/

/* in nsga2_functions.c */

int read_ini(void);
int read_var(void);
void write_sel(void);
void write_arc(void);
int check_sel(void);
int check_arc(void);

/* in nsga2_io.c */

int read_pop(char *filename, pop *pp, int size, int dim);
void write_pop(char *filename, pop *pp, int size);
int check_file(char *filename);

#endif /* NSGA2_H */
