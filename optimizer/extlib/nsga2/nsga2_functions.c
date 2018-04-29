/*========================================================================
  PISA  (www.tik.ee.ethz.ch/pisa/)
  ========================================================================
  Computer Engineering (TIK)
  ETH Zurich
  ========================================================================
  NSGA2

  Implements most functions.

  file: nsga2_functions.c
  author: Marco Laumanns, laumanns@tik.ee.ethz.ch

  revision by: Stefan Bleuler, bleuler@tik.ee.ethz.ch
  last change: $date$
  ========================================================================
*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "nsga2.h"

/* common parameters */
int alpha;  /* number of individuals in initial population */
int mu;     /* number of individuals selected as parents */
int lambda; /* number of offspring individuals */
int dim;    /* number of objectives */


/* local parameters from paramfile*/
int seed;   /* seed for random number generator */
int tournament;  /* parameter for tournament selection */


/* other variables */
char cfgfile[FILE_NAME_LENGTH];  /* 'cfg' file (common parameters) */
char inifile[FILE_NAME_LENGTH];  /* 'ini' file (initial population) */
char selfile[FILE_NAME_LENGTH];  /* 'sel' file (parents) */
char arcfile[FILE_NAME_LENGTH];  /* 'arc' file (archive) */
char varfile[FILE_NAME_LENGTH];  /* 'var' file (offspring) */


/* population containers */
pop *pp_all = NULL;
pop *pp_new = NULL;
pop *pp_sel = NULL;


/* NSGA2 internal global variables */
int *copies;
int **front;
double *dist;



/*-----------------------| initialization |------------------------------*/

void initialize(char *paramfile, char *filenamebase)
/* Performs the necessary initialization to start in state 0. */
{
    FILE *fp;
    int result;
    char str[CFG_ENTRY_LENGTH];

#if (NDEBUG)
    /*
      'result' is used in assert() only. Without below statement
       we get a warning if we compile with NDEBUG (assert() disabled)
    */
    (void)result;
#endif
    /* reading parameter file with parameters for selection */
    fp = fopen(paramfile, "r");
    assert(fp != NULL);

    result = fscanf(fp, "%s", str);
    assert(strcmp(str, "seed") == 0);
    result = fscanf(fp, "%d", &seed);

    result = fscanf(fp, "%s", str);
    assert(strcmp(str, "tournament") == 0);
    result = fscanf(fp, "%d", &tournament); /* fscanf() returns EOF
					       if reading fails. */
    assert(result != EOF); /* no EOF, parameters correctly read */

    fclose(fp);

    srand(seed); /* seeding random number generator */

    sprintf(varfile, "%svar", filenamebase);
    sprintf(selfile, "%ssel", filenamebase);
    sprintf(cfgfile, "%scfg", filenamebase);
    sprintf(inifile, "%sini", filenamebase);
    sprintf(arcfile, "%sarc", filenamebase);

    /* reading cfg file with common configurations for both parts */
    fp = fopen(cfgfile, "r");
    assert(fp != NULL);

    result = fscanf(fp, "%s", str);
    assert(strcmp(str, "alpha") == 0);
    result = fscanf(fp, "%d", &alpha);
    assert(alpha > 0);

    result = fscanf(fp, "%s", str);
    assert(strcmp(str, "mu") == 0);
    result = fscanf(fp, "%d", &mu);
    assert(mu > 0);

    result = fscanf(fp, "%s", str);
    assert(strcmp(str, "lambda") == 0);
    result = fscanf(fp, "%d", &lambda);
    assert(lambda > 0);

    result = fscanf(fp, "%s", str);
    assert(strcmp(str, "dim") == 0);
    result = fscanf(fp, "%d", &dim);
    assert(result != EOF); /* no EOF, 'dim' correctly read */
    assert(dim > 0);

    fclose(fp);

    /* create individual and archive pop */
    pp_all = create_pop(alpha + lambda, dim);
    pp_sel = create_pop(mu, dim);
}


/*-------------------| memory allocation functions |---------------------*/

void* chk_malloc(size_t size)
/* Wrapper function for malloc(). Checks for failed allocations. */
{
    void *return_value = malloc(size);
    if(return_value == NULL)
	PISA_ERROR("Selector: Out of memory.");
    return (return_value);
}


pop* create_pop(int maxsize, int dim)
/* Allocates memory for a population. */
{
    int i;
    pop *pp;

    assert(dim >= 0);
    assert(maxsize >= 0);

    pp = (pop*) chk_malloc(sizeof(pop));
    pp->size = 0;
    pp->maxsize = maxsize;
    pp->ind_array = (ind**) chk_malloc(maxsize * sizeof(ind*));

    for (i = 0; i < maxsize; i++)
	pp->ind_array[i] = NULL;

    return (pp);
}


ind* create_ind(int dim)
/* Allocates memory for one individual. */
{
    ind *p_ind;

    assert(dim >= 0);

    p_ind = (ind*) chk_malloc(sizeof(ind));

    p_ind->index = -1;
    p_ind->fitness = -1;
    p_ind->f = (double*) chk_malloc(dim * sizeof(double));
    return (p_ind);
}


void free_memory()
/* Frees all memory. */
{
   free_pop(pp_sel);
   complete_free_pop(pp_all);
   free_pop(pp_new);
   pp_sel = NULL;
   pp_all = NULL;
   pp_new = NULL;
}


void free_pop(pop *pp)
/* Frees memory for given population. */
{
   if(pp != NULL)
   {
      free(pp->ind_array);
      free(pp);
   }
}


void complete_free_pop(pop *pp)
/* Frees memory for given population and for all individuals in the
   population. */
{
   int i = 0;
   if (pp != NULL)
   {
      if(pp->ind_array != NULL)
      {
         for (i = 0; i < pp->size; i++)
         {
            if (pp->ind_array[i] != NULL)
            {
               free_ind(pp->ind_array[i]);
               pp->ind_array[i] = NULL;
            }
         }

         free(pp->ind_array);
      }

      free(pp);
   }
}


void free_ind(ind *p_ind)
/* Frees memory for given individual. */
{
    assert(p_ind != NULL);

    free(p_ind->f);
    free(p_ind);
}



/*-----------------------| selection functions|--------------------------*/

void selection()
{
    int i;
    int size;

    /* Join offspring individuals from variator to population */
    mergeOffspring();

    size = pp_all->size;

    /* Create internal data structures for selection process */
    /* Vectors */
    copies = (int*) chk_malloc(size * sizeof(int));
    dist = (double*) chk_malloc(size * sizeof(double));

    /* Matrices */
    front = (int**) chk_malloc(size * sizeof(int*));
    for (i = 0; i < size; i++)
    {
	front[i] = (int*) chk_malloc(size * sizeof(int));
    }

    /* Calculates NSGA2 fitness values for all individuals */
    calcFitnesses();

    /* Calculates distance cuboids */
    calcDistances();

    /* Performs environmental selection
       (truncates 'pp_all' to size 'alpha') */
    environmentalSelection();

    /* Performs mating selection
       (fills mating pool / offspring population pp_sel */
    matingSelection();

    /* Frees memory of internal data structures */
    free(copies);
    copies = NULL;
    free(dist);
    dist = NULL;
    for (i = 0; i < size; i++)
    {
	free(front[i]);
    }
    free(front);
    front = NULL;
    return;
}


void mergeOffspring()
{
    int i;

    assert(pp_all->size + pp_new->size <= pp_all->maxsize);

    for (i = 0; i < pp_new->size; i++)
    {
	pp_all->ind_array[pp_all->size + i] = pp_new->ind_array[i];
    }

    pp_all->size += pp_new->size;

    free_pop(pp_new);
    pp_new = NULL;
}


void calcFitnesses()
{
    int i, j, l;
    int size;
    int num;
    int *d;
    int *f;

    size = pp_all->size;
    d = (int*) chk_malloc(size * sizeof(int));
    f = (int*) chk_malloc(size * sizeof(int));

    /* initialize fitness and strength values */
    for (i = 0; i < size; i++)
    {
        pp_all->ind_array[i]->fitness = 0;
	d[i] = 1;
	f[i] = 1;
	copies[i] = 0;
    }

    /* calculate strength values */
    num = size;
    for (l = 0; l < size; l++)
    {
	/* find next front */
	for (i = 0; i < size; i++)
	{
	    d[i] = 0;
	    if (f[i] != 0)
            {
		for (j = 0; j < i && d[i] == 0; j++)
		    if (f[j] != 0)
			if (dominates(pp_all->ind_array[j], pp_all->ind_array[i]))
			    d[i] = 1;
		for(j = i+1; j < size && d[i] == 0; j++)
                    if (f[j] != 0)
                        if (dominates(pp_all->ind_array[j], pp_all->ind_array[i]))
                            d[i] = 1;
	    }
	}

	/* extract front */
	for (i = 0; i < size; i++)
	{
	    if (f[i] != 0 && d[i] == 0)
            {
                pp_all->ind_array[i]->fitness = l;
                f[i] = 0;
                num--;
		front[l][copies[l]] = i;
		copies[l] += 1;
            }
        }

	if (num == 0)
	    break;
    }

    free(d);
    d = NULL;
    free(f);
    f = NULL;
    return;
}


void calcDistances()
{
    int i, j, l, d;
    int size = pp_all->size;
    double dmax = PISA_MAXDOUBLE / (dim + 1);

    for (i = 0; i < size; i++)
    {
	dist[i] = 1;
    }

    for (l = 0; l < size; l++)
    {
	for (d = 0; d < dim; d++)
	{
	    /* sort accorting to d-th objective */
	    for (i = 0; i < copies[l]; i++)
	    {
		int min_index = -1;
		int min = i;
		for (j = i + 1; j < copies[l]; j++)
		{
		    if (pp_all->ind_array[front[l][j]]->f[d] <
			pp_all->ind_array[front[l][min]]->f[d])
			min = j;
		}
		min_index = front[l][min];
		front[l][min] = front[l][i];
		front[l][i] = min_index;
	    }

	    /* add distances */
	    for (i = 0; i < copies[l]; i++)
	    {
		if (i == 0 || i == copies[l] - 1)
		    dist[front[l][i]] += dmax;
		else
		{
		    dist[front[l][i]] +=
			pp_all->ind_array[front[l][i+1]]->f[d] -
			pp_all->ind_array[front[l][i-1]]->f[d];
		}
	    }
	}
    }
}


void environmentalSelection()
{
    int i, j;
    int size = pp_all->size;


    for (i = 0; i < size; i++)
    {
	pp_all->ind_array[i]->fitness += 1.0 / dist[i];
    }

    for (i = 0; i < alpha; i++)
    {
	ind *p_min;
	int min = i;
	for (j = i + 1; j < size; j++)
	{
	    if (pp_all->ind_array[j]->fitness <
		pp_all->ind_array[min]->fitness)
		min = j;
	}
	p_min = pp_all->ind_array[min];
	pp_all->ind_array[min] = pp_all->ind_array[i];
	pp_all->ind_array[i] = p_min;
    }

    for (i = alpha; i < size; i++)
    {
       free_ind(pp_all->ind_array[i]);
       pp_all->ind_array[i] = NULL;
    }

    pp_all->size = alpha;

    return;
}


void matingSelection()
/* Fills mating pool 'pp_sel' */
{
    int i, j;

    for (i = 0; i < mu; i++)
    {
	int winner = irand(pp_all->size);

	for (j = 1; j < tournament; j++)
	{
	    int opponent = irand(pp_all->size);
	    if (pp_all->ind_array[opponent]->fitness
		< pp_all->ind_array[winner]->fitness || winner == opponent)
	    {
		winner = opponent;
	    }
	}
	pp_sel->ind_array[i] = pp_all->ind_array[winner];
    }
    pp_sel->size = mu;
}


void select_initial()
/* Performs initial selection. */
{
    selection();
}


void select_normal()
/* Performs normal selection.*/
{
    selection();
}


int dominates(ind *p_ind_a, ind *p_ind_b)
/* Determines if one individual dominates another.
   Minimizing fitness values. */
{
    int i;
    int a_is_worse = 0;
    int equal = 1;

     for (i = 0; i < dim && !a_is_worse; i++)
     {
	 a_is_worse = p_ind_a->f[i] > p_ind_b->f[i];
          equal = (p_ind_a->f[i] == p_ind_b->f[i]) && equal;
     }

     return (!equal && !a_is_worse);
}


int is_equal(ind *p_ind_a, ind *p_ind_b)
/* Determines if two individuals are equal in all objective values.*/
{
     int i;
     int equal = 1;

     for (i = 0; i < dim; i++)
	 equal = (p_ind_a->f[i] == p_ind_b->f[i]) && equal;

     return (equal);
}


int irand(int range)
/* Generate a random integer. */
{
    int j;
    j=(int) ((double)range * (double) rand() / (RAND_MAX+1.0));
    return (j);
}


/*--------------------| data exchange functions |------------------------*/

int read_ini()
{
    int i;
    pp_new = create_pop(alpha, dim);

    for (i = 0; i < alpha; i++)
	pp_new->ind_array[i] = create_ind(dim);
    pp_new->size = alpha;

    return (read_pop(inifile, pp_new, alpha, dim));
}


int read_var()
{
    int i;

    pp_new = create_pop(lambda, dim);

    for (i = 0; i < lambda; i++)
	pp_new->ind_array[i] = create_ind(dim);

    pp_new->size = lambda;
    return (read_pop(varfile, pp_new, lambda, dim));
}


void write_sel()
{
    write_pop(selfile, pp_sel, mu);
}


void write_arc()
{
     write_pop(arcfile, pp_all, pp_all->size);
}


int check_sel()
{
     return (check_file(selfile));
}


int check_arc()
{
     return (check_file(arcfile));
}