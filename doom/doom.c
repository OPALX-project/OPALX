/* #define _TREE_ */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <ctype.h>
#include <time.h>
#ifdef _CERN_
#include "/afs/cern.ch/user/h/hansg/public/doom/dev/doomex.h"
#include "/afs/cern.ch/user/h/hansg/public/doom/dev/doom.h"
#endif
#ifndef _CERN_
#include "doomex.h"
#include "doom.h"
#endif

/*---------------------------------------------------------------------*
*                                                                      *
*                           CERN                                       *
*                                                                      *
*     European Organization for Nuclear Research                       *
*                                                                      *
*     Program name: DOOM: Mad Object Oriented Database                 *
*                                                                      *
*     Author and contact:   Hans GROTE                                 *
*                           SL Division                                *
*                           CERN                                       *
*                           CH-1211 GENEVA 23                          *
*                           SWITZERLAND                                *
*                      Tel. [041] (022) 767 49 61                      *
*                           Hans.Grote@cern.ch                         *
*                                                                      *
*     Copyright  CERN,  Geneva  1990  -  Copyright  and  any   other   *
*     appropriate  legal  protection  of  this  computer program and   *
*     associated documentation reserved  in  all  countries  of  the   *
*     world.                                                           *
*                                                                      *
*     Organizations collaborating with CERN may receive this program   *
*     and documentation freely and without charge.                     *
*                                                                      *
*     CERN undertakes no obligation  for  the  maintenance  of  this   *
*     program,  nor responsibility for its correctness,  and accepts   *
*     no liability whatsoever resulting from its use.                  *
*                                                                      *
*     Program  and documentation are provided solely for the use  of   *
*     the organization to which they are distributed.                  *
*                                                                      *
*     This program  may  not  be  copied  or  otherwise  distributed   *
*     without  permission. This message must be retained on this and   *
*     any other authorized copies.                                     *
*                                                                      *
*     The material cannot be sold. CERN should be  given  credit  in   *
*     all references.                                                  *
*                                                                      *
*---------------------------------------------------------------------*/

/* element attribute names and positions */

const int el_att_num = 70;

const int el_att_pos[] = {
10, 7, 16, 2, 3, 2, 3, 9, 4, 4,
5, 11, 8, 12, 13, 14, 32, 33, 34, 35,
36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
46, 47, 48, 49, 50, 51, 15, 16, 52, 53,
17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
27, 28, 29, 30, 31, 7, 0, 5, 11, 17,
8, 1, 9, 14, 15, 10, 6, 1, 12, 13
};

char* el_att_names[] = {
"ANGLE","BETRF","CHARGE","E1","E2","EX","EY","FINT","FREQ","H1",
"H2","HARMON","HGAP","K0","K0S","K1","K10","K10S","K11","K11S",
"K12","K12S","K13","K13S","K14","K14S","K15","K15S","K16","K16S",
"K17","K17S","K18","K18S","K19","K19S","K1S","K2","K20","K20S",
"K2S","K3","K3S","K4","K4S","K5","K5S","K6","K6S","K7",
"K7S","K8","K8S","K9","K9S","KS","L","LAG","LRAD","NPART",
"PG","RHOINV","SHUNT","SIGX","SIGY","TFILL","TILT","VOLT","XSIZE","YSIZE"
};

/* optics summary attribute names and positions */

const int opt_att_num = 17;

const int opt_att_pos[] = {
1, 11, 12, 0, 13, 15, 14, 16, 2, 3, 
4, 5, 6, 9, 7, 10, 8
};

char* opt_att_names[] = {
"ALFA","BXMAX","BYMAX","DELTAP","DXMAX","DXRMS",
"DYMAX","DYRMS","GAMMATR","QX",
"QY","XIX","XIY","XMAX","XRMS","YMAX","YRMS"
};

/*further declarations */

static FILE* db_stream = NULL;

static int curr_off = 0;                 /* current file offset */
#ifdef _TREE_
static int debug_flag = 3;               /* debug print levels */
#endif
#ifndef _TREE_
static int debug_flag = 0;               /* debug print levels */
#endif
static int force_old  = 0;               /* if not 0, require existing d.b. */
static int twiss_count = 0;              /* delta counter for twiss tables */
static int newdb_flag = 0;               /* if 1 new d.b., else 0 */
static int swap_flag = 0;                /* if != 0 swap */
static int no_update  = 0;               /* if 1, d.b. not updated */
static int prev_kept  = 0;               /* no. of old objects kept */
static int prev_freed = 0;               /* no. of old objects freed */

static const int buff_size = 10000000;    /* output buffer size  */
static const int max_alloc = 20000000;   /* maximum single memory 
                                            allocation in char.s  */
static const int c_table_size = 10000;   /* length of table work array */
static const int o_list_numb = 2;        /* number of system lists */
static const int o_list_count = 20;      /* number of system sublists  */

static const int major_version = 1;
static const int minor_version = 15;
static const int time_vault = 852076800;  /* sec.s from 1.1.1970 to 1.1.1997 */
static const double micro_step = 0.00001; /* time microstep */

static int major_time = 0;                /* time in sec.s since 1.1.1997 */
static int minor_time = 0;                /* used for time microsteps */
static int reopen = 0;                    /* 1 if d.b. re-opened */
static char info[] = "info-00-00";        /* d.b. info - directly saved */
static int info_int[INFO_BASE + N_SYSTEM * INFO_STEP];
static const int info_size = INFO_BASE + N_SYSTEM * INFO_STEP;
/*
   info_int contains the basic system information.

   0 = length of info_int in bytes
   1 = major version number (must agree between program and d.b.file)
   2 = minor version number (#/10 must agree between program and d.b.file) 
   3 = major time
   4 = minor time
   5 = current dynamic code

  INFO_BASE+0 .. +7: offset, length (bytes), int, double, char size,
                                    int, double, char occupation
            for system bank zero
  INFO_BASE+INFO_STEP+0 .. +7: ditto system bank one etc.
*/

static char* system_0[] = 
       {"system-00", /* list of objects - int list: index to list */
        "system-01", /* int list: major time of last object load from d.b.*/
        "system-02", /* int list: minor time of last object load from d.b.*/
        "system-03", /* int list: file offset in bytes */
        "system-04"};/* int list: size in bytes */
static char system_10[] = "system-10"; /* list of environment variables */
                       /* int list: index to list */
static char system_11[] = "system-11"; /* list of environment values */
                       /* int list: index to list */
static char *pc,                         /* utility character pointer  */
       *buffer;                          /* input/output buffer  */

static char *cdummy, *char_tok;
static int idummy[1]; 
static double* ddummy;

static struct object*** o_list_sys;      /* system lists */

/* internal set variable list */
const int n_int_var = 3;
char* l_int_var[] = {
"debug_flag", "force_old", "twiss_count",
};
int* p_int_var[] = {&debug_flag, &force_old, &twiss_count};

/* environment variable list and defaults */
const int environ_size = 4;
char* environ_vars[] = {
"SEQUENCE", "DUMMY", "USE", "SUB_TABLE"
};
char* environ_defs[] = {
"free", "DUMMY", "UNNAMED_USE", "none"
};

/* start of code  */

void clean_obj(struct object* p) /* sets integer, d.p., pointers to zero,
                                    and sets their counter to zero */
{
  int j;

  for (j = 0; j < p->l_int; j++)  p->a_int[j] = 0;
  p->c_int = 0;
  for (j = 0; j < p->l_dble; j++)  p->a_dble[j] = 0;
  p->c_dble = 0;
  for (j = 0; j < p->l_obj; j++)  p->p_obj[j] = NULL;
  p->c_obj = 0;
}

int db_close()
{
  int ret;

  ret = fclose(db_stream); db_stream = NULL;
  return ret;
}

int db_exists(char* lkey) /*  */
{
#ifdef _TREE_
       puts("++++++ enter db_exists");
#endif
  return (n_list_pos(lkey, o_list_sys[0][0]) < 0 ? 1 : 0);
}

void db_get /* gets one record from the database - returns 0 failure, 1 OK */
            (char* lkey,      /* input: full record key */
             int* counts,     /* input: max. length of control, integer,
                                 d.p., char arrays;
                                 output: actual occupations */
             int* a_cont,     /* output: control array */
             int* a_int,      /* output: integer array */
             double* a_dble,  /* output: d.p. array */
             char* a_char)    /* output: string */
{
  struct object* l00 = o_list_sys[0][0];
  struct object* l03 = o_list_sys[0][3];
  struct object* l04 = o_list_sys[0][4];
  int lp, index, err, offset, size;

#ifdef _TREE_
       puts("++++++ enter db_get");
#endif
  lp = n_list_pos(lkey, l00);
  if (lp >= 0)
    {
     printf("<<< DOOM >>> fatal - db_get undefined key: %s\n", lkey);
     exit(1);
    }
  index = l00->a_int[-lp-1];
  offset = l03->a_int[index];
  size = l04->a_int[index];
  if ((err = get_text(offset, size, buffer)) != size)
   {
    if (err < 0)  printf("DOOM fatal - error %d when reading record: %s\n", 
                         err, lkey);
    else printf("DOOM fatal - size %d (actual %d) when reading record: %s\n", 
                err, size, lkey);
    exit(1);
   }
  rec_get(counts, a_cont, a_int, a_dble, a_char);
}

void db_get_info()    /* gets info record from d.b. - always at offset 0 */
{
  int j, err, size;

#ifdef _TREE_
       puts("++++++ enter db_get_info");
#endif
  size = INT_SIZE * info_size;
  if ((err = get_text(0, size, (char*) info_int)) != size) 
   {
    if (err < 0)  printf("DOOM fatal - error %d when reading record: %s\n", 
                         err, info);
    else  printf("DOOM fatal - size %d (actual %d) when reading record: %s\n",
                 err, size, info);
    printf("           - not a DOOM database ?");
    exit(1);
   }
  if(info_int[0] != size)
   {
    swap_flag = 1; swap_int(info_size, info_int);
    if(info_int[0] != size)
     {
      printf("DOOM fatal - wrong length %d in record: %s\n", err, info);
      exit(1);
     }
   }
  if (info_int[1] != major_version || info_int[2]/10 != minor_version/10)
   {
    printf("<<< DOOM >>> fatal - program   version number = %d.%d\n",
    major_version, minor_version);
    printf("                     d.b. file version number = %d.%d\n",
    info_int[1], info_int[2]);
    exit(1);
   }
  major_time = info_int[3]; minor_time = info_int[4];
  if (debug_flag > 3)
    {
     puts("+++ info record +++");
     for (j = 0; j < info_size; j++)
       {
        if (j%10 == 0) printf("\n");
        printf("%d%s", info_int[j], "  ");
       }
     printf("\n");
    }
}

void db_get_system()    /* gets system records from d.b. -  */
{
  int j, k, size, aux, err, counts[] = {CONT_OBJ, 0, 0, 0};
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter db_get_system");
#endif
  for (k = 0; k < N_SYSTEM; k++)
    {
     aux = INFO_BASE + INFO_STEP * k;
     p = o_list_sys[0][k] = make_obj(system_0[k],
                        info_int[aux+2], info_int[aux+3],
                        info_int[aux+4], info_int[aux+2]);
     counts[1] = info_int[aux+5]; 
     counts[2] = info_int[aux+6]; 
     counts[3] = info_int[aux+7];
     size = db_length(counts);
     if (info_int[aux+1] != size)
      {
       printf("DOOM fatal - system bank %d length %d should be %d\n",
             k, info_int[aux+1], size);
       exit(1);
      }
     if ((err = get_text(info_int[aux], size, buffer)) != size) 
      {
       printf("DOOM fatal - db_get size error reading record: %s\n", 
              system_0[k]);
       printf("           - size read = %d should be %d\n",
              err, size);
       exit(1);
      }
     rec_get(counts, &p->ma_time, p->a_int, p->a_dble, p->a_char);
     for (j = 0; j < p->l_obj; j++) p->p_obj[j] = NULL;
     if (k == 0)  make_names(p);
    }
  curr_off = info_int[INFO_BASE];
}

int db_is_open() {return (db_stream == NULL ? 0 : 1);}

int db_length(int* counts)
{
  return (INT_SIZE * (REC_CONT + counts[0] + counts[1]) 
          + DOUBLE_SIZE * counts[2] + counts[3]);
}

void db_list()
{
}

int db_open(char* filename, int flag)
{
  size_t l_written, size;

#ifdef _TREE_
       puts("++++++ enter db_open");
#endif
  size = INT_SIZE * info_size;
  if (flag == 0)  /* open existing d.b. */
    {
     if ((db_stream = fopen(filename, "r+")) == NULL)  return 0;
    }
  else            /* create new d.b. */
    {
     if ((db_stream = fopen(filename, "w+")) == NULL)  return 0;
     if (fseek(db_stream, 0, SEEK_SET) != 0)  return 0;
     size = INT_SIZE * info_size;
     if ((l_written = fwrite((char*) info_int, 1L, size, db_stream)) 
	 < size) return 0; /* just fill the space for later */
     curr_off = INT_SIZE * info_size;
    }
  return 1;
}

int db_put /* puts one record into the database - returns 0 failure, 1 OK */
            (char* lkey,     /* input: full record key */
             int* counts,    /* input: length of control, integer,
                                 d.p., char arrays */
             int* a_cont,    /* input: control array */
             int* a_int,     /* input: integer array */
             double* a_dble, /* input: d.p. array */
             char* a_char)   /* input: string */
{
  int lp, index, size;
  struct object* l00 = o_list_sys[0][0];
  struct object* l03 = o_list_sys[0][3];
  struct object* l04 = o_list_sys[0][4];

#ifdef _TREE_
       puts("++++++ enter db_put");
#endif
  if (no_update == 0)
    {
     lp = n_list_pos(lkey, l00);
     if (lp >= 0)
       {
        printf("<<< DOOM >>> fatal - db_put undefined key: %s", lkey);
        exit(1);
       }
     index = l00->a_int[-lp-1];
     size = db_length(counts);
     if (size > buff_size)
       {
        printf("%s%d\n", "DOOM fatal - db_put buffer overflow at: ",
               buff_size);
        printf("%s%s\n", "             record key: ", lkey);
        exit(1);
       }
     if (l03->a_int[index] == 0 
         || l04->a_int[index] < size ) /* requires new record */
       {
        l04->a_int[index] = size;
        l03->a_int[index] = curr_off;
        curr_off += size;
       }
     rec_put(counts, a_cont, a_int, a_dble, a_char);
     return (put_text(l03->a_int[index], size, buffer) == size ? 1 : 0);
    }
  else return 0;
}

void db_put_info()    /* puts info record into d.b. -  */
{
  int size;

#ifdef _TREE_
       puts("++++++ enter db_put_info");
#endif
  size = INT_SIZE * info_size;
  info_int[0] = size;
  info_int[1] = major_version; info_int[2] = minor_version;
  info_int[3] = major_time; info_int[4] = minor_time;
  if (swap_flag) swap_int(info_size, info_int);
  if (put_text(0, size, (char*) info_int) != size) 
   {
    printf("%s%s\n",
    "DOOM fatal - db_put failure to store record: ", info);
    exit(1);
   }
}

void db_put_system()    /* puts system records into d.b. -  */
{
  int k, size, aux, err, counts[] = {CONT_OBJ, 0, 0, 0};
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter db_put_system");
#endif
  for (k = 0; k < N_SYSTEM; k++)
    {
     aux = INFO_BASE + INFO_STEP * k;
     p = o_list_sys[0][k];
     counts[1] = p->c_int; counts[2] = p->c_dble; counts[3] = p->c_char;
     size = db_length(counts);
     info_int[aux]   = curr_off;
     info_int[aux+1] = size;
     info_int[aux+2] = p->l_int; 
     info_int[aux+3] = p->l_dble; 
     info_int[aux+4] = p->l_char;
     info_int[aux+5] = p->c_int; 
     info_int[aux+6] = p->c_dble; 
     info_int[aux+7] = p->c_char;
     rec_put(counts, &p->ma_time, p->a_int, p->a_dble, p->a_char);
     if ((err = put_text(curr_off, size, buffer)) != size) 
      {
       printf("%s%s\n",
       "DOOM fatal - db_put failure to store record: ", system_0[k]);
       exit(1);
      }
     curr_off += size;
    }
}

void delete_obj(struct object* p)  /* deletes object + sys entry */
{
  int lp, index;

#ifdef _TREE_
       puts("++++++ enter delete_obj");
#endif
  if (p != NULL)
    {
     lp = n_list_pos(p->key, o_list_sys[0][0]);
     if (p->a_int != NULL)  free(p->a_int);
     if (p->a_dble != NULL)  free(p->a_dble);
     if (p->a_char != NULL)  free(p->a_char);
     if (p->p_obj != NULL)  free(p->p_obj);
     if (p->names != NULL)  free(p->names);
     free(p); p = NULL;
     if (lp < 0)  /* object exists in list */
       {
        index = o_list_sys[0][0]->a_int[-lp-1];
        o_list_sys[0][0]->p_obj[index] = NULL;
       }
    }
}

void dmp_i_list(int c, int* l) 
{
  int j;
  for (j = 0; j < c; j++) printf("%d\n", l[j]);
}

/* C start doom_close: Close database */
void doom_close()
/* C end comment */
/* F start doom_close

subroutine doom_close

Closes the database currently open.

F end comment */
{
  int k, j, tma, tmi, index, saved = 0;
  struct object *p, *lsys[N_SYSTEM];

#ifdef _TREE_
       puts("++++++ enter doom_close");
#endif
  if (db_is_open() == 0)
    {
     puts("+++DOOM - doom_close: db not open"); return;
    }
  else if (no_update == 0 
           && o_list_sys[0][0] != NULL) 
    {
     for (k = 0; k < N_SYSTEM; k++) lsys[k] = o_list_sys[0][k];
     for (j = 0; j < lsys[0]->c_obj; j++)
       {
        index = lsys[0]->a_int[j];
        p = lsys[0]->p_obj[index];
        for (k = 0; k < N_SYSTEM; k++)
	  {
           if (p == lsys[k]) p = NULL;
	  }
        if (p != NULL)
	  {
           tma = p->ma_time; tmi = p->mi_time;
           if(tma > lsys[1]->a_int[index] ||
              (tma == lsys[1]->a_int[index] && tmi  > lsys[2]->a_int[index]))
	     {
	      lsys[1]->a_int[index] = tma; lsys[2]->a_int[index] = tmi;
              put_obj(p);
              saved++;
             }
	  }
       }
     db_put_system();
     db_put_info(); /* must be last - time + system inf. saved there */
     printf("\n%s%d\n\n", " <<< DOOM >>> number of objects saved: ", saved);
    }
   if (db_close() != 0)
    {
     puts("+++DOOM - doom_close: closure failure");
    }
}

void doom_blow                 /* replaces '\0' in string by blanks */
            (char* a_char,     /* input/output string */
             int*  l_in,       /* input length */
             int*  l_item,     /* length of one item */
             int*  l_out)      /* output length in multiples of l_item */
{
  int i, j, k, k0, n = 0, loc = *l_in;

#ifdef _TREE_
       puts("++++++ enter doom_blow");
#endif
  if (a_char[loc-1] != '\0')
    {
     puts("+++doom_blow fatal - bad compressed string");
     exit(1);
    }
  for (i = 0; i < loc; i++)
    {
     if (a_char[i] == '\0')  n++;
    }
  *l_out = n; k0 = n * *l_item; k = k0;
  for (i = loc-2; i>=0; i--)
    {
     if (a_char[i] == '\0')
       {
        n = k0 - k; k0 = k0 - *l_item;
        for (j = 0; j < n; j++)
	  {
           a_char[k0+j] = a_char[k++];
          }
        for (j = n; j < *l_item; j++)
	  {
           a_char[k0+j] = ' ';
	  }
        k = k0;
       }
     else    a_char[--k] = a_char[i];
    }
  n = k0 - k; k0 = k0 - *l_item;
  if (k0 != 0)
    {
     printf("%s%d\n", "doom_blow fatal - internal error - k0 = ", k0);
     exit(1);
    }
  for (j = 0; j < n; j++)
    {
     a_char[j] = a_char[k++];
    }
  for (j = n; j < *l_item; j++)
    {
     a_char[j] = ' ';
    }
}

void doom_compress             /* replaces blanks in string by one '\0' */
            (char* a_char,     /* input/output string */
             int*  l_in,       /* input length */
             int*  l_out)      /* output length */
{
  int i, loc = 0, k = 0;

#ifdef _TREE_
       puts("++++++ enter doom_compress");
#endif
  for (i = 0; i < *l_in; i++)
    {
     if (a_char[i] != ' ' && a_char[i] != '\x9')
       {
        a_char[loc++] = a_char[i]; k = 0;
       }
     else if (k == 0)
       {
        a_char[loc++] = '\0'; k = 1;
       }
    }
  if (a_char[loc-1] != '\0')  a_char[loc++] = '\0';
  *l_out = loc;
}

void doom_comp_max             /* replaces blanks in string by one '\0',
                                  inserts as well at max. count/item */
            (char* a_char,     /* input/output string */
             int*  l_in,       /* input length */
             int*  l_item,     /* max. item length */
             int*  l_out)      /* output length */
{
  int i, c = 0, loc = 0, k = 0;

#ifdef _TREE_
       puts("++++++ enter doom_comp_max");
#endif
  for (i = 0; i < *l_in; i++)
    {
     if (a_char[i] != ' ' && a_char[i] != '\x9')
       {
        a_char[loc++] = a_char[i]; k = 0; c++;
       }
     else
       {
        if (k++ == 0)
          {
           a_char[loc++] = '\0'; c = 0;
          }
        else if (c == *l_item)
          {
           a_char[loc++] = '\0'; c = 0;
          }
       }
    }
  if (a_char[loc-1] != '\0')  a_char[loc++] = '\0';
  *l_out = loc;
}

void doom_dump(  /*  dumps all objects starting with a letter */
              int* flag) /* 0: names, 1: partial, 2: full */
{
  int j, index;
  struct object* p=0;

  printf("+++ DOOM d.b. dump +++  objects total: %d\n\n", 
         o_list_sys[0][0]->c_obj);
  for (j = 0; j < o_list_sys[0][0]->c_obj; j++)
    {
     index = o_list_sys[0][0]->a_int[j];
     if (o_list_sys[0][0]->p_obj[index] == NULL)
        p = doom_fetch(o_list_sys[0][0]->names[index]);
     if (p != NULL)
       {
              puts("+++ ");
        if (*flag == 0) 
           printf("%d  %s\n", j, p->key);
        else if(*flag == 1)
           prt_obj(p);
        else if(*flag > 1)
           prt_obj_full(p);
       }
    }
}


/* C start doom_dyncode: produce record code dynamically */
int doom_dyncode()          /* returns each time a new integer code */
/* C end comment */
/* F start doom_dyncode

integer function doom_dyncode()

Returns each time a new integer code

F end comment */
{
  return ++info_int[5];
}

/* C start doom_exist: Does object exist ? */
int doom_exist              /* checks existence in list and d.b.
                               returns -1 not found, 0 not mod., 1 mod. */
            (char* name,    /* input: (full) d.b. object name */
             int* key_list) /* input: key list for composed name */
/* C end comment */
/* F start doom_exist

integer function doom_exist(name, keylist)

Checks the existence of an object in the database.

Input        (character)          name
Input        (integer array)      keylist

Returns -1 if object not found, 0 if not modified, 1 if modified

F end comment */
{
  char rec_key[KEY_LENGTH], lkey[KEY_LENGTH];
  int lp, index;

#ifdef _TREE_
       puts("++++++ enter doom_exist");
#endif
  mycpy(lkey, name);
  make_key(lkey, key_list, rec_key);
  lp = n_list_pos(rec_key, o_list_sys[0][0]);
  if (lp < 0)  /* object exists */
    {
     index = o_list_sys[0][0]->a_int[-lp-1];
     return ((o_list_sys[0][1]->a_int[index] < 
              o_list_sys[0][0]->p_obj[index]->ma_time ||
              (o_list_sys[0][1]->a_int[index] == 
               o_list_sys[0][0]->p_obj[index]->ma_time &&
               o_list_sys[0][2]->a_int[index] < 
               o_list_sys[0][0]->p_obj[index]->mi_time)
              ) ? 1 : 0
             );
    }
  else
    {
     return -1;
    }
}

/* C start doom_exist_env: Does env$object exist ? */
int doom_exist_env          /* checks existence in list and d.b.
                               returns -1 not found, 0 not mod., 1 mod. */
            (char* env,     /* input: precede with environment value+_ */
             char* name,    /* input: (full) d.b. object name */
             int* key_list) /* input: key list for composed name */
/* C end comment */
/* F start doom_exist_env

integer function doom_exist_env(env, name, keylist)

Checks the existence of an object in the database with preceding env label

Input        (character)          env      environment label preceding name
Input        (character)          name
Input        (integer array)      keylist

Returns -1 if object not found, 0 if not modified, 1 if modified

F end comment */
{
  char rec_key[KEY_LENGTH], t_key[KEY_LENGTH], lkey[KEY_LENGTH],
       env_name[KEY_LENGTH], env_value[KEY_LENGTH];
  int lp, index;

#ifdef _TREE_
       puts("++++++ enter doom_exist_env");
#endif
  mycpy(lkey, name);
  if (*env != ' ' && *env != '\0')
    {
     mycpy(env_name, env);
     doom_getenv(env_name, env_value);
     make_key(lkey, key_list, t_key);
     make_dollar_key(env_value, t_key, rec_key);
    }
  else make_key(lkey, key_list, rec_key);
  lp = n_list_pos(rec_key, o_list_sys[0][0]);
  if (lp < 0)  /* object exists */
    {
     index = o_list_sys[0][0]->a_int[-lp-1];
     return ((o_list_sys[0][1]->a_int[index] < 
              o_list_sys[0][0]->p_obj[index]->ma_time ||
              (o_list_sys[0][1]->a_int[index] == 
               o_list_sys[0][0]->p_obj[index]->ma_time &&
               o_list_sys[0][2]->a_int[index] < 
               o_list_sys[0][0]->p_obj[index]->mi_time)
             ) ? 1 : 0
            );
    }
  else
    {
     return -1;
    }
}

/* C start doom_fetch: Fetch an object */
struct object* doom_fetch( /* fetches complete object, enters in list */
       char* rec_key)      /* input: full object key */
/* C end comment */
{
  struct object* p = NULL;
  char lkey[KEY_LENGTH];
  int lp, index;

#ifdef _TREE_
       puts("++++++ enter doom_fetch");
#endif
  mycpy(lkey, rec_key);
  lp = n_list_pos(lkey, o_list_sys[0][0]);
  if (lp < 0)  /* object exists */
    {
     index = o_list_sys[0][0]->a_int[-lp-1];
     p = o_list_sys[0][0]->p_obj[index];
     if (p == NULL)
       {
        if (db_is_open() == 0)
          {
           puts("+++DOOM - doom_fetch: db not open"); 
          }
        else
          {
           if (db_exists(lkey))  p = get_obj(lkey);
           else
	     {
              printf
              (" <<< DOOM >>> fatal: object %s in list, not in d.b.\n", 
               rec_key);
              exit(1);
	     }
           o_list_sys[0][0]->p_obj[index] = p;
           o_list_sys[0][1]->a_int[index] = p->ma_time;
           o_list_sys[0][2]->a_int[index] = p->mi_time;
	  }
       }
    }
  return p;
}

/* C start doom_galign: Get alignment errors */
void doom_galign         /* gets alignment errors for one element */
        (char* name,     /* input: element name */
         int*  occ,      /* input: occurrence count */
         int*  nerr,     /* input/output: max/returned number of errors,
                            or = -1 if element not found */
         double* al_err) /* alignment errors */
/* C end comment */
/* F start doom_galign

subroutine doom_galign(name, occur, nerr, errors)

Returns the alignment errors for one element in a used sequence

Input        (character)          name
Input        (integer)            occur(rence count)
Input/Output (integer)     nerr - on input maximum number of errors requested
                                - on output number of errors returned,
                                  or = -1 if errors not found
Output (d.p. array)               (alignment) errors

F end comment */
{
  int j, key_list[3] = {2, ALIGN_CODE, 0};
  char rec_key[KEY_LENGTH], t_key[KEY_LENGTH], lkey[KEY_LENGTH], 
       sequ_name[KEY_LENGTH], 
       type[] = "ALIGN_ERROR";
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter doom_galign");
#endif
  key_list[2] = *occ;
  mycpy(lkey, name);
  doom_getenv("USE", sequ_name);
  make_key(lkey, key_list, t_key);
  make_dollar_key(sequ_name, t_key, rec_key);
  p = doom_fetch(rec_key);
  if (p == NULL) *nerr = -1;
  else if (strcmp(p->obj_type, type) != 0)
    {
     printf("%s%s\n%s%s%s%s\n", " <<< DOOM >>> warning: object ", rec_key,
     " has type: ", p->obj_type, " and not type: ", type);
     *nerr = -1;
    }
  else
    {
     if (*nerr > p->c_dble)  *nerr = p->c_dble;
     for (j = 0; j < *nerr; j++) al_err[j] = p->a_dble[j];
    }
}

/* C start doom_gbase: Get object basic type name */
void doom_gbase       /* returns basic type of object */
     (char* name,     /* input: element name */
      int* key_list,  /* input: key_list */
      int* nitem,     /* input: length of type character variable */
      char* type,     /* output: base type name */
      int* ntype)     /* output: act. number of characters,
                         or = -1 if element not in d.b. */
/* C end comment */
/* F start doom_gbase

subroutine doom_gbase(name, key_list, nitem, type, ntype)

Returns the basic type name of an object.

Input        (character)          name
Input        (integer array)      key_list
Input        (integer)            nitem: length of type variable
Output       (character)          (basic) type
Output       (integer)            ntype: length of parent, or -1 if no obj.

F end comment */
{
  char rec_key[KEY_LENGTH], lkey[KEY_LENGTH];
  struct object* p;
  int one = 1, n;

#ifdef _TREE_
       puts("++++++ enter doom_gbase");
#endif
  mycpy(lkey, name);
  make_key(lkey, key_list, rec_key);
  p = doom_fetch(rec_key);
  if (p == NULL) *ntype = -1;
  else
    {
     strcpy(type, p->base_name); *ntype = strlen(type);
     n = *ntype + 1;
     doom_blow(type, &n, nitem, &one);
    }
}

/* C start doom_gcorrect: Get corrector settings */
void doom_gcorrect        /* gets corrector values for one element */
        (char* name,     /* input: element name */
         int*  occ,      /* input: occurrence count */
         int*  bits,     /* output: flag bit string */
         double* corr)   /* output: corrector settings */
/* C end comment */
/* F start doom_gcorrect

subroutine doom_gcorrect(name, occur, bits, settings)

Returns the corrector settings in a used sequence

Input        (character)          name
Input        (integer)            occur(rence count)
Output       (integer)            bits (MAD8 flag bit string)
Output       (d.p. array)         settings (two values)

F end comment */
{
  int j, key_list[3] = {2, CORR_CODE, 0};
  char rec_key[KEY_LENGTH], t_key[KEY_LENGTH], lkey[KEY_LENGTH], 
       sequ_name[KEY_LENGTH], 
       type[] = "CORR_SETTING";
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter doom_gcorrect");
#endif
  key_list[2] = *occ;
  mycpy(lkey, name);
  doom_getenv("USE", sequ_name);
  make_key(lkey, key_list, t_key);
  make_dollar_key(sequ_name, t_key, rec_key);
  p = doom_fetch(rec_key);
  if (p == NULL) *bits = -1;
  else if (strcmp(p->obj_type, type) != 0)
    {
     printf("%s%s\n%s%s%s%s\n", " <<< DOOM >>> warning: object ", rec_key,
     " has type: ", p->obj_type, " and not type: ", type);
     *bits = -1;
    }
  else
    {
     *bits = p->a_int[0];
     for (j = 0; j < 2; j++) corr[j] = p->a_dble[j];
    }
}

/* C start doom_gdirect: Get directory */
void doom_gdirect     /* gets DIRECTORY */
     (char* name,    /* input:  name */
      int*  n_name,   /* input/output: max/number of names */
      int*  nitem,   /* input:  no. characters/key */
      char* a_char)  /* output: char array for keys */
/* C end comment */
/* F start doom_gdirect

subroutine doom_gdirect(name, n_name, nitem, keys)

Returns a DOOM directory

Input        (character)          name (main: DIRECTORY)
Input/Output (integer)            n_name (max/actual number of keys)
Input        (integer)            nitem (no. characters/key)
Output       (character array)    keys

F end comment */
{
  char rec_key[KEY_LENGTH], type[] = "DIRECTORY";
  struct object* p;
  int j;

#ifdef _TREE_
       puts("++++++ enter doom_gdirect");
#endif
  mycpy(rec_key, name);
  p = doom_fetch(rec_key);
  if (p == NULL) *n_name = -1;
  else if (strcmp(p->obj_type, type) != 0)
    {
     printf("%s%s\n%s%s%s%s\n", " <<< DOOM >>> warning: object ", rec_key,
     " has type: ", p->obj_type, " and not type: ", type);
     *n_name = -1;
    }
  else
    {
     for (j = 0; j < p->c_char; j++) a_char[j] = p->a_char[j];
     doom_blow(a_char, &p->c_char, nitem, n_name);
    }
}

/* C start doom_gelcode: Get element ISP code (MAD-8) */
int doom_gelcode     /* returns element ISP code (MAD-8), or zero */
     (char* name)     /* input: element name */
/* C end comment */
/* F start doom_gelcode

integer function doom_gelcode(name)

Returns the element ISP code (MAD-8), or zero

Input        (character)          name

F end comment */
{
  char rec_key[KEY_LENGTH];
  int j;
  struct object* p;
#ifdef _TREE_
       puts("++++++ enter doom_gelcode");
#endif
  mycpy(rec_key, name);
  p = doom_fetch("isp_table");
  if (p == NULL)  p = make_isp_table();
  for (j = 0; j < p->c_obj; j++)
    {
     if(strcmp(rec_key, p->names[j]) == 0) return j+1;
    }
  return 0;
}
  
/* C start doom_gelement: Get element parameter values */
void doom_gelement    /* gets parameters for one element */
     (char* name,     /* input: element name */
      int*  npar,     /* input/output: max/act. number of parameters,
                         or = -1 if element not in d.b. */
      double* el_par) /* output: element parameters */
/* C end comment */
/* F start doom_gelement

subroutine doom_gelement(name, npar, elpar)

Returns the element parameter values (see as well doom_gfelem)

Input        (character)          name
Input/Output (integer)     npar - on input maximum number of par. requested
                                - on output number of parameters returned,
                                  or = -1 if element not found
Output       (d.p. array)         elpar - element parameters

F end comment */
{
  char rec_key[KEY_LENGTH], type[] = "ELEMENT";
  struct object* p;
  int j;

#ifdef _TREE_
       puts("++++++ enter doom_gelement");
#endif
  mycpy(rec_key, name);
  p = doom_fetch(rec_key);
  if (p == NULL) *npar = -1;
  else if (strcmp(p->obj_type, type) != 0)
    {
     printf("%s%s\n%s%s%s%s\n", " <<< DOOM >>> warning: object ", rec_key,
     " has type: ", p->obj_type, " and not type: ", type);
     *npar = -1;
    }
  else
    {
     if (*npar > p->c_dble) *npar = p->c_dble;
     for (j = 0; j < *npar ; j++) el_par[j] = p->a_dble[j];
    }
}

/* C start doom_gfelem: Get full element parameters */
void doom_gfelem      /* gets full parameters for one element */
     (char* name,     /* input: element name */
      int*  npar,     /* input/output: max/act. number of parameters,
                         or = -1 if element not in d.b. */
      int*  nchar,   /* input/output: max/act. length of char array */
      int*  isvflg,  /* output: expr. flags */
      double* el_par, /* output: element parameters */
      char* a_char)  /* output: char array for expressions */
/* C end comment */
/* F start doom_gfelem

subroutine doom_gfelem(name, npar, nchar, iexpr, elpar, ex_string)

Returns the (full) element parameters (see as well doom_gelement)

Input        (character)          name
Input/Output (integer)     npar - on input maximum number of par. requested
                                - on output number of parameters returned,
                                  or = -1 if element not found
Input/Output (integer)            nchar (max/actual length of ex_string)
Output       (integer array)      iexpr (expression flags: 1 value, 2 expr.)
Output       (d.p. array)         elpar - element parameters
Output       (character)          ex_string (expression string)

F end comment */
{
  char rec_key[KEY_LENGTH], type[] = "ELEMENT";
  struct object* p;
  int j, nlp;

#ifdef _TREE_
       puts("++++++ enter doom_gfelem");
#endif
  mycpy(rec_key, name);
  p = doom_fetch(rec_key);
  if (p == NULL) *npar = -1;
  else if (strcmp(p->obj_type, type) != 0)
    {
     printf("%s%s\n%s%s%s%s\n", " <<< DOOM >>> warning: object ", rec_key,
     " has type: ", p->obj_type, " and not type: ", type);
     *npar = -1;
    }
  else
    {
     nlp = *npar < p->c_int ? *npar : p->c_int;
     for (j = 0; j < nlp; j++) isvflg[j] = p->a_int[j];
     if (*npar > p->c_dble)  *npar = p->c_dble;
     for (j = 0; j < *npar; j++) el_par[j] = p->a_dble[j];
     if (*nchar > p->c_char)  *nchar = p->c_char;
     for (j = 0; j < *nchar; j++) 
       a_char[j] = p->a_char[j] == '\0' ? '|' : p->a_char[j];
    }
}

/* C start doom_getenv: Get environment variable */
void doom_getenv(  /* returns current value of environment variable */
     char* variable, /* input */
     char* value)    /* output (\0 terminated) */
/* C end comment */
{
  char lkey[KEY_LENGTH];
  int lp, index;

#ifdef _TREE_
       puts("++++++ enter doom_getenv");
#endif
  mycpy(lkey, variable);
  lp = n_list_pos(lkey, o_list_sys[1][0]);
  if (lp < 0)
    {
     index = o_list_sys[1][0]->a_int[-lp-1]; 
     strcpy(value, o_list_sys[1][1]->names[index]);
    }
  else strcpy(value, "        ");
}

/* C start doom_getvar: Get value of internal variable */
void doom_getvar(  /* get DOOM internal variable (ad hoc) */
                 char* name,  /* input */
                 int*  value) /* output */
/* C end comment */
/* F start doom_getvar

subroutine doom_getvar(name, value)

Gets internal variable value

Input        (character)          name  (must match an existing variable)
Output       (int)                value (integer value of name)

F end comment */
{
  char lkey[KEY_LENGTH];
  int j;

#ifdef _TREE_
       puts("++++++ enter doom_getvar");
#endif
  mycpy(lkey, name);
  for (j = 0; j < n_int_var; j++)
    {
     if (strcmp(lkey, l_int_var[j]) == 0) *value = *p_int_var[j];
    }
}

/* C start doom_gfield: Get field errors */
void doom_gfield     /* gets field errors for one element */
     (char* name,    /* input: element name */
      int*  occ,     /* input: occurrence count */
      int*  nerr,    /* input/output: max/returned number of errors,
                             or = -1 if element not found */
      double* f_err) /* output: field  errors */
/* C end comment */
/* F start doom_gfield

subroutine doom_gfield(name, occur, nerr, errors)

Returns the multipole field errors.

Input        (character)          name
Input        (integer)            occurrence count
Input/Output (integer)     nerr - on input maximum number of errors requested
                                - on output number of errors returned,
                                  or = -1 if errors not found
Output       (d.p. array)         (field) errors

F end comment */
{
  int j, key_list[3] = {2, FIELD_CODE, 0};
  char rec_key[KEY_LENGTH], t_key[KEY_LENGTH], lkey[KEY_LENGTH], 
       sequ_name[KEY_LENGTH], 
       type[] = "FIELD_ERROR";
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter doom_gfield");
#endif
  key_list[2] = *occ;
  mycpy(lkey, name);
  doom_getenv("USE", sequ_name);
  make_key(lkey, key_list, t_key);
  make_dollar_key(sequ_name, t_key, rec_key);
  p = doom_fetch(rec_key);
  if (p == NULL) *nerr = -1;
  else if (strcmp(p->obj_type, type) != 0)
    {
     printf("%s%s\n%s%s%s%s\n", " <<< DOOM >>> warning: object ", rec_key,
     " has type: ", p->obj_type, " and not type: ", type);
     *nerr = -1;
    }
  else
    {
     if (*nerr > p->c_dble)  *nerr = p->c_dble;
     for (j = 0; j < *nerr; j++) f_err[j] = p->a_dble[j];
    }
}

/* C start doom_gline:   Get next line from line buffer object */
void doom_gline          /* gets  next line from line buffer object */
        (char* name,     /* input: object name */
         char* line,     /* output: next line if length > 0 */
         int* length)    /* output: last non-blank in line, or 0 for e.o.f 
                                    -1 if object does not exist */
/* C end comment */
/* F start doom_gline

subroutine doom_gline(name, line, length)

Returns the next line from a line buffer object

Input        (character)          (object) name
Output       (character)          (next) line (if length > 0)
Output       (int)                last non-blank in line, or 0 for e.o.f,
                                  -1 if object does not exist

F end comment */
{
  char rec_key[KEY_LENGTH];
  struct object* p;
  int j;

#ifdef _TREE_
       puts("++++++ enter doom_gline");
#endif
  mycpy(rec_key, name);
  p = doom_fetch(rec_key);
  if (p == NULL) *length = -1;
  else
    {
     if (p->a_int[1] == p->a_int[0])  *length = 0;
     else
       {
        strcpy(line, &p->a_char[p->a_int[2]]);
        *length = strlen(line);
        for (j = *length; j < LINE_LENGTH; j++) line[j] = ' ';
        p->a_int[1]++; p->a_int[2] += *length + 1;
       }
    }
}

/* C start doom_gmonitor: Get monitor settings */
void doom_gmonitor        /* gets monitor values for one element */
        (char* name,     /* input: element name */
         int*  occ,      /* input: occurrence count */
         int*  bits,     /* output: flag bit string */
         double* corr)   /* output: monitor settings */
/* C end comment */
/* F start doom_gmonitor

subroutine doom_gmonitor(name, occur, bits, settings)

Returns the monitor settings in a used sequence

Input        (character)          name
Input        (integer)            occur(rence count)
Output       (integer)            bits (MAD8 flag bit string)
Output       (d.p. array)         settings (two values)

F end comment */
{
  int j, key_list[3] = {2, MON_CODE, 0};
  char rec_key[KEY_LENGTH], t_key[KEY_LENGTH], lkey[KEY_LENGTH], 
       sequ_name[KEY_LENGTH], 
       type[] = "MON_SETTING";
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter doom_gmonitor");
#endif
  key_list[2] = *occ;
  mycpy(lkey, name);
  doom_getenv("USE", sequ_name);
  make_key(lkey, key_list, t_key);
  make_dollar_key(sequ_name, t_key, rec_key);
  p = doom_fetch(rec_key);
  if (p == NULL) *bits = -1;
  else if (strcmp(p->obj_type, type) != 0)
    {
     printf("%s%s\n%s%s%s%s\n", " <<< DOOM >>> warning: object ", rec_key,
     " has type: ", p->obj_type, " and not type: ", type);
     *bits = -1;
    }
  else
    {
     *bits = p->a_int[0];
     for (j = 0; j < 2; j++) corr[j] = p->a_dble[j];
    }
}

/* C start doom_gparam: Get parameter */
void doom_gparam(   /* gets parameter from database */
            char* name,    /* input */
            int*  exflag,  /* output: flag: 1 value, > 1 string */
            double* value, /* output: value if exflag = 1 */
            int* l_char,   /* i/o: max/length of string in array a_char */
            char* a_char)  /* output: string if exflag > 1 */
/* C end comment */
/* F start doom_gparam

subroutine doom_gparam(name, exflag, value, iexpr, ex_string)

Returns a parameter.

Input        (character)          name
Output       (integer)            exflag: 1 value, >1 expression
Output       (d.p.)               value
Input/Output (integer)            iexpr - max/actual length of ex_string
Output       (character)          ex_string - expression as read

F end comment */
{
  char rec_key[KEY_LENGTH], type1[] = "CONSTANT", type2[] = "VARIABLE";
  struct object* p;
  int j;

#ifdef _TREE_
       puts("++++++ enter doom_param");
#endif
  mycpy(rec_key, name);
  p = doom_fetch(rec_key);
  if (p == NULL) *l_char = -1;
  else if (strcmp(p->obj_type, type1) != 0 && strcmp(p->obj_type, type2) != 0)
    {
     printf
     (" <<< DOOM >>> warning: object %s has type: %s and not type: %s/%s", 
      rec_key, p->obj_type, type1, type2);
     *l_char = -1;
    }
  else
    {
     *exflag = p->a_int[0]; *value = p->a_dble[0];
     if (*l_char > p->c_char)  *l_char = p->c_char;
     for (j = 0; j < *l_char; j++) a_char[j] = p->a_char[j];  
    }
}

/* C start doom_gparent: Get object parent */
void doom_gparent     /* returns parent name of object */
     (char* name,     /* input: element name */
      int* key_list,  /* input: key_list */
      int* nitem,     /* input: length of parent character variable */
      char* parent,   /* output: parent name */
      int*  npar)     /* output: act. number of characters,
                         or = -1 if element not in d.b. */
/* C end comment */
/* F start doom_gparent

subroutine doom_gparent(name, key_list, nitem, parent, nchar)

Returns the parent name of an object.

Input        (character)          name
Input        (integer array)      key_list
Input        (integer)            nitem: length of parent variable
Output       (character)          parent
Output       (integer)            nchar: length of parent, or -1 if no obj.

F end comment */
{
  char rec_key[KEY_LENGTH], lkey[KEY_LENGTH];
  struct object* p;
  int one = 1, n;

#ifdef _TREE_
       puts("++++++ enter doom_gparent");
#endif
  mycpy(lkey, name);
  make_key(lkey, key_list, rec_key);
  p = doom_fetch(rec_key);
  if (p == NULL) *npar = -1;
  else
    {
     strcpy(parent, p->par_name); *npar = strlen(p->par_name);
     n = *npar + 1;
     doom_blow(parent, &n, nitem, &one);
    }
}

/* C start doom_gpos: Get attribute position (starting at 1) */
int doom_gpos              /* returns attribute position in current list */
         (char* name) /* input: attribute name */
/* C end comment */
/* F start doom_gsequ

integer function doom_gpos(name)

Returns the position of an attribute in the current list (see
doom_setenv), or -1 if not found.

Input        (character)          (attribute) name

F end comment */
{
  int n;
  char rec_key[KEY_LENGTH], lkey[KEY_LENGTH], lname[KEY_LENGTH];
  struct object* p;
  int key_list[2] = {1, TABLE_HEAD};

#ifdef _TREE_
       puts("++++++ enter doom_gpos");
#endif
  doom_getenv("SUB_TABLE", lkey);
  if (*lkey != ' ')
    {
     make_key(lkey, key_list, rec_key);
     p = doom_fetch(rec_key);
     if (p != NULL)
       {
        mycpy(lname, name);
        if ((n = var_list_pos(lname, p)) < 0)  return p->a_int[-n-1]+1;
       }
    }
  return -1;
}

/* C start doom_gsequ: Get expanded sequence */
void doom_gsequ            /* returns sequence table if exists */
         (char* sequ_name, /* input: sequence name */
          int*  max_item,  /* input/output: max./actual no. of elements,
                                 or -1 if not found */
          int*  nitem,     /* no. of characters/name */
          char* names,     /* element names */
          int*  occ,       /* array occ(max_item): occ. counts */
          double* pos)     /* array pos(max_item): centre position */
/* C end comment */
/* F start doom_gsequ

subroutine doom_gsequ(name, nelement, nitem, names, occur, pos)

Returns the list of elements in a sequence. A sequence currently in USE
will be stored automatically at the end of the special MAD.

Input        (character)          (sequence) name
Input/Output (integer) nelement - on input maximum requested
                                - on output actual number in list,
                                  or = -1 if sequence not found
Input        (integer)            nitem = no. of characters/name in names
Output       (character array)    (element) names
Output       (integer array)      occur(rence count of element)
Output       (d.p. array)         centre position

The routine returns a list of elements, and for each element the name, the
occurrence count, the centre position. The arrays must be declared
as follows (minimum sizes):

parameter (ndim = 1 + nelement + 2 * (no. of sublines))
integer occur(ndim)
character *(nitem) names(ndim)
double precision pos(ndim)

F end comment */
{
  char rec_key[KEY_LENGTH], lkey[KEY_LENGTH], type[] = "EX_SEQUENCE";
  struct object* p;
  int j, n;
  int icleng = *nitem * *max_item;
  int key_list[2] = {1, SEQU_CODE};

#ifdef _TREE_
       puts("++++++ enter doom_gsequ");
#endif
  mycpy(lkey, sequ_name);
  make_key(lkey, key_list, rec_key);
  p = doom_fetch(rec_key);
  if (p == NULL) *max_item = -1;
  else if (strcmp(p->obj_type, type) != 0)
    {
     printf("%s%s\n%s%s%s%s\n", " <<< DOOM >>> warning: object ", rec_key,
     " has type: ", p->obj_type, " and not type: ", type);
     *max_item = -1;
    }
  else
    {
     if (icleng > p->c_char)  icleng = p->c_char;
     if (*max_item > p->c_int)  *max_item = p->c_int;
     for (j = 0; j < icleng; j++)  names[j] = p->a_char[j];
     doom_blow(names, &icleng, nitem, &n);
     if (n != *max_item)
     printf("%s%s\n", 
            "doom_gsequ warning - user buffers too small, sequence cut: ",
            rec_key);
     for (j = 0; j < *max_item; j++)
       {
        occ[j]   = p->a_int[j];
        pos[j]   = p->a_dble[j];
       } 
    }
}

/* C start doom_gstring: Get string */
void doom_gstring           /* returns string */
         (char* name,       /* input: key name */
          int* key_list,    /* input */
          int* sleng,       /* input/output: max. length (cut or blank filled)
                                 or -1 if not found */
          char* string)
/* C end comment */
/* F start doom_gstring

subroutine doom_gstring(name, key_list, nchar, string)

Returns a string.

Input        (character)          name
Input        (integer array)      key_list
Input/Output (integer)            nchar - max/actual length of string
Output       (character)          string

F end comment */
{
  char rec_key[KEY_LENGTH], l_key[KEY_LENGTH], type[] = "STRING";
  struct object* p;
  int j, ml = *sleng;

#ifdef _TREE_
       puts("++++++ enter doom_gstring");
#endif
  mycpy(l_key, name);
  make_key(l_key, key_list, rec_key);
  p = doom_fetch(rec_key);
  if (p == NULL) *sleng = -1;
  else if (strcmp(p->obj_type, type) != 0)
    {
     printf("%s%s\n%s%s%s%s\n", " <<< DOOM >>> warning: object ", rec_key,
     " has type: ", p->obj_type, " and not type: ", type);
     *sleng = -1;
    }
  else
    {
     if (ml > p->c_char)  ml = p->c_char;
     for (j = 0; j < ml; j++)
       {
        string[j] = p->a_char[j];
       }
     for (j = ml; j < *sleng; j++) string[j] = ' ';
    } 
}

/* C start doom_gsumm: Get summary record */
void doom_gsumm            /* returns summary table if exists */
         (char* summ_name, /* input: summary table name */
          int*  max_item,  /* input/output: max./actual no. of elements,
                                 or -1 if not found */
          double* table)
/* C end comment */
/* F start doom_gsumm 

subroutine doom_gsumm(name, nchar, table)

Returns a summary table.

Input        (character)          name
Input/Output (integer)            ntable - max/actual length of table
Output       (d.p.)               table (of d.p. values)

F end comment */
{
  char rec_key[KEY_LENGTH], type[] = "SUMMARY";
  struct object* p;
  int j;

#ifdef _TREE_
       puts("++++++ enter doom_gsumm");
#endif
  mycpy(rec_key, summ_name);
  p = doom_fetch(rec_key);
  if (p == NULL) *max_item = -1;
  else if (strcmp(p->obj_type, type) != 0)
    {
     printf("%s%s\n%s%s%s%s\n", " <<< DOOM >>> warning: object ", rec_key,
     " has type: ", p->obj_type, " and not type: ", type);
     *max_item = -1;
    }
  else
    {
     if (*max_item > p->c_dble)  *max_item = p->c_dble;
     for (j = 0; j < *max_item; j++)
       {
        table[j]   = p->a_dble[j];
       } 
    }
}

/* C start doom_gsurv: Get survey table */
void doom_gsurv            /* returns survey table if exists */
         (char* surv_name, /* input: survey name */
          int*  max_item,  /* input/output: max./actual no. of elements,
                                 or -1 if not found */
          int*  nitem,     /* no. of characters/name */
          char* names,     /* element names */
          int*  occ,       /* array occ(max_item): occ. counts */
          double* pos)     /* array pos(max_item,7):
                              s, x, y, z, theta, phi, psi */
/* C end comment */
/* F start doom_gsurv

subroutine doom_gsurv(name, nelem, nitem, elements, occur, pos)

Returns a survey table.

Input        (character)          name (MAD default: SURVEY)
Input/Output (integer)            nelem - max/actual length of table,
                                  or -1 if not found
Output       (character array)    elements (names)
Output       (integer array)      occur(rence counts)
Output       (d.p.array)          table (of d.p. values)

F end comment */
{
  char rec_key[KEY_LENGTH], type[] = "SURVEY";
  struct object* p;
  int j, n;
  int icleng = *nitem * *max_item;

#ifdef _TREE_
       puts("++++++ enter doom_gsurv");
#endif
  mycpy(rec_key, surv_name);
  p = doom_fetch(rec_key);
  if (p == NULL) *max_item = -1;
  else if (strcmp(p->obj_type, type) != 0)
    {
     printf("%s%s\n%s%s%s%s\n", " <<< DOOM >>> warning: object ", rec_key,
     " has type: ", p->obj_type, " and not type: ", type);
     *max_item = -1;
    }
  else
    {
     if (icleng > p->c_char)  icleng = p->c_char;
     if (*max_item > p->c_int)  *max_item = p->c_int;
     for (j = 0; j < icleng; j++)  names[j] = p->a_char[j];
     doom_blow(names, &icleng, nitem, &n);
     if (n != *max_item)
     printf("%s%s\n", 
            "doom_gsurv warning - user buffers too small, survey cut: ",
            rec_key);
     for (j = 0; j < *max_item; j++)
       {
        occ[j]   = p->a_int[j];
       } 
     for (j = 0; j < 7 * *max_item; j++)
       {
        pos[j]   = p->a_dble[j];
       } 
    }
}

/* C start  doom_gtbody: Get table (TWISS, SURVEY etc.) */
void doom_gtbody /* Get table (TWISS, SURVEY etc.) */       
         (char* tabl_name, /* input: table name */
          double* deltap,  /* input: deltap value of twiss table:
                              the table with the nearest deltap is returned */
          char* sequ_name, /* output: sequence name */
          int*  n_column,  /* input/output: max/actual no. of columns */
          int*  n_row,     /* input/output: max/actual no. of rows, or -1 */
          int* name_l,     /* input: no. of characters/name */
          char* names,     /* output: element names */
          int*  occ,       /* output: occ. counts */
          double* optics_t) /* output: table (see twiss_att_names) */
/* C end comment */
/* F start doom_gtbody

subroutine doom_gtbody(tabl_name, deltap, sequ_name, n_column, n_row,
                       name_l, names, occur, optics_t)

Retrieves a table (TWISS, SURVEY etc.)

Input        (character)          tabl_name
Input        (double)             deltap value (the table with nearest deltap
                                                is returned)
Output       (character)          sequ_name
Input/Output (integer)            max/actual n_column
Input/Output (integer)            max/actual n_row, or -1 if not found
Input        (integer)            name_l = no. of characters/name in names
Output       (character array)    (element) names
Output       (integer array)      occur(ence count of element)
Output       (d.p. array)         optics_t

The arrays must be declared as follows:

integer occur(n_row)
character *(nch) names(n_row)
double precision twiss_v(n_column, n_row)

F end comment */
{
  int icleng, ncl, j, jn, n, step, ndelta, posf,
      nrow, suml = MAX_DELTA*OPTICS_SUMM;
  double dn, summ[MAX_DELTA*OPTICS_SUMM];
  char rec_key[KEY_LENGTH], lkey[KEY_LENGTH], type[] = "SUB_TABLE";
  struct object* p;
  int key_list[3] = {2, TABLE_BODY, 0};

#ifdef _TREE_
       puts("++++++ enter doom_gtbody");
#endif
  doom_gthead(tabl_name, name_l, sequ_name, &ndelta, &nrow, &posf, 
  &suml, summ);
  if (ndelta <= 0)
    {
     *n_row = -1; return;
    }
  step = suml / ndelta;
  jn = 0; dn = summ[0];
  for (j = 1; j < ndelta; j++)
    {
     if(fabs(*deltap - summ[step*j]) < fabs(*deltap - dn))
       {
        jn = j; dn = summ[step*j];
       }
    }
  key_list[2] = jn;
  mycpy(lkey, tabl_name);
  make_key(lkey, key_list, rec_key);
  p = doom_fetch(rec_key);
  if (p == NULL) *n_row = -1;
  else if (strcmp(p->obj_type, type) != 0)
    {
     printf("%s%s\n%s%s%s%s\n", " <<< DOOM >>> warning: object ", rec_key,
     " has type: ", p->obj_type, " and not type: ", type);
     *n_row = -1;
    }
  else
    {
     ncl = p->c_dble / p->c_int;
     if (*n_row > p->c_int)  *n_row = p->c_int;
     if ((icleng = *name_l * (*n_row)) < p->c_char)  icleng = p->c_char;
     jn = icleng > p->c_char ? p->c_char : icleng;
     for (j = 0; j < jn; j++)  names[j] = p->a_char[j];
     doom_blow(names, &jn, name_l, &n);
     if (n != *n_row)
     printf("%s%s\n", 
            "doom_gtbody warning - user buffers too small, Twiss table cut: ",
            rec_key);
     for (j = 0; j < *n_row; j++)
       {
        occ[j]   = p->a_int[j];
       } 
     n = 0;
     for (jn = 0; jn < *n_row; jn++)
       {
        for (j = 0; j < *n_column; j++)
	    optics_t[jn*(*n_column)+j] = p->a_dble[n+j];
        n+= ncl;
       } 
     if (*n_column > ncl)  *n_column = ncl;
     delete_obj(p);
    }
}

/* C start doom_gthead: Get table master record */
void doom_gthead
            (char* tabl_name, /* input: table name */
             int* name_l,     /* input: character length of sequ_name */
             char* sequ_name, /* output: sequence name */
             int* n_delta,    /* output: # of delta values*/
             int* n_row,      /* output: Twiss table length (per delta) */
             int* pos_flag,   /* output: element position flag:
                                 1 start, 2 centre, 3 end */
             int* sum_leng,   /* input/output: max/actual total length  
                                               of summary record */
             double* sum_rec) /* output: summary record */
/* C end comment */
/* F start doom_gthead

subroutine doom_gthead(tabl_name, name_l, n_deltap, deltap,
                        n_row, pos_flag, sum_leng, sum_rec)
Returns the table master record.

Input        (character)         tabl_name
Input        (integer)           name_l = no. of characters/name in names
Output       (character)         sequ_name
Output       (integer)           actual n_delta (# of delta values) or -1
Output       (integer)           actual n_row = table length per delta or -1
Output       (integer)           pos_flag = element position flag:
                                 1 start, 2 centre, 3 end
Input/Output (integer)           max/actual sum_leng = total summary length
Output       (d.p. array)        sum_rec = summary record

F end comment */
{
  char rec_key[KEY_LENGTH], lcp[KEY_LENGTH], type[] = "TABLE_HEADER";
  int j, k, nc, key_list[2] = {1, TABLE_HEAD};
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter doom_gthead");
#endif
  mycpy(lcp, tabl_name);
  make_key(lcp, key_list, rec_key);
  p = doom_fetch(rec_key);
  if (p == NULL)
    {
     *n_delta = -1; *n_row = -1;
    }
  else if (strcmp(p->obj_type, type) != 0)
    {
     printf("%s%s\n%s%s%s%s\n", " <<< DOOM >>> warning: object ", rec_key,
     " has type: ", p->obj_type, " and not type: ", type);
     *n_row = -1;
     *n_delta = -1;
    }
  else
    {
     nc = p->c_int - OPTICS_INTL;
     *n_delta  = p->a_int[nc];
     *n_row = p->a_int[nc+1];
     k = p->a_int[nc] * p->a_int[nc+2];
     *pos_flag = p->a_int[nc+3];
     if (*sum_leng > k) *sum_leng = k;
     for (j = 0; j < *sum_leng; j++)
       {
        sum_rec[j] = p->a_dble[j]; 
       }
     strcpy(sequ_name, p->par_name);
     j = strlen(p->par_name) + 1;
     doom_blow(sequ_name, &j, name_l, &k);
     doom_setenv("SUB_TABLE", lcp);
    }
}

/* C start doom_gtime: Get object time of creation or mod. */
void  doom_gtime            /* returns -1. if not found, else mod. time */
            (char* name,    /* input: (full) d.b. object name */
             int* key_list, /* input: key list for composed name */
             double* time)  /* output */
/* C end comment */
/* F start doom_gtime

subroutine doom_gtime(name, key_list, time)

Returns the time of last modification for an object.

Input        (character)          name
Input        (integer)            key_list
Output       (d.p.)               time (sec.s since 1.1.1997 + microsteps)

F end comment */
{
  char rec_key[KEY_LENGTH], lkey[KEY_LENGTH];
  int lp, index;

#ifdef _TREE_
       puts("++++++ enter doom_gtime");
#endif
  mycpy(lkey, name);
  make_key(lkey, key_list, rec_key);
  lp = n_list_pos(rec_key, o_list_sys[0][0]);
  if (lp < 0)  /* object exists */
    {
     index = o_list_sys[0][0]->a_int[-lp-1];
     if (o_list_sys[0][0]->p_obj[index] == NULL)
       *time = (double) o_list_sys[0][1]->a_int[index] +
               (double) o_list_sys[0][2]->a_int[index] * micro_step;
     else
       *time = (double) o_list_sys[0][0]->p_obj[index]->ma_time +
               (double) o_list_sys[0][0]->p_obj[index]->mi_time * micro_step;
    }
  else *time = (double) -1;
}

/* C start doom_gtime_env: Get env$object time */
void doom_gtime_env         /* returns -1 if not found, else mod. time */
            (char* env,     /* input: precede with environment value+$ */
             char* name,    /* input: */
             int* key_list, /* input: key list for composed name */
             double* time)  /* output */
/* C end comment */
/* F start doom_gtime_env

subroutine doom_gtime_env(env, name, key_list, time)

Returns the time of last modification for an object with a composed name.

Input        (character)          env    (environment variable)
Input        (character)          name
                        --------> the object key is then env_name...
Input        (integer)            key_list
Output       (d.p.)               time (sec.s since 1.1.1997 + microsteps)

F end comment */
{
  char rec_key[KEY_LENGTH], t_key[KEY_LENGTH], lkey[KEY_LENGTH],
       env_name[KEY_LENGTH], env_value[KEY_LENGTH];
  int lp, index;

#ifdef _TREE_
       puts("++++++ enter doom_gtime_env");
#endif
  mycpy(lkey, name);
  if (*env != ' ' && *env != '\0')
    {
     mycpy(env_name, env);
     doom_getenv(env_name, env_value);
     make_key(lkey, key_list, t_key);
     make_dollar_key(env_value, t_key, rec_key);
    }
  else make_key(lkey, key_list, rec_key);
  lp = n_list_pos(rec_key, o_list_sys[0][0]);
  if (lp < 0)  /* object exists */
    {
     index = o_list_sys[0][0]->a_int[-lp-1];
     if (o_list_sys[0][0]->p_obj[index] == NULL)
       *time = (double) o_list_sys[0][1]->a_int[index] +
               (double) o_list_sys[0][2]->a_int[index] * micro_step;
     else
       *time = (double) o_list_sys[0][0]->p_obj[index]->ma_time +
               (double) o_list_sys[0][0]->p_obj[index]->mi_time * micro_step;
    }
  else *time = (double) -1;
}

/* C start doom_gtnames: Get optics table column names */
void doom_gtnames
            (char* tabl_name, /* input: table name */
             int* name_l,     /* input: character length of sequ_name */
             int* n_names,    /* input/output: max/actual # names, or -1 */
             char* names)     /* output: column names in order of occ. */
/* C end comment */
/* F start doom_gtnames

subroutine doom_gtnames(tabl_name, name_l, n_names, names)
Returns the optics table column names

Input        (character)         tabl_name
Input        (integer)           name_l = no. of characters/name in names
Input/Output (integer)           n_names = max/actual no. of names, or -1
Output       (character array)   names

F end comment */
{
  char rec_key[KEY_LENGTH], lcp[KEY_LENGTH], type[] = "TABLE_HEADER";
  int i, j, k, nnt, key_list[2] = {1, TABLE_HEAD};
  char *c1, *c2;
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter doom_gtnames");
#endif
  mycpy(lcp, tabl_name);
  make_key(lcp, key_list, rec_key);
  p = doom_fetch(rec_key);
  if (p == NULL) *n_names = -1;
  else if (strcmp(p->obj_type, type) != 0)
    {
     printf("%s%s\n%s%s%s%s\n", " <<< DOOM >>> warning: object ", rec_key,
     " has type: ", p->obj_type, " and not type: ", type);
     *n_names = -1;
    }
  else
    {
     nnt = p->c_int - OPTICS_INTL;
     if (*n_names > nnt)  *n_names = nnt;
     for (j = 0; j < p->c_char; j++)  cdummy[j] = p->a_char[j];
     doom_blow(cdummy, &j, name_l, &k);
     for (k = 0; k < nnt; k++)
       {
        if ((i = p->a_int[k]) < *n_names)
	  {
           c1 = names+i*(*name_l); c2 = cdummy+k*(*name_l);
           for (j = 0; j < *name_l; j++) *c1++ = *c2++;
	  }
       }
    }
}

/* C start doom_gtrack: Get track master record */
void doom_gtrack    /* gets track survival table */
     (int*  locc,   /* input/output: max/number of list items,
                                     or = -1 if not found */
      int*  list)   /* output: list of turn numbers, and surviving particles */
/* C end comment */
/* F start doom_gtrack

subroutine doom_gtrack(length, list)

Returns the table containing the list of turns with results saved.

Input/Output (integer)            length = maximum/actual no. turns in list
Output       (integer array)      (2,*): list (of turn numbers,
                                  and # of surving particle)

F end comment */
{
  char rec_key[KEY_LENGTH], name[] = "TRACK_SURVIVAL", type[] = "SURVIVAL";
  int j, key_list[3] = {1, TRACK_SV};
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter doom_gtrack");
#endif
  make_key(name, key_list, rec_key);
  p = doom_fetch(rec_key);
  if (p == NULL) *locc = -1;
  else if (strcmp(p->obj_type, type) != 0)
    {
     printf("%s%s\n%s%s%s%s\n", " <<< DOOM >>> warning: object ", rec_key,
     " has type: ", p->obj_type, " and not type: ", type);
     *locc = -1;
    }
  else
    {
     if (*locc > p->c_int/2)  *locc = p->c_int/2;
     for (j = 0; j < p->c_int; j++)  list[j] = p->a_int[j];
    }
}

/* C start doom_gturn: Get one turn record */
void doom_gturn(     /* return the track data of one turn */
      int* turn,     /* input: turn number */
      int* number,   /* input/output: max/number of particles,
                             or -1 if not found */
      int* list,     /* output: particle list */
      double* track) /* output: track coordinates (6,...) */
/* C end comment */
/* F start doom_gturn

subroutine doom_gturn(turn, npart, list, coords)

Returns the particle numbers and their coordinates for a specific turn.

Input        (integer)            turn (number)
Input/Output (integer)      npart - input: max. number of particles requested
                                  - output: number of particles returned
Output       (integer array)      list (of particle numbers)
Output       (d.p. array)         coords = particle coordinates

Declarations required:

integer list(..)
double precision coords(6,..)

F end comment */
{
  struct object* p;
  char rec_key[KEY_LENGTH], lkey[KEY_LENGTH];
  int j, key_list[3] = {1, TRACK_RECORD};
  char turn_num[16], name[] = "TURN", type[] = "TURN";

#ifdef _TREE_
       puts("++++++ enter doom_gturn");
#endif
  sprintf(turn_num, "%d", *turn);
  make_dollar_key(name, turn_num, lkey);
  make_key(lkey, key_list, rec_key);
  p = doom_fetch(rec_key);
  if (p == NULL) *number = -1;
  else if (strcmp(p->obj_type, type) != 0)
    {
     printf("%s%s\n%s%s%s%s\n", " <<< DOOM >>> warning: object ", rec_key,
     " has type: ", p->obj_type, " and not type: ", type);
     *number = -1;
    }
  else
    {
     if (*number > p->c_int)  *number = p->c_int;
     for (j = 0; j < *number; j++)  list[j] = p->a_int[j];
     for (j = 0; j < 6 * *number; j++)  track[j] = p->a_dble[j];
     delete_obj(p);
    }
}

/* C start doom_gtype: Get object type */
void doom_gtype       /* returns type of object */
     (char* name,     /* input: element name */
      int* key_list,  /* input: key_list */
      int* nitem,     /* input: length of parent character variable */
      char* type,     /* output: type */
      int* ntype)     /* output: act. number of characters,
                         or = -1 if element not in d.b. */
/* C end comment */
/* F start doom_gtype

subroutine doom_gtype(name, key_list, nitem, type, ntype)

Returns the type of an object.

Input        (character)          name
Input        (integer array)      key_list
Input        (integer)            nitem: length of type variable
Output       (character)          type
Output       (integer)            ntype: length of parent, or -1 if no obj.

F end comment */
{
  char rec_key[KEY_LENGTH], lkey[KEY_LENGTH];
  struct object* p;
  int one = 1, n;

#ifdef _TREE_
       puts("++++++ enter doom_gtype");
#endif
  mycpy(lkey, name);
  make_key(lkey, key_list, rec_key);
  p = doom_fetch(rec_key);
  if (p == NULL) *ntype = -1;
  else
    {
     strcpy(type, p->obj_type); *ntype = strlen(p->obj_type);
     n = *ntype + 1;
     doom_blow(type, &n, nitem, &one);
    }
}

void doom_list()
{
  if (db_is_open() == 0)     puts("+++DOOM - db not open");
  else                       db_list();
}

int doom_lpos(int* number, int* c_int, int* list)
{
  int mid, low = 0, last = 0, high = *c_int - 1;
  while (low <= high)
    {
     mid = (low + high) / 2;
     if (*number < list[mid])  {high = mid - 1; last = mid;}
     else if (*number > list[mid]) {low  = mid + 1; last = low;}
     else              return (-mid-1);
    }
    return last;
}

void doom_memory()
{
  int amount = 250000;
  int* ip;

#ifdef _TREE_
       puts("++++++ enter doom_memory");
#endif
   ip = (int*) malloc(amount * sizeof(int));
   printf("<<< DOOM >>> current upper memory: %d\n", (int) ip);
   free(ip);
}

int doom_nadd( /* Adds an item to a string list (fixed item size) */
    char* carg, /* item to be added */
    char* s_list, /* input/output; list */
    int* c_list,  /* input/output: no. of items in s_list (augmented) */
    int* item)    /* input: item length */
{
  int lp, j, k = *c_list;

#ifdef _TREE_
       puts("++++++ enter doom_nadd");
#endif
  lp = doom_npos(carg, s_list, c_list, item);
  if (lp > -1)
    {
     for (j = *c_list; j > lp; j--)  
          strncpy(&s_list[*item * j], &s_list[*item * (j-1)], *item);
     strncpy(&s_list[*item * lp], carg, *item);
     k++;
    }
  return k;
}

int doom_newdb()
{
  return newdb_flag;
}

int doom_npos(  /* Returns list position if < 0 or list place if >= 0 */
    char* carg,   /* input: string to be looked up */
    char* s_list, /* input: list of character items */
    int*  c_list, /* input: no. of items in s_list */
    int* item)    /* input: item length (fixed) */
{
  int comp, mid, low = 0, last = 0, high = *c_list - 1;

#ifdef _TREE_
       puts("++++++ enter doom_npos");
#endif
  while (low <= high)
    {
     mid = (low + high) / 2;
     comp = strncmp(carg, &s_list[*item * mid], *item);
     if (comp < 0)      {high = mid - 1; last = mid;}
     else if (comp > 0) {low  = mid + 1; last = low;}
     else               return (-mid-1);
    }
    return last;
}

/* C start doom_open: Open database */
void doom_open(char* filename)
/* C end comment */
/* F start doom_open

subroutine doom_open(file_name)

Opens the database if existing, else creates and opens a new one.

Input        (character)          file_name (of the database)

F end comment */
{

  FILE *fpnt;
  int j;

#ifdef _TREE_
       puts("++++++ enter doom_open");
#endif
  if (cdummy == NULL) cdummy = (char*) malloc(C_DUM_SIZE);
  if (char_tok == NULL)   char_tok   = (char*) malloc(C_DUM_SIZE);
  if (buffer == NULL)  buffer = (char*) malloc(buff_size);
  myscpy(cdummy, filename);
  printf("\n%s%d%s%d\n%s%s\n\n", 
         " <<< DOOM >>> program version ", major_version, ".",
         minor_version,
         "          >>> opening data base file: ", cdummy);
  fpnt = fopen(cdummy, "r");  /* check for existence, set flag */
  if (fpnt == NULL)
    {
     if (force_old != 0)
       {
        printf
        (" <<< DOOM >>> fatal: no pre-existing database\n"); 
        exit(1);
       }
     newdb_flag = 1;
    }
  else
    {
     newdb_flag = 0; fclose(fpnt);
    }
  if (db_open(cdummy, newdb_flag) == 0)
    {
     printf("%s\n", "          >>> error when opening db");
     exit(1);
    }
  else
    {
     if (o_list_sys == NULL)  o_list_sys 
     = (struct object***) calloc(o_list_numb, sizeof(struct object**));
     for (j = 0; j < o_list_numb; j++)
       {
        if (o_list_sys[j] == NULL)  o_list_sys[j] 
        = (struct object**) calloc(o_list_count, sizeof(struct object*));
       }
     if (newdb_flag == 0)
       {
        open_old();
        if (reopen > 0) 
            printf("objects freed, kept from d.b. in memory: %d %d\n",
                   prev_freed, prev_kept);
       }
     else open_new(); /* create d.b. from scratch */
    }
  make_volatile();
}

/* C start doom_palign: Put alignment errors */
void doom_palign     /* puts alignment errors for one element */
     (char* name,    /* input: element name */
      int*  occ,     /* input: occurrence count */
      int*  nerr,    /* input: number of errors */
     double* al_err) /* input: alignment errors */
/* C end comment */
/* F start doom_palign

subroutine doom_palign(name, occur, nerr, errors)

Stores the alignment errors.

Input        (character)          name
Input        (integer)            occur(ence count)
Input        (integer)            nerr = number of errors present
Input        (d.p. array)         (alignment) errors

F end comment */
{
  int j, key_list[3] = {2, ALIGN_CODE, 0};
  char rec_key[KEY_LENGTH], t_key[KEY_LENGTH], lkey[KEY_LENGTH], 
       sequ_name[KEY_LENGTH], 
       type[] = "ALIGN_ERROR";
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter doom_palign");
#endif
  key_list[2] = *occ;
  mycpy(lkey, name);
  doom_getenv("USE", sequ_name);
  make_key(lkey, key_list, t_key);
  make_dollar_key(sequ_name, t_key, rec_key);
  p = doom_fetch(rec_key);
  if (p != NULL)
    {
     if (*nerr > p->l_dble) 
       { 
        grow_obj(p, 0, *nerr, 0, 0);
        p->c_dble = *nerr;
       }
     for (j = 0; j < *nerr; j++) p->a_dble[j] = al_err[j];
    }
  else
    {
     p = make_obj(rec_key, 0, *nerr, 0, 0);
     strcpy(p->par_name, sequ_name);
     strcpy(p->obj_type, type);
     fill_obj(p, 0, *nerr, 0, idummy, al_err, cdummy);
    }
  doom_save(p);
}

/* C start doom_pcorrect: Put corrector settings */
void doom_pcorrect    /* puts corrector setting for one element */
     (char* name,    /* input: element name */
      int*  occ,     /* input: occurrence count */
      int*  bits,    /* input: bit pattern */
      double* corr)  /* input: corrector settings */
/* C end comment */
/* F start doom_pcorrect

subroutine doom_pcorrect(name, occur, bits, settings)

Stores the corrector settings in a used sequence

Input        (character)          name
Input        (integer)            occur(rence count)
Input        (integer)            bits (MAD8 flag bit string)
Input        (d.p. array)         settings (two values)

F end comment */
{
  int j, key_list[3] = {2, CORR_CODE, 0};
  char rec_key[KEY_LENGTH], t_key[KEY_LENGTH], lkey[KEY_LENGTH], 
       sequ_name[KEY_LENGTH],
       type[] = "CORR_SETTING";
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter doom_pcorrect");
#endif
  key_list[2] = *occ;
  mycpy(lkey, name);
  doom_getenv("USE", sequ_name);
  make_key(lkey, key_list, t_key);
  make_dollar_key(sequ_name, t_key, rec_key);
  p = doom_fetch(rec_key);
  if (p != NULL)
    {
     p->a_int[0] = *bits;
     for (j = 0; j < 2; j++)  p->a_dble[j] = corr[j];
    }
  else
    {
     p = make_obj(rec_key, 1, 2, 0, 0);
     strcpy(p->par_name, sequ_name);
     strcpy(p->obj_type, type);
     fill_obj(p, 1, 2, 0, bits, corr, cdummy);
    }
  doom_save(p);
}

/* C start doom_pdirect: Put directory */
void doom_pdirect     /* puts directory object */
     (char* name,    /* input:  name */
      int*  n_name,   /* input: number of names */
      int*  nitem,   /* input: (fixed) length of names */
      char* a_char)  /* input: char array for names */
/* C end comment */
/* F start doom_pdirect

subroutine doom_pdirect(name, n_name, nitem, keys)

Stores a DOOM directory

Input        (character)          name (main: DIRECTORY)
Input        (integer)            n_name (number of keys)
Input        (integer)            nitem (no. characters/key)
Output       (character array)    keys

F end comment */
{
  char rec_key[KEY_LENGTH], type[] = "DIRECTORY";
  struct object* p;
  int c_char = *nitem * *n_name;

#ifdef _TREE_
       puts("++++++ enter doom_pdirect");
#endif
  if (c_char > 0)
    {
     mycpy(rec_key, name);
     p = doom_fetch(rec_key);
     if (p != NULL)
       {
        if (c_char > p->l_char || *n_name > p->l_obj) 
           { 
            grow_obj(p, 0, 0, c_char, *n_name);
           }
       }
     else
       {
        p = make_obj(rec_key, 0, 0, c_char, *n_name);
        strcpy(p->par_name, type);
        strcpy(p->obj_type, type);
       }
     strncpy(p->a_char, a_char, c_char);
     doom_comp_max(p->a_char, &c_char, nitem, &p->c_char);
     p->c_obj = p->l_obj;
     doom_save(p);
    }
}

/* C start doom_pelement: Put element parameter values */
void doom_pelement    /* puts parameter values in existing element */
     (char* name,     /* input: element name */
      int*  npar,     /* input: number of parameters */
      double* el_par) /* input: element parameters */
/* C end comment */
/* F start doom_pelement

subroutine doom_pelement(name, npar, elpar)

Stores element parameter values (see as well doom_pfelem)

Input        (character)          name
Input        (integer)            npar - number of parameters in elpar
Input        (d.p. array)         elpar - element parameters

Only the first npar parameters are replaced, the rest remains unchanged.

F end comment */
{
  char rec_key[KEY_LENGTH];
  struct object* p;
  int j, nlp;

#ifdef _TREE_
       puts("++++++ enter doom_pelement");
#endif
  mycpy(rec_key, name);
  p = doom_fetch(rec_key);
  if (p != NULL)
    {
     nlp = p->c_dble < *npar ? p->c_dble : *npar;
     for (j = 0; j < nlp ; j++)  p->a_dble[j] = el_par[j];
     for (j = 0; j < p->c_int; j++)  
       {
        if (p->a_int[j] > 1) p->a_int[j] = 1; /* expr. flag */
       }
     doom_save(p);
    }
}

/* C start doom_pfelem: Put full element parameters */
void doom_pfelem     /* puts full parameters for one element */
     (char* name,    /* input: element name */
      int*  isp,     /* input: MAD element type (1=drift etc.) */
      char* parent,  /* input: parent name */
      int*  npar,    /* input: number of parameters */
      int*  nchar,   /* input: length of char array */
      int*  isvflg,  /* input: expr. flags */
     double* el_par, /* input: element parameters */
      char* a_char)  /* input: char array for expressions */
/* C end comment */
/* F start doom_pfelem

subroutine doom_pfelem(name, isp, parent, npar, nchar, iexpr, elpar, ex_string)

Stores the (full) element parameters (see as well doom_pelement)

Input        (character)          name
Input        (integer)            isp - MAD element type (1=drift etc.)
Input        (character)          parent
Input        (integer)            npar - number of parameters
Input        (integer)            nchar (length of ex_string)
Input        (integer array)      iexpr (expression flags: 1 value, 2 expr.)
Output       (d.p. array)         elpar - element parameters
Output       (character)          ex_string (expression string,
                                  seperator: '|')

F end comment */
{
  char rec_key[KEY_LENGTH], type[] = "ELEMENT";
  struct object* p;
  struct object* pt;

#ifdef _TREE_
       puts("++++++ enter doom_pfelem");
#endif
  mycpy(rec_key, name);
  p = doom_fetch(rec_key);
  if (p != NULL)  
    {
     if (*nchar > p->l_char)  grow_obj(p, 0, 0, *nchar, 0);
     if (*npar  > p->l_int)   grow_obj(p, *npar, *npar, 0, 0);
    }
  else p = make_obj(rec_key, *npar, *npar, *nchar, *npar);
  mycpy(p->par_name, parent);
  strcpy(p->obj_type, type);
  pt = doom_fetch("isp_table");
  if (pt == NULL)  pt = make_isp_table();
  if (pt == NULL || pt->c_obj == 0)
     {
      printf
      (" <<< DOOM >>> fatal: NULL or empty MAD type table\n"); 
      exit(1);
     }
  if (*isp < 1 || *isp > pt->c_obj)
     {
      printf
      (" <<< DOOM >>> fatal: (pfelem) illegal MAD type %d\n", 
       *isp);
      exit(1);
     }
  strcpy(p->base_name, pt->names[*isp-1]);
  fill_obj(p, *npar, *npar, *nchar, isvflg, el_par, a_char);
  p->c_obj = p->l_obj;
  doom_save(p);
}

/* C start doom_pfield: Put field errors */
void doom_pfield     /* puts field errors for one element */
     (char* name,    /* input: element name */
      int*  occ,     /* input: occurrence count */
      int*  nerr,    /* input: number of errors */
      double* f_err) /* input: field errors */
/* C end comment */
/* F start doom_pfield

subroutine doom_pfield(name, occur, nerr, errors)

Stores the multipole field errors.

Input        (character)          name
Input        (character)          name
Input        (integer)            occur(ence count)
Input        (integer)            nerr = number of errors present
Input        (d.p. array)         (field) errors

F end comment */
{
  int j, key_list[3] = {2, FIELD_CODE, 0};
  char rec_key[KEY_LENGTH], t_key[KEY_LENGTH], lkey[KEY_LENGTH], 
       sequ_name[KEY_LENGTH],
       type[] = "FIELD_ERROR";
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter doom_pfield");
#endif
  key_list[2] = *occ;
  mycpy(lkey, name);
  doom_getenv("USE", sequ_name);
  make_key(lkey, key_list, t_key);
  make_dollar_key(sequ_name, t_key, rec_key);
  p = doom_fetch(rec_key);
  if (p != NULL)
    {
     if (*nerr > p->l_dble) 
       { 
        grow_obj(p, 0, *nerr, 0, 0);
        p->c_dble = *nerr;
       }
     for (j = 0; j < *nerr; j++)  p->a_dble[j] = f_err[j];
    }
  else
    {
     p = make_obj(rec_key, 0, *nerr, 0, 0);
     strcpy(p->par_name, sequ_name);
     strcpy(p->obj_type, type);
     fill_obj(p, 0, *nerr, 0, idummy, f_err, cdummy);
    }
  doom_save(p);
}

/* C start doom_pline:   Put next line into line buffer object */
void doom_pline          /* puts  next line into line buffer object */
        (char* name,     /* input: object name */
         char* line,     /* input: next line if length */
         int* length)    /* last non-blank in line */
/* C end comment */
/* F start doom_pline

subroutine doom_pline(name, line, length)

Puts the next line into a line buffer object

Input        (character)          (object) name
Input        (character)          (next) line
Input        (int)                length = last non-blank in line
F end comment */
{
  char rec_key[KEY_LENGTH];
  struct object* p;
  int lp, index;

#ifdef _TREE_
       puts("++++++ enter doom_pline");
#endif
  if (*length > 0)
    {
     mycpy(rec_key, name);
     lp = n_list_pos(rec_key, o_list_sys[0][0]);
     if (lp < 0)  /* object exists */
       {
        index = o_list_sys[0][0]->a_int[-lp-1];
        if ((p = o_list_sys[0][0]->p_obj[index]) == NULL)
          {
           printf
           (" <<< DOOM >>> fatal: object %s in list, not in d.b.\n", 
            rec_key);
           exit(1);
          }
       }
     else
       {
        p = make_obj(rec_key, 3, 0, 100 * LINE_LENGTH, 0);
        p->c_int = p->l_int;
        doom_save(p);
       }
     if (p->c_char+LINE_LENGTH+1 > p->l_char) 
         grow_obj(p, 0, 0, 2*p->l_char, 0);
     strncpy(&p->a_char[p->c_char], line, *length);
     p->c_char += *length;
     p->a_char[p->c_char++] = '\0';
     p->a_int[0]++;
    }
}

void doom_plink(        /* saves a link area */
             char* name,
             int* link_start,
             int* link_end)
{
  int k, zero = 0, c_n = 0, clint = link_end - link_start + 1;

#ifdef _TREE_
       puts("++++++ enter doom_plink");
#endif
  for (k = 0; k < clint; k++)
    {
     if (link_start[k] != 0)
       {
        c_n++; 
       } 
    }
  doom_store(name, &zero, &clint, &zero, &zero, link_start, ddummy, cdummy);
  if (debug_flag > 0)
    {
     printf("%s%s%s%d\n", "DOOM - link area saved:  ", name,
            "  count: ", c_n);
    }
}

/* C start doom_pmonitor: Put monitor settings */
void doom_pmonitor    /* puts monitor setting for one element */
     (char* name,    /* input: element name */
      int*  occ,     /* input: occurrence count */
      int*  bits,    /* input: bit pattern */
      double* corr)  /* input: monitor settings */
/* C end comment */
/* F start doom_pmonitor

subroutine doom_pmonitor(name, occur, bits, settings)

Stores the monitor settings in a used sequence

Input        (character)          name
Input        (integer)            occur(rence count)
Input        (integer)            bits (MAD8 flag bit string)
Input        (d.p. array)         settings (two values)

F end comment */
{
  int j, key_list[3] = {2, MON_CODE, 0};
  char rec_key[KEY_LENGTH], t_key[KEY_LENGTH], lkey[KEY_LENGTH], 
       sequ_name[KEY_LENGTH],
       type[] = "MON_SETTING";
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter doom_pmonitor");
#endif
  key_list[2] = *occ;
  mycpy(lkey, name);
  doom_getenv("USE", sequ_name);
  make_key(lkey, key_list, t_key);
  make_dollar_key(sequ_name, t_key, rec_key);
  p = doom_fetch(rec_key);
  if (p != NULL)
    {
     p->a_int[0] = *bits;
     for (j = 0; j < 2; j++)  p->a_dble[j] = corr[j];
    }
  else
    {
     p = make_obj(rec_key, 1, 2, 0, 0);
     strcpy(p->par_name, sequ_name);
     strcpy(p->obj_type, type);
     fill_obj(p, 1, 2, 0, bits, corr, cdummy);
    }
  doom_save(p);
}

/* C start doom_pparam: Put parameter */
void doom_pparam(   /* puts parameter or constant in database */
            char* name,    /* input */
            int*  flag,    /* input: flag mod 10: 1 value, > 1 string
                                     flag/10: 1 = constant, 2 = variable */
            double* value, /* input: value if exflag = 1 */
            int* l_char,   /* input: length of string in array a_char */
            char* a_char)  /* input: string if exflag > 1 */
/* C end comment */
/* F start doom_pparam

subroutine doom_pparam(name, type, exflag, value, lexpr, expr)

Stores a parameter.

Input        (character)          name
Input        (integer)            exflag: 1 value, >1 expression
Input        (d.p.)               value (if exflag = 1)
Input        (integer)            lexpr - length of expr(ession) if exflag > 1
Input        (character)          expr - expression as read

F end comment */
{
  char rec_key[KEY_LENGTH]; 
  struct object* p;
  int vexpr, lv = 1, vtype = *flag/10, lmch = *l_char == 0 ? 1 : *l_char;

  vexpr = *flag - 10 * vtype;
#ifdef _TREE_
       puts("++++++ enter doom_pparam");
#endif
  mycpy(rec_key, name);
  p = doom_fetch(rec_key);
  if (p != NULL)  
    {
     if (lmch > p->l_char)  grow_obj(p, 0, 0, lmch, 0);
    }
  else p = make_obj(rec_key, lv, lv, lmch, lv);
  p->a_char[0] = '\0';
  if (vtype == 1)
    {
     strcpy(p->par_name, "REAL_CONSTANT");
     strcpy(p->obj_type, "CONSTANT");
    }
  else
    {
     strcpy(p->par_name, "REAL_VARIABLE");
     strcpy(p->obj_type, "VARIABLE");
    }
  fill_obj(p, lv, lv, *l_char, &vexpr, value, a_char);
  p->c_char = lmch;
  p->c_obj = 1;
  doom_save(p);
}

/* C start doom_pparent: Put object parent */
void doom_pparent     /* puts parent name into object */
     (char* name,     /* input: object name */
      int* key_list,  /* input: key_list */
      char* parent)   /* input: parent name */
/* C end comment */
/* F start doom_pparent

subroutine doom_pparent(name, key_list, parent)

Stores the parent name of an object.

Input        (character)          name
Input        (integer array)      key_list
Input       ( character)          parent

F end comment */
{
  char rec_key[KEY_LENGTH], lkey[KEY_LENGTH];
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter doom_pparent");
#endif
  mycpy(lkey, name);
  make_key(lkey, key_list, rec_key);
  p = doom_fetch(rec_key);
  if (p != NULL) mycpy(p->par_name, parent);
}

/* C start  doom_prtab: Prepare name table from Fortran */
void doom_prtab
         (char* name_list,  /* input: list of names */
          int*  n_names,    /* input: number of names */
          int*  name_l,     /* input: number of characters/name */
          char* prep_list,  /* output: list in DOOM format (ordered, \0 sep */
          int* positions)   /* output: positions in value table
                               The initial order is assumed to be 0, 1, 2, */
/* C end comment */
/* F start doom_prtab

subroutine doom_prtab(name_list, n_names, nch, prep_list, positions)

Prepares a Fortran name table for usage in DOOM 

Input        (character)          name_list
Input        (integer)            n_names
Input        (integer)            nch = no. of characters/name in names
Output       (character array)    prep_list
Output       (integer array)      positions

F end comment */
{
  int i, j, k = 1, temp, out_length, in_length = *n_names * (*name_l);
  struct object* p;
  char *c = prep_list;
  p = make_obj("temp_object", *n_names, 0, in_length, *n_names);
  fill_obj(p, 0, 0, in_length, idummy, ddummy, name_list);
  doom_compress(p->a_char, &in_length, &out_length);
  for (j = 0; j < *n_names; j++) p->a_int[j] = j;
  p->c_int = p->c_obj = p->l_int;
  make_names(p);
  while (k != 0)    /* sort loop */
    {
     k = 0;
     for (i = 0; i < p->c_obj; i++)
       {
        for (j = i+1; j < p->c_obj; j++)
	  {
           if (strcmp(p->names[p->a_int[i]], p->names[p->a_int[j]]) > 0)
	     {
              temp = p->a_int[i]; p->a_int[i] = p->a_int[j]; 
              p->a_int[j] = temp; k = 1;
	     }
	  }
       }
    }
  for (j = 0; j < p->c_obj; j++)
    {
     strcpy(c, p->names[p->a_int[j]]); 
     c+= strlen(p->names[p->a_int[j]]) + 1;
     positions[j] = p->a_int[j];
    }
  delete_obj(p);
}

/* C start  doom_psequ: Put expanded sequence */
void doom_psequ            /* stores sequence table */
         (char* sequ_name, /* input: sequence name */
          int*  n_item,    /* input: actual no. of elements */
          int*  name_l,    /* no. of characters/name */
          char* names,     /* element names */
          int*  occ,       /* array occ(max_item): occ. counts */
          double* pos)     /* array pos(max_item): centre pos. */
/* C end comment */
/* F start doom_psequ

subroutine doom_psequ(name, nelement, nch, names, occur, pos)

Stores a USEd sequence. A sequence currently in USE
will be stored automatically at the end of the special MAD.

Input        (character)          (sequence) name
Input        (integer)            nelement
Input        (integer)            nch = no. of characters/name in names
Input        (character array)    (element) names
Output       (integer array)      occur(ence count of element)
Output       (d.p. array)         centre position

Remark: the last element must give the end position (marker..)

The arrays must be declared as follows (minimum sizes):

integer occur(nelement)
character *(nch) names(nelement)
double precision pos(nelement)

F end comment */
{
  int idleng = *n_item, icleng = *name_l * *n_item;
  char rec_key[KEY_LENGTH], lkey[KEY_LENGTH], type[] = "EX_SEQUENCE";
  struct object* p;
  int key_list[2] = {1, SEQU_CODE};

#ifdef _TREE_
       puts("++++++ enter doom_psequ");
#endif
  mycpy(lkey, sequ_name);
  make_key(lkey, key_list, rec_key);
  p = doom_fetch(rec_key);
  if (p != NULL)
    {
     if (*n_item > p->l_int)  grow_obj(p, *n_item, 0, 0, *n_item);
     if (idleng > p->l_dble)  grow_obj(p, 0, idleng, 0, 0);
     if (icleng > p->l_char)  grow_obj(p, 0, 0, icleng, 0);
    }
  else
    {
     p = make_obj(rec_key, *n_item, idleng, icleng, *n_item);
     doom_getenv("USE", p->par_name);
     strcpy(p->obj_type, type);
    }
  fill_obj(p, *n_item, idleng, icleng, occ, pos, names);
  doom_comp_max(p->a_char, &icleng, name_l, &p->c_char);
  p->c_obj = p->l_obj;
  doom_save(p);
}

/* C start doom_pstring: Put string */
void doom_pstring           /* stores string */
         (char* name,       /* input: key name */
          int* key_list,    /* input */
          int* sleng,       /* input: string length */
          char* string)     /* input */
/* C end comment */
/* F start doom_pstring

subroutine doom_pstring(name, key_list, nchar, string)

Stores a string.

Input        (character)          name
Input        (integer array)      key_list
Input        (integer)            nchar - length of string
Input        (character)          string

F end comment */
{
  char rec_key[KEY_LENGTH], l_key[KEY_LENGTH], type[] = "STRING";
  struct object* p;
  int j = *sleng;

#ifdef _TREE_
       puts("++++++ enter doom_pstring");
#endif
  mycpy(l_key, name);
  make_key(l_key, key_list, rec_key);
  strncpy(cdummy, string, j);
  cdummy[j] = '\0'; j++;
  p = doom_fetch(rec_key);
  if (p == NULL)
    {
     p = make_obj(rec_key, 1, 0, j, 0);
     strcpy(p->par_name, type);
     strcpy(p->obj_type, type);
     p->a_int[0] = 1;
     p->c_int = 1;
     strcpy(p->a_char, cdummy);
     p->c_char = --j;
    }
  else
    {
     if (j > p->l_char)  grow_obj(p, 0, 0, j, 0);
     strcpy(p->a_char, cdummy);
     p->c_char = --j;
    }
  doom_save(p);
}

/* C start doom_psumm: Put summary record */
void doom_psumm            /* stores summary tables = array of double */
         (char* surv_name, /* input: summary name */
          int*  n_item,    /* input: actual no. items */
          double* table)
/* C end comment */
/* F start doom_psumm

subroutine doom_psumm(name, ntable, table)

Stores a summary table.

Input        (character)          name
Input        (integer)            ntable - length of table
Input        (d.p.)               table (of d.p. values)

F end comment */
{
  char rec_key[KEY_LENGTH], type[] = "SUMMARY";
  struct object* p;
  int j;

#ifdef _TREE_
       puts("++++++ enter doom_psumm");
#endif
  mycpy(rec_key, surv_name);
  p = make_obj(rec_key, 0, *n_item, 0, 0);
  strcpy(p->par_name, type);
  strcpy(p->obj_type, type);
  for (j = 0; j < *n_item; j++) p->a_dble[j] = table[j];
  p->c_dble = *n_item;
  doom_save(p);
}

/* C start doom_psurv: Put survey table */
void doom_psurv            /* stores survey table */
         (char* surv_name, /* input: survey name */
          int*  n_item,    /* input: actual no. of elements */
          int*  name_l,    /* no. of characters/name */
          char* names,     /* element names */
          int*  occ,       /* array occ(max_item): occ. counts */
          double* pos)     /* array pos(max_item,7):
                              s, x, y, z, theta, phi, psi */
/* C end comment */
/* F start doom_psurv

subroutine doom_psurv(name, nelem, nch, elements, occur, pos)

Stores a survey table.

Input        (character)          name (MAD default: SURVEY)
Input        (integer)            nelem - length of table
Input        (character array)    elements (names)
Input        (integer array)      occur(rence counts)
Input        (d.p.array)          table (of d.p. values)

F end comment */
{
  int idleng = 7 * *n_item, icleng = *name_l * *n_item;
  char rec_key[KEY_LENGTH], type[] = "SURVEY";
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter doom_psurv");
#endif
  mycpy(rec_key, surv_name);
  p = make_obj(rec_key, *n_item, idleng, icleng, *n_item);
  strcpy(p->par_name, type);
  strcpy(p->obj_type, type);
  fill_obj(p, *n_item, idleng, icleng, occ, pos, names);
  doom_comp_max(p->a_char, &icleng, n_item, &p->c_char);
  p->c_obj = p->l_obj;
  doom_save(p);
}

/* C start  doom_ptbody: Put table (TWISS, SURVEY etc.)*/
void doom_ptbody           /* stores a table (TWISS, SURVEY etc.) */
         (char* tabl_name, /* input: table name */
          char* sequ_name, /* input: sequence name */
          int*  n_column,  /* input: no. of columns in table */
          int*  n_row,     /* input: no. of rows = no. of positions */
          int* name_l,     /* input: no. of characters/name */
          char* names,     /* input: element names */
          int*  occ,       /* input: occ. counts */
          double* optics_t) /* input: table (see twiss_att_names) */
/* C end comment */
/* F start doom_ptbody

subroutine doom_ptbody(tabl_name, sequ_name, n_column, n_row,
                       name_l, names, occur, optics_t)

Stores a table (TWISS, SURVEY etc.)

Input        (character)          tabl_name
Input        (character)          sequ_name
Input        (integer)            n_column
Input        (integer)            n_row
Input        (integer)            name_l = no. of characters/name in names
Input        (character array)    (element) names
Input        (integer array)      occur(ence count of element)
Input        (d.p. array)         optics_t
F end comment */
{
  int idleng = *n_column * (*n_row), icleng = *name_l * (*n_row);
  char rec_key[KEY_LENGTH], lkey[KEY_LENGTH], type[] = "SUB_TABLE";
  struct object* p;
  int key_list[3] = {2, TABLE_BODY, 0};

#ifdef _TREE_
       puts("++++++ enter doom_ptbody");
#endif
  key_list[2] = twiss_count++;
  mycpy(lkey, tabl_name);
  make_key(lkey, key_list, rec_key);
  p = make_obj(rec_key, *n_row, idleng, icleng, *n_row);
  strcpy(p->obj_type, type);
  mycpy(p->par_name, sequ_name);
  fill_obj(p, *n_row, idleng, icleng, occ, optics_t, names);
  doom_compress(p->a_char, &icleng, &p->c_char);
  p->c_obj = p->l_obj;
  doom_save(p);
}

/* C start doom_pthead: Stores table header (TWISS, SURVEY etc.) */
void doom_pthead    /*  Stores table header (TWISS, SURVEY etc.) */
            (char* tabl_name, /* input: table name */
             char* sequ_name, /* input: sequence name */
             int* n_delta,    /* input: # of delta values*/
             int* n_row,      /* input: # of rows per table (per delta) */
             int* sum_leng,   /* input: length of summ. record (per delta) */
             double* sum_rec, /* input: summary record */
             int* n_col,      /* input: # of columns */
             int* pos_flag,   /* input: element position flag:
				 1 start, 2 centre, 3 end */
             char* col_names, /* input: column names (compressed) */
             int* col_pos)    /* input: column positions in table */
/* C end comment */
/* F start doom_pthead

subroutine doom_pthead(tabl_name, sequ_name, n_delta, n_row, 
                        sum_leng, sum_rec, n_col, pos_flag, 
                        col_names, col_pos)

Stores table header (TWISS, SURVEY etc.)

Input        (character)         tabl_name
Input        (character)         sequ_name
Input        (integer)           n_delta (# of delta values)
Input        (integer)           n_row  = # of rows per delta
Input        (integer)           sum_leng = summary record length per delta
Input        (d.p. array)        sum_rec = summary record
Input        (integer)           n_col = # of columns
Input        (integer)           pos_flag = element position flag:
                                 1 start, 2 centre, 3 end
Input        (character)         col_names (compressed!) column names
Input        (integer array)     col_pos = column positions (of col_names)

The summary record contains 17 values per delta, i.e. a convenient declaration
would be
         double precision sum_rec(17,25)

(25 is the maximum number of deltas MAD-8 admits)
F end comment */
{
  char rec_key[KEY_LENGTH], lcp[KEY_LENGTH], type[] = "TABLE_HEADER";
  char* cs = col_names;
  int sum_l = *n_delta * (*sum_leng), key_list[2] = {1, TABLE_HEAD},
      j = 0, nch, nc = *n_col;
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter doom_pthead");
#endif
  mycpy(lcp, tabl_name);
  make_key(lcp, key_list, rec_key);
  while (j < nc)
    {
     if (*cs++ == '\0') j++;
    }
  if (j != nc)
    {
     printf("%s%s\n%s %d  %s %d\n", " <<< DOOM >>> warning: object ", 
     rec_key, " no. of columns:", nc, "no. of column names:", j);
    }
  nch = cs - col_names;
  p = make_obj(rec_key, OPTICS_INTL + nc, sum_l, nch, j);
  strcpy(p->obj_type, type);
  mycpy(p->par_name, sequ_name);
  fill_obj(p, nc, sum_l, nch, col_pos, sum_rec, col_names);
  p->a_int[nc] = *n_delta;
  p->a_int[nc+1] = *n_row;
  p->a_int[nc+2] = *sum_leng;
  p->a_int[nc+3] = *pos_flag;
  p->c_int = p->l_int;
  p->c_obj = p->l_obj;
  make_names(p);
  doom_save(p);
}

/* C start doom_ptrack: Put track master record */
void doom_ptrack    /* puts track survival table */
            (int*  locc,   /* input: number of list items */
             int*  list)   /* input: list address */
/* C end comment */
/* F start doom_ptrack

subroutine doom_ptrack(length, list)

Stores the table containing the list of turns with results saved.

Input        (integer)           length = maximum/actual no. turns in list
Input        (integer array)     (2,*): list (of turn numbers,
                                   and # of surving particle)

F end comment */
{
  char rec_key[KEY_LENGTH], name[] = "TRACK_SURVIVAL", type[] = "SURVIVAL";
  int key_list[3] = {1, TRACK_SV}, occ = 2 * *locc;
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter doom_ptrack");
#endif
  make_key(name, key_list, rec_key);
  p = make_obj(rec_key, occ, 0, 0, 0);
  strcpy(p->par_name, type);
  strcpy(p->obj_type, type);
  fill_obj(p, occ, 0, 0, list, ddummy, cdummy);
  doom_save(p);
}

/* C start doom_pturn: Put one turn record */
void doom_pturn(    /* saves the track data of one turn */
     int* turn,     /* input: turn number */
     int* number,   /* input: number of particles */
     int* list,     /* input: particle list */
     double* track) /* input: track coordinates (6,...) */
/* C end comment */
/* F start doom_pturn

subroutine doom_pturn(turn, npart, list, coords)

Stores the particle numbers and their coordinates for a specific turn.

Input        (integer)            turn (number)
Input        (integer)            npart = number of particles in list
Input        (integer array)      list (of particle numbers)
Input        (d.p. array)         coords = particle coordinates

Declarations required:

integer list(..)
double precision coords(6,..)

F end comment */
{
  struct object* p;
  int docc = 6 * *number;
  char rec_key[KEY_LENGTH], lkey[KEY_LENGTH];
  int key_list[3] = {1, TRACK_RECORD};
  char turn_num[16], name[] = "TURN", type[] = "TURN";

#ifdef _TREE_
       puts("++++++ enter doom_pturn");
#endif
  sprintf(turn_num, "%d", *turn);
  make_dollar_key(name, turn_num, lkey);
  make_key(lkey, key_list, rec_key);
  p = make_obj(rec_key, *number, docc, 0, 0);
  strcpy(p->par_name, type);
  strcpy(p->obj_type, type);
  fill_obj(p, *number, docc, 0, list, track, cdummy);
  doom_save(p);
}

/* C start doom_ptype: Put object type */
void doom_ptype       /* puts type into object */
     (char* name,     /* input: object name */
      int* key_list,  /* input: key_list */
      char* type)     /* input: type */
/* C end comment */
/* F start doom_ptype

subroutine doom_ptype(name, key_list, type)

Stores the parent name of an object.

Input        (character)          name
Input        (integer array)      key_list
Input        (character)          type

F end comment */
{
  char rec_key[KEY_LENGTH], lkey[KEY_LENGTH];
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter doom_ptype");
#endif
  mycpy(lkey, name);
  make_key(lkey, key_list, rec_key);
  p = doom_fetch(rec_key);
  if (p != NULL) mycpy(p->obj_type, type);
}

/* C start doom_recall: Get (any) array collection */
void doom_recall(    /* recalls an object from memory or the d.b.
		      if it does not exist, all length = -1 */
     char* name,     /* input */
     int*  key_list, /* input: key token list: number + codes */
     int* lint,      /* input/output: max/length of integer array */
     int* ldble,     /* input/output: max/length of double array */
     int* lchar,     /* input/output: max/length of char array */
     int* a_int,     /* output: integer array */
     double* a_dble, /* output: d.p. array */
     char* a_char)   /* output: string */
/* C end comment */
/* F start doom_recall

subroutine doom_recall(name, keylist,
                       lint, ldble, lchar, aint, adble, achar)

Returns any type of database object.

Input        (character)          name
Input        (integer array)      keylist
Input/Output (integer)            lint: max./actual number of integers in aint
Input/Output (integer)            ldble: max./actual number of values in adble
Input/Output (integer)            lchar: max./actual number of char.s in achar
Output       (integer array)      aint
Output       (d.p. array)         adble
Output       (character)          achar

F end comment */
{
  char rec_key[KEY_LENGTH], lkey[KEY_LENGTH];
  struct object* p;
  int lp, index, j;

#ifdef _TREE_
       puts("++++++ enter doom_recall");
#endif
  mycpy(lkey, name);
  make_key(lkey, key_list, rec_key);
  lp = n_list_pos(rec_key, o_list_sys[0][0]);
  if (lp < 0)  /* object exists */
    {
     index = o_list_sys[0][0]->a_int[-lp-1];
     p = o_list_sys[0][0]->p_obj[index]; 
     if (p == NULL)
       {
        if (db_is_open() == 0)
          {
           puts("+++DOOM - doom_recall: db not open"); 
           *lint = -1; *ldble = -1; *lchar = -1; return;
          }
        else
          {
           if (db_exists(rec_key))
             {
              p = get_obj(rec_key);
             }
           else
	     {
              *lint = -2; *ldble = -2; *lchar = -2; return;
	     }
           o_list_sys[0][0]->p_obj[index] = p;
	  }
       }
    }
  else
    {
     *lint = -1; *ldble = -1; *lchar = -1; return;
    }
  if (*lint > p->c_int) *lint = p->c_int; 
  if (*ldble > p->c_dble) *ldble = p->c_dble; 
  if (*lchar > p->c_char) *lchar = p->c_char;
  for (j = 0; j < *lint;  j++)   a_int[j] = p->a_int[j];
  for (j = 0; j < *ldble; j++)  a_dble[j] = p->a_dble[j];
  for (j = 0; j < *lchar; j++)  a_char[j] = p->a_char[j];
}

/* C start doom_rline:     rewinds line buffer object */
void doom_rline          /* "rewinds" line buffer object */
        (char* name)     /* input: object name */
/* C end comment */
/* F start doom_rline

subroutine doom_rline(name)

Rewinds a line buffer object

Input        (character)          (object) name

F end comment */
{
  char rec_key[KEY_LENGTH];
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter doom_rline");
#endif
  mycpy(rec_key, name);
  p = doom_fetch(rec_key);
  if (p != NULL)
    {
     if (p->c_int > 2)  p->a_int[1] = p->a_int[2] = 0;
    }
  else
    {
     printf("<<< DOOM >>> warning - file object: %s not found\n", rec_key);
     exit(1);
    }
}

/* C start doom_robj: Get object (from Fortran) */
void doom_robj       /* restores object (for call from fortran) */
     (char* name,    /* input:  name */
      int* key_list, /* input: key list for composed name */
      char* parent,  /* output: */
      char* type,    /* output: */
      int*  n_int,   /* input/output: max/act. number of integers */
      int*  n_dble,  /* input/output: max/act. number of doubles */
      int*  n_name,  /* input/output: max/number of names in char array */
      int*  ncname,  /* input: (fixed)number of char.s per name */
      int*  a_int,   /* output: integer array */
      double*  a_dble,/* output: double array */
      char* a_char)  /* output: char array for names (blank padded) */
/* C end comment */
/* F start doom_robj

subroutine doom_robj(name, keylist, parent, type, nint, ndble, nname,
                     nch, aint, adble, achar)

Restores an object (for Fortran calls)

Input        (character)          (object) name
Input        (integer array)      keylist   
Output       (character * nch)    parent 
Output       (character * nch)    type  
Input/Output (integer)            nint = # integers
Input/Output (integer)            ndble = # d.p.
Input/Output (integer)            nname = # names
Input        (integer)            nch = char./character variable
Output       (integer array)      aint (integers)
Output       (d.p. array)         adble (d.p. values)
Output       (character array)    achar = names (blank padded)

declaration for achar:

            character * (nch)  achar(nname)
 
F end comment */
{
  char lkey[KEY_LENGTH], rec_key[KEY_LENGTH], tmp[KEY_LENGTH];
  struct object* p;
  int j, c_all, nloc;
  char* tt;

#ifdef _TREE_
       puts("++++++ enter doom_robj");
#endif
  mycpy(lkey, name);
  mycpy(tmp, type);
  make_key(lkey, key_list, rec_key);
  p = doom_fetch(rec_key);
  if (p == NULL)
    {
     *n_int = -1; *n_dble = -1; *n_name = -1;
    }
  else
    {
     if (strcmp(p->obj_type, tmp) != 0)
       {
        printf("%s%s\n%s%s%s%s\n", " <<< DOOM >>> warning: object ",
        rec_key, " has type: ", p->obj_type, " and not type: ", tmp);
       }
     if (*n_int > p->c_int) *n_int = p->c_int; 
     if (*n_dble > p->c_dble) *n_dble = p->c_dble; 
     if (*n_name > p->c_obj) *n_name = p->c_obj; 
     for (j = 0; j < *n_int; j++)  a_int[j] = p->a_int[j];
     for (j = 0; j < *n_dble; j++)  a_dble[j] = p->a_dble[j];
     c_all = *ncname * p->c_obj;
     if (c_all > max_alloc)
       {
        printf(" <<< DOOM >>> fatal:\n memory allocation request = %d%s%d\n",
               c_all, " above limit = ", max_alloc);
        exit(1);
       }
     tt = (char*) malloc(c_all);
     for (j = 0; j < c_all; j++)  tt[j] = p->a_char[j];
     doom_blow(tt, &p->c_char, ncname, &nloc);
     if (*n_name > nloc)  *n_name = nloc;
     c_all = *ncname * *n_name;
     for (j = 0; j < c_all; j++)  a_char[j] = tt[j];
     free(tt);
    }
}

/* C start doom_save: Put object (inverse of doom_fetch) */
void doom_save(  /* puts object into list to be saved at end.
                    If it exists, the previous object is deleted
                    if it has a different pointer */
      struct object* psave) /* input: object to be saved */
/* C end comment */
{
  int lp, index, j;

#ifdef _TREE_
       puts("++++++ enter doom_save");
#endif
  if (psave!= NULL)
  {
   my_time();
   psave->ma_time = major_time; psave->mi_time = minor_time;
   for (j = 0; j < psave->c_char; j++)
    {
     if (psave->a_char[j] == '|') psave->a_char[j] = '\0';
    }
   lp = n_list_pos(psave->key, o_list_sys[0][0]);
   if (lp < 0)  /* object exists */
    {
     index = o_list_sys[0][0]->a_int[-lp-1];
     if (o_list_sys[0][0]->p_obj[index] != psave)
       { 
        delete_obj(o_list_sys[0][0]->p_obj[index]);
        o_list_sys[0][0]->p_obj[index] = psave;
       }
    }
   else
    {
     o_list_add(lp, psave);
    }
   if (debug_flag > 3) printf("doom_save major + minor time = %d %d\n", 
       psave->ma_time, psave->mi_time); 
  }
}

/* C start doom_setenv: Set value of environment variable */
void doom_setenv(  /* sets new value of environment variable */
     char* variable, /* input */
     char* value)    /* input */
/* C end comment */
/* F start doom_setenv

subroutine doom_setenv(env, value)

Sets environment variable (only valid while job is running)

Input        (character)          env  (must match an existing variable)
Input        (character)          value (character value of env)

F end comment */
{
  char l1key[KEY_LENGTH], l2key[KEY_LENGTH];
  int lp, index;

#ifdef _TREE_
       puts("++++++ enter doom_setenv");
#endif
  mycpy(l1key, variable);
  lp = n_list_pos(l1key, o_list_sys[1][0]);
  if (lp < 0)
    {
     mycpy(l2key, value);
     index = o_list_sys[1][0]->a_int[-lp-1]; 
     n_list_replace(index, l2key, o_list_sys[1][1]);
    }
}

/* C start doom_setvar: Set value of internal variable */
void doom_setvar(  /* set DOOM internal variable (ad hoc) */
                 char* name,
                 int*  value)
/* C end comment */
/* F start doom_setvar

subroutine doom_setvar(name, value)

Sets internal variable (only valid while job is running)

Input        (character)          name  (must match an existing variable)
Input        (character)          value (integer value of name)

F end comment */
{
  char lkey[KEY_LENGTH];
  int j;

  mycpy(lkey, name);
  for (j = 0; j < n_int_var; j++)
    {
     if (strcmp(lkey, l_int_var[j]) == 0) *p_int_var[j] = *value;
    }
}

/* C start doom_sobj: Put object (from Fortran) */
void doom_sobj       /* stores object (from fortran) */
     (char* name,    /* input:  name */
      int* key_list, /* input: key list for composed name */
      char* parent,  /* input: */
      char* type,    /* input: */
      int*  n_int,   /* input: number of integers */
      int*  n_dble,  /* input: number of doubles */
      int*  n_name,  /* input: number of names in char array */
      int*  ncname,  /* input: (fixed)number of char.s per name */
      int*  a_int,   /* input: integer array */
      double*  a_dble,/* input: double array */
      char* a_char)  /* input: char array for names */
/* C end comment */
/* F start doom_sobj

subroutine doom_sobj(name, keylist, parent, type, nint, ndble, nname,
                     nch, aint, adble, achar)

Stores an object (for Fortran calls)

Input        (character)          (object) name
Input        (integer array)      keylist  
Input        (character * nch)    parent 
Input        (character * nch)    type    
Input        (integer)            nint = # integers
Input        (integer)            ndble = # d.p.
Input        (integer)            nname = # names
Input        (integer)            nch = char./character variable
Input        (integer array)      aint (integers)
Input        (d.p. array)         adble (d.p. values)
Input        (character array)    achar = names

declaration for achar:

            character * (nch)  achar(nname)
 
F end comment */
{
  char lkey[KEY_LENGTH], rec_key[KEY_LENGTH], tmp[KEY_LENGTH];
  struct object* p;
  int c_char = *ncname * *n_name;

#ifdef _TREE_
       puts("++++++ enter doom_sobj");
#endif
  mycpy(lkey, name);
  mycpy(tmp, type);
  make_key(lkey, key_list, rec_key);
  p = doom_fetch(rec_key);
  if (p != NULL)
    {
     if (strcmp(p->obj_type, tmp) != 0)
       {
        printf("%s%s\n%s%s%s%s\n", " <<< DOOM >>> warning: object ",
        rec_key, " has type: ", p->obj_type, " and not type: ", tmp);
       }
     if (*n_int > p->l_int || *n_dble > p->l_dble
         || c_char > p->l_char || *n_name > p->l_obj)
        { 
         grow_obj(p, *n_int, *n_dble, c_char, *n_name);
        }
    }
  else
    {
     p = make_obj(rec_key, *n_int, *n_dble, c_char, *n_name);
     strcpy(p->obj_type, tmp);
    }
  fill_obj(p, *n_int, *n_dble, c_char, a_int, a_dble, a_char);
  p->c_obj = p->l_obj;
  mycpy(tmp, parent);
  strcpy(p->par_name, tmp);
  doom_compress(p->a_char, &c_char, &p->c_char);
  doom_save(p);
}

void doom_snoup(  /* set no_update flag */
                 int* flag)
{
  no_update = *flag;
}

/* C start doom_store: Put (any) array collection */
void doom_store(     /* creates an object in memory and in the d.b.
                        If it exists, it is replaced. */
     char* name,     /* input */
     int*  key_list, /* input: key token list: number + codes */
     int* lint,      /* input: length of integer array */
     int* ldble,     /* input: length of double array */
     int* lchar,     /* input: length of char array */
     int* a_int,     /* input: integer array */
     double* a_dble, /* input: d.p. array */
     char* a_char)   /* input: string */
/* C end comment */
/* F start doom_store

subroutine doom_store(name, keylist,
                      lint, ldble, lchar, aint, adble, achar)

Stores any type of database object.

Input        (character)          name
Input        (integer array)      keylist
Input        (integer)            lint = number of integers in aint
Input        (integer)            ldble = number of d.p. values in adble
Input        (integer)            lchar = number of characters in achar
Input        (integer array)      aint
Input        (d.p. array)         adble
Input        (character)          achar

F end comment */
{
  char rec_key[KEY_LENGTH], lkey[KEY_LENGTH];
  struct object* p;
  int lp, index, j;

#ifdef _TREE_
       puts("++++++ enter doom_store");
#endif
  mycpy(lkey, name);
  make_key(lkey, key_list, rec_key);
  lp = n_list_pos(rec_key, o_list_sys[0][0]);
  if (lp < 0)  /* object exists */
    {
     index = o_list_sys[0][0]->a_int[-lp-1];
     p = o_list_sys[0][0]->p_obj[index];
     if (p == NULL)
       {
        if (db_is_open() == 0)
          {
           puts("+++DOOM - doom_store: db not open"); exit(1);
          }
        else
          {
           if (db_exists(rec_key))
            {
             p = get_obj(rec_key);
	    }
           else
            {
             printf("%s%s%s\n", 
             "DOOM fatal - doom_store:",
             " object in system bank, not in db: ", rec_key);
             exit(1);
            }
	  }
        o_list_sys[0][0]->p_obj[index] = p;
       }
     if (p->l_int < *lint || p->l_dble < *ldble || p->l_char < *lchar)
     grow_obj(p, *lint, *ldble, *lchar, 0);
    }
  else
    {
     p = make_obj(rec_key, *lint, *ldble, *lchar, 0);
     o_list_add(lp, p);
    }
  my_time();
  p->ma_time = major_time; p->mi_time = minor_time;
  if (debug_flag > 3) printf("doom_store major + minor time = %d %d\n",
      p->ma_time, p->mi_time); 
  for (j = 0; j < *lint; j++)   p->a_int[j] = a_int[j];
  for (j = 0; j < *ldble; j++)  p->a_dble[j] = a_dble[j];
  for (j = 0; j < *lchar; j++)  p->a_char[j] = a_char[j];
  p->c_int = *lint; p->c_dble = *ldble; p->c_char = *lchar;
}

void doom_swap(     /* swaps one int if flag set, else does nothing */
               int* itoswap)   /* integer to swap */
{
  if (swap_flag) swap_int(1, itoswap);
}

/* C start doom_time: Get current DOOM time */
void doom_time(double* r_time)    /* external time in d.p. */
/* C end comment */
/* F start doom_time

subroutine doom_time(time)

Returns current DOOM time used to stamp objects

Output   time  (d.p.)      time (sec.s since 1.1.97 + microsteps)

F end comment */
{
  *r_time =  (double) major_time + (double) minor_time * micro_step;
}

void doom_wtree(   /* recursive listing of tree banks */
             int* iq,   /* zebra iq */
             int* lq,   /* zebra lq */
             int* lb)   /* top bank address as value */

{
  int j, lbl;

#ifdef _TREE_
       puts("++++++ enter doom_wtree");
#endif
  if (*lb > 0)
    {
     lbl = *lb;
     printf("%s%d%s%d%s%d\n", "   ns: ", iq[lbl-3],
            "   nl: ", iq[lbl-4] - iq[lbl-3], "   nd: ", iq[lbl-2]);
     for (j = 0; j <= iq[*lb-3]; j++)
          doom_wtree(iq, lq, &lq[--lbl]);
    }
}

void fill_obj( /* copies integer, double, and char arrays into object */
              struct object* p, /* input: object pointer */
              int c_i, /* input: no. of integers */
              int c_d, /* input: no. of doubles */
              int c_c, /* input: length of char array */
              const int* a_int, /* input: integer array */
              const double* a_dble, /* input: double array */
              const char* a_char) /* input: char array */
{
  int j;

  p->c_int  = c_i > p->l_int ? p->l_int : c_i; 
  p->c_dble = c_d > p->l_dble ? p->l_dble : c_d; 
  p->c_char = c_c > p->l_char ? p->l_char : c_c;
  for (j = 0; j < p->c_int; j++)  p->a_int[j] = a_int[j];
  for (j = 0; j < p->c_dble; j++)  p->a_dble[j] = a_dble[j];
  for (j = 0; j < p->c_char; j++)  p->a_char[j] = a_char[j];
}

int flat_list(    /* flattens char** array, returns length */
              char* list[],    /* input: list */
              int length,      /* input: length of list */
              char*  f_list)   /* flat list */
{
  int j, n = 0;
  char* c1; 
  char* c2; 

#ifdef _TREE_
       puts("++++++ enter flat_list");
#endif
  c1 = f_list;
  for (j = 0; j < length; j++)
    {
     c2 = list[j];
     while (*c2 != '\0') {*c1++ = *c2++; n++;}
     *c1++ = '\0'; n++;
    }
  return n;
}

void free_obj(struct object* p)  /* deletes object, but not sys entry */
{

#ifdef _TREE_
       puts("++++++ enter free_obj");
#endif
  if (p != NULL)
    {
     if (p->a_int != NULL)  free(p->a_int);
     if (p->a_dble != NULL)  free(p->a_dble);
     if (p->a_char != NULL)  free(p->a_char);
     if (p->p_obj != NULL)  free(p->p_obj);
     if (p->names != NULL)  free(p->names);
     free(p); p = NULL;
    }
}

struct object* get_obj(            /* gets an object from db */
                       char* key)  /* object key */
{
  int counts[4] = {CONT_OBJ, 0, 0, 0};
  char rec_key[KEY_LENGTH];
  struct object* p;
  int c_all;

#ifdef _TREE_
       puts("++++++ enter get_obj");
#endif
  mycpy(rec_key, key);
  p = make_obj(rec_key, 0, 0, 0, 0);
/* first get array sizes */
  if (start_get(rec_key, 0, counts, &p->ma_time, idummy, 
      ddummy, cdummy) == 0)
    {
     printf(" <<< DOOM >>> fatal: cannot retrieve record %s\n", rec_key);
     exit(1);
    }
/* create arrays with correct sizes */
  if (p->l_int > 0)
    {
     c_all = p->l_int * sizeof(int);
     if (c_all > max_alloc)
       {
        printf(" <<< DOOM >>> fatal:\n memory allocation request = %d%s%d\n",
               c_all, " above limit = ", max_alloc);
        exit(1);
       }
     p->a_int = (int*) malloc(c_all);
    }
  if (p->l_dble > 0) p->a_dble = (double*) malloc(p->l_dble * sizeof(double));
  if (p->l_char > 0) p->a_char = (char*) malloc(p->l_char);
  if (p->l_obj > 0)
    {
     p->p_obj = (struct object**) calloc(p->l_obj, sizeof(struct object*));
     p->names = (char**) calloc(p->l_obj, sizeof(char*));
    }
  counts[0] = 0;
  counts[1] = p->c_int; counts[2] = p->c_dble; 
  counts[3] = p->c_char;
  rec_get(counts, &p->ma_time, p->a_int, p->a_dble, p->a_char);
  if (p->l_obj > 0)  make_names(p);
  return p;
}

struct object* get_obj_list()
{
  return o_list_sys[0][0];
}

int get_text(    /* gets one text record from file */
             int offset,      /* in bytes */
             int size,        /* in bytes */
             char* io_buff)
{
  size_t l_read, l_asked;

#ifdef _TREE_
       puts("++++++ enter get_text");
#endif

  if (fseek(db_stream, (long) offset, SEEK_SET) != 0) return -4;  
  l_asked = size;
  l_read = fread(io_buff, 1L, l_asked, db_stream);
  if (feof(db_stream) != 0)  return -3;
  if (ferror(db_stream) != 0)  return -2;
  return (int) l_read;
}

void grow_obj(   /* augments object arrays to new specified values */
              struct object* p,
              int new_i,  /* new int array length */
              int new_d,  /* new double array length */
              int new_c,  /* new char array length */
              int new_o)  /* new object array length */
{
  struct object** p_obj;
  char* a_char;
  double* a_dble;
  int* a_int;
  char** names;
  int j, c_all, cshift = 0;

#ifdef _TREE_
       puts("++++++ enter grow_obj");
  printf("grow-obj %s with %d %d %d %d\n", p->key, new_i, new_d, new_c, new_o);
#endif
  if (p->l_int < new_i)
    {
     a_int = p->a_int;
     p->l_int = new_i;
     c_all = p->l_int * sizeof(int);
     if (c_all > max_alloc)
      {
       printf(" <<< DOOM >>> fatal:\n memory allocation request = %d%s%d\n",
              c_all, " above limit = ", max_alloc);
       exit(1);
      }
     p->a_int = (int*) malloc(c_all);
     for (j = 0; j < p->c_int; j++) p->a_int[j] = a_int[j];
     free(a_int);
    }
  if (p->l_dble < new_d)
    {
     a_dble = p->a_dble;
     p->l_dble = new_d;
     c_all = p->l_dble * sizeof(double);
     if (c_all > max_alloc)
      {
       printf(" <<< DOOM >>> fatal:\n memory allocation request = %d%s%d\n",
              c_all, " above limit = ", max_alloc);
       exit(1);
      }
     p->a_dble = (double*) malloc(c_all);
     for (j = 0; j < p->c_dble; j++) p->a_dble[j] = a_dble[j];
     free(a_dble);
    }
  if (p->l_char < new_c)
    {
     cshift = 1;
     a_char = p->a_char;
     p->l_char = new_c;
     c_all = p->l_char;
     if (c_all > max_alloc)
      {
       printf(" <<< DOOM >>> fatal:\n memory allocation request = %d%s%d\n",
              c_all, " above limit = ", max_alloc);
       exit(1);
      }
     p->a_char = (char*) malloc(c_all);
     for (j = 0; j < p->c_char; j++) p->a_char[j] = a_char[j];
     free(a_char);
    }
  if (p->l_obj < new_o)
    {
     p_obj = p->p_obj;
     p->l_obj = new_o;
     c_all = p->l_obj * sizeof(struct object*);
     if (c_all > max_alloc)
      {
       printf(" <<< DOOM >>> fatal:\n memory allocation request = %d%s%d\n",
              c_all, " above limit = ", max_alloc);
       exit(1);
      }
     p->p_obj = (struct object**) calloc(p->l_obj, sizeof(struct object*));
     for (j = 0; j < p->c_obj; j++) p->p_obj[j] = p_obj[j];
     names = p->names;
     p->names = (char**) calloc(p->l_obj, sizeof(char*));
     for (j = 0; j < p->c_obj; j++) p->names[j] = names[j];
     free(p_obj); free(names);
    }
  if (cshift > 0 && p->l_obj > 0)  make_names(p);
}

int i_list_add(int number, struct object* p)
{
  int lp, j;

#ifdef _TREE_
       puts("++++++ enter i_list_add");
#endif
  lp = i_list_pos(number, p);
  if (lp > -1)
    {
     if (p->c_int == p->l_int)  grow_obj(p, 2*p->l_int, 0, 0, 0);
     for (j = p->c_int; j > lp; j--)  p->a_int[j] = p->a_int[j-1];
     p->a_int[lp] = number;
     p->c_int++;
    }
  return lp;
}

int i_list_pos(int number, struct object* p)
{
  int mid, low = 0, last = 0, high = p->c_int - 1;

#ifdef _TREE_
       puts("++++++ enter i_list_pos");
#endif
  while (low <= high)
    {
     mid = (low + high) / 2;
     if ( number < p->a_int[mid])  {high = mid - 1; last = mid;}
     else if (number > p->a_int[mid]) {low  = mid + 1; last = low;}
     else              return (-mid-1);
    }
    return last;
}

void make_dollar_key        /* joins two strings with a '_' sep. */
          (const char* string1,
           const char* string2,
           char* out_string)
{
  int i = 0, j = 0;

#ifdef _TREE_
       puts("++++++ enter make_dollar_key");
#endif
  while (string1[i] != '\0') out_string[j++] = string1[i++];
  out_string[j++] = '$'; i = 0; 
  while (string2[i] != '\0') out_string[j++] = string2[i++];
  out_string[j++] = '\0'; 
}
          
struct object* make_isp_table() /* creates the element/isp table from dict */
{
  char dict_name[] = "MAD-8-DICT";
  struct object* p;
  char line[100];
  char* cp[100];
  int n, length, count = 0;
  char* tp = cdummy;

#ifdef _TREE_
       puts("++++++ enter make_isp_table");
#endif
  doom_rline(dict_name);
  do
    {
     doom_gline(dict_name, line, &length);
     cp[0] = strtok(line," =:,\n"); n = 0;
     while (cp[n] != NULL)  cp[++n] = strtok(NULL, " =:,\n");
     if (n > 5 && strcmp(cp[1], "keyword") == 0 && *cp[3] == '5')
       {
        upper_conv(cp[0]); count++;
	strcpy(tp, cp[0]); tp += strlen(tp)+1; /* element key */
       }
    } while (length > 0);
  length = tp-cdummy;
  p = make_obj("isp_table", 0, 0, length, count); 
  fill_obj(p, 0, 0, length, idummy, ddummy, cdummy);
  p->c_obj = count;
  make_names(p);
  doom_save(p);
  return p;
}

void make_key               /* converts name and list into key */
          (const char* pn,        /* name */
           const int *add,  /* count + numbers */
           char* pret)      /* key */
{
  int j, ladd;
  char c[8];
  const char *pin;
  char *pout;

  pin = pn;
  pout = pret;
  while(*pin) *pout++ = *pin++;
  for (j = 1; j <= add[0]; j++)
    {
     ladd = add[j] < 0 ? -add[j] : add[j];
     *pout++ = '@';
     sprintf(c, "%d", ladd); pin = c;
     while (*pin) *pout++ = *pin++;
    }
  *pout = '\0';
  if (pout - pret > KEY_LENGTH - 1)
    {
     printf("DOOM fatal - name + keys > KEY_LENGTH for %s\n", pret);
     exit(1);
    }
}

void make_names( /* splits char array into names */
                struct object* p)
{
  int j, k = 0;
  char* c = p->a_char;

#ifdef _TREE_
       puts("++++++ enter make_names");
#endif
  if (p->l_obj > 0)
    {
     for (j = 0; j < p->c_char; j++)
       {
        if (p->a_char[j] == '\0')
          {
	   p->names[k++] = c;  c = &p->a_char[j+1];
           if (k > p->c_obj)
             {
              printf("%s%s\n",
              "DOOM make_names warning: wrong name count in object ",
              p->key);
              printf("%s%d%s%d\n", "stopped at: ", k, "  object count: ", 
              p->c_obj); 
              break;
             }
          }
       }
    }
}

struct object* make_obj(   /* creates a new object */
               char* key,
               int vlint,       /* length of integer array */
               int vldble,      /* length of double array */
               int vlchar,      /* length of char array */
               int vlpobj)      /* length of object pointer array */
{
  struct object* p;

#ifdef _TREE_
       puts("++++++ enter make_object");
#endif
  p = (struct object*)  calloc(1, sizeof(struct object));
  mycpy(p->key, key);
  if ((p->l_int = vlint) > 0) 
       p->a_int = (int*) malloc(p->l_int * sizeof(int));
  if ((p->l_dble = vldble) > 0)
       p->a_dble = (double*) malloc(p->l_dble * sizeof(double));
  if ((p->l_char = vlchar) > 0)
       p->a_char = (char*) malloc(p->l_char);
  if ((p->l_obj  = vlpobj) > 0)
    {
     p->p_obj = (struct object**) calloc(p->l_obj, sizeof(struct object*));
     p->names = (char**) calloc(p->l_obj, sizeof(char*));
    }
  p->parent = NULL;
  my_time();
  p->ma_time = major_time; p->mi_time = minor_time;
#ifdef _TREE_
       printf("++++++ exit make_object, time: %d %d\n", 
              p->ma_time, p->mi_time);
#endif
  return p;
}

int put_text(    /* puts one text record into file */
             int offset,      /* in bytes */
             int size,        /* in bytes */
             char* io_buff)
{
  size_t l_written, l_asked;

#ifdef _TREE_
       puts("++++++ enter put_text");
#endif

  if (fseek(db_stream, (long) offset, SEEK_SET) != 0) return -4;  
  l_asked = size;
  l_written = fwrite(io_buff, 1L, l_asked, db_stream);
  if (l_written < l_asked)  return -2;
  return (int) l_written;
}

void make_volatile()  /* makes tables which are not saved */
{
  int sys10_l[4];
  int j, lp;

  sys10_l[0] = environ_size;
  sys10_l[1] = 0;
  sys10_l[2] = KEY_LENGTH*environ_size;
  sys10_l[3] = environ_size;

#ifdef _TREE_
       puts("++++++ enter make_volatile");
#endif
/* environment variables */
  o_list_sys[1][0] = make_obj(system_10,
         sys10_l[0], sys10_l[1], sys10_l[2], sys10_l[3]);
  o_list_sys[1][1] = make_obj(system_11,
         sys10_l[0], sys10_l[1], sys10_l[2], sys10_l[3]);
  for (j = 0; j < environ_size; j++)
    {
     lp = n_list_pos(environ_vars[j], o_list_sys[1][0]);
     if (lp >= 0)
       {
        n_list_add(lp, environ_vars[j], o_list_sys[1][0]);
        n_list_add(lp, environ_defs[j], o_list_sys[1][1]);
       }
    }
}

void mycpy(char* sout, char* sin)
{
  char *p, *q;
  int l = 1;

  p = sin;  q = sout;
  while (*p > ' ' && *p <= '~' && l < KEY_LENGTH)  
      {
       *q++ = *p++;  l++;
      }
  *q = '\0';
}

void myscpy(char* sout, char* sin)
{
  char *p, *q;
  int l = 1;

  p = sin;  q = sout;
  while (*p > ' ' && *p <= '~' && l < C_DUM_SIZE)  
      {
       *q++ = *p++;  l++;
      }
  *q = '\0';
}

void my_time()       /* sets major_time and minor_time =
                        time in sec.s + microsteps */
{
  time_t utime;
  int ltime;

  time(&utime); ltime = (int) utime - time_vault;
  if (ltime > major_time)
    {
     major_time = ltime; minor_time = 0;
    }
  else minor_time++;
}

void n_list_add( /* adds a name to a list */
             int lp, /* list position (from n_list_pos) */
             char* str, /* string to be added */
             struct object* list)
{
  int j, sl, c, ac;

#ifdef _TREE_
       puts("++++++ enter n_list_add");
#endif
  if (lp > -1)
    {
     if (list->c_obj == list->l_obj)  
         grow_obj(list, 2*list->l_int, 0, 0, 2*list->l_obj);
     sl = strlen(str) + 1; c = list->c_char;
     if ((ac = c + sl) >= list->l_char)  grow_obj(list, 0, 0, 2*ac, 0);
     list->names[list->c_obj] = &list->a_char[c];
     if(lp < list->c_obj)
       {
        for (j = list->c_obj; j > lp; j--) list->a_int[j] = list->a_int[j-1];
        }
     strcpy(list->names[list->c_obj], str);
     list->a_int[lp] = list->c_int;
     list->c_obj = ++list->c_int; list->c_char += sl;
    }
}

int n_list_pos( /* gives -pos-1 in indexed list, or place to go if absent */
               char* str, /* lookup string */
               struct object* list)
{
  int pos, low, mid, high, last;
  low = 0; last = 0; high = list->c_obj - 1;

#ifdef _TREE_
       puts("++++++ enter n_list_pos");
#endif
  while (low <= high)
    {
     mid = (low + high) / 2;
     if ((pos = strcmp(str, list->names[list->a_int[mid]])) < 0)
                       {high = mid - 1; last = mid;}
     else if (pos > 0) {low  = mid + 1; last = low;}
     else              return (-mid-1);
    }
    return last;
}

void n_list_replace( /* replaces a name in a list */
             int index, /* real list position >= 0 (index from n_list_pos) */
             char* p, /* replacement string */
             struct object* list)
{
  int j, slold, slnew, c, d, lpc;

#ifdef _TREE_
       puts("++++++ enter n_list_replace");
#endif
  if (index > -1)
    {
     slold = strlen(list->names[index]) + 1; 
     slnew = strlen(p) + 1;
     if (slold != slnew)
       {
        d = slnew - slold; c = list->c_char; 
        if (c+d > list->l_char) grow_obj(list, 0, 0, 2*(c+d), 0);
        if (index+1 < list->c_obj)
	  {
           lpc = &list->a_char[c] - list->names[index+1];  
           shift_char(list->names[index+1], d, lpc);
           for (j = index+1; j < list->c_obj; j++) list->names[j] += d;
	  }
        list->c_char += d;
       }
     strcpy(list->names[index], p);
    }
}

int non_zero(int* list, int length)
{
  int i, nz = 0;

  for (i = 0; i < length; i++) nz += list[i] == 0 ? 0 : 1;
  return nz;
}

void o_list_add(                    /* adds an object to list 0.0 */
                 int lp,            /* list position */
                 struct object* p)
{
  int j, k, sl, space;
  int c, ac;
  struct object* l00 = o_list_sys[0][0];
  struct object* lsys[N_SYSTEM];

#ifdef _TREE_
       puts("++++++ enter o_list_add");
#endif
  if (lp > -1)
    {
     if (debug_flag > 2) prt_obj(p);
     for (k = 0; k < N_SYSTEM; k++)  lsys[k] = o_list_sys[0][k];
     if (l00->c_obj == l00->l_obj) 
       {
        space = 2*l00->l_obj; 
        grow_obj(o_list_sys[0][0], space, 0, 0, space); 
        for (k = 1; k < N_SYSTEM; k++) grow_obj(lsys[k], space, 0, 0, 0); 
       }
     sl = strlen(p->key) + 1; c = l00->c_char; ac = c + sl;
     if (ac >= l00->l_char)
       {
        grow_obj(o_list_sys[0][0], 0, 0, 2*ac, 0);
       }
     l00->names[l00->c_obj] = &l00->a_char[c]; /* add name to end */
     if (lp < l00->c_obj)  /* split index */
       {
        for (j = l00->c_obj; j > lp; j--)
         {
          l00->a_int[j] = l00->a_int[j-1];
         }
       }
     l00->a_int[lp] = l00->c_int;
     l00->p_obj[l00->c_int] = p;
     strcpy(l00->names[l00->c_int], p->key);
     for (k = 1; k < N_SYSTEM; k++)  lsys[k]->a_int[l00->c_int] = 0;
     l00->c_obj = ++l00->c_int; l00->c_char += sl;
     for (k = 1; k < N_SYSTEM; k++)  lsys[k]->c_int = l00->c_int;
    }
#ifdef _TREE_
       puts("++++++ exit o_list_add");
#endif
}

void open_new()
{
  int j, k, shft = 6*N_SYSTEM;
  int sys00_l[4] = {O_LIST_SIZE, 0, 16*O_LIST_SIZE, O_LIST_SIZE};
  int sys01_l[4] = {O_LIST_SIZE, 0, 0, 0};

#ifdef _TREE_
       puts("++++++ enter open_new");
#endif
  o_list_sys[0][0] = make_obj(system_0[0],
         sys00_l[0], sys00_l[1], sys00_l[2], sys00_l[3]);
  for (k = 1; k < N_SYSTEM; k++)
    {
     o_list_sys[0][k] = make_obj(system_0[k],
                        sys01_l[0], sys01_l[1], sys01_l[2], sys01_l[3]);
    }
  for (j = 0; j < INFO_BASE + N_SYSTEM * INFO_STEP; j++)  info_int[j] = 0;
  info_int[5] = DYNAMIC_CODE;
  info_int[shft] = major_version;
  info_int[shft+1] = minor_version;
  info_int[shft+2] = 0; /* current last count for dead records */
}

void open_old()
{
  int j, k, lp, index, index_old;

#ifdef _TREE_
       puts("++++++ enter open_old");
#endif
  db_get_info();
  db_get_system();
  if (o_list_sys[0][10] != NULL)       /* previous list exists */
    {
     reopen = 1;
     for (j = 0; j < o_list_sys[0][0]->c_obj; j++)
      { 
       index = o_list_sys[0][0]->a_int[j];
       lp = 
       n_list_pos(o_list_sys[0][0]->names[index], o_list_sys[0][10]);
       if (lp < 0)  /* object exists in old list */
        {
         index_old = o_list_sys[0][10]->a_int[-lp-1];
         if (o_list_sys[0][10]->p_obj[index_old] != NULL)
                                     /* object is still in memory */
  	{
           if (o_list_sys[0][1]->a_int[index] ==
               o_list_sys[0][11]->a_int[index_old] &&
               o_list_sys[0][2]->a_int[index] ==
               o_list_sys[0][12]->a_int[index_old]) /* not modified */
  	  {
             o_list_sys[0][0]->p_obj[index] = 
             o_list_sys[0][10]->p_obj[index_old];
             prev_kept++;
  	  }
           else
  	  {
             free_obj(o_list_sys[0][10]->p_obj[index_old]);
             prev_freed++;
            }
  	}
        }
      }
     for (k = 0; k < N_SYSTEM; k++) free_obj(o_list_sys[0][10+k]);
    }
  else  reopen = 0;
}

void prt_obj(struct object* p)
{
  printf("object: %s  parent: %s type: %s ma+mi time: %d %d\n", 
         p->key, p->par_name, p->obj_type, p->ma_time, p->mi_time);
  printf("int l/c: %d %d dp l/c: %d %d char l/c: %d %d obj l/c: %d %d\n",
         p->l_int, p->c_int, p->l_dble, p->c_dble,
         p->l_char, p->c_char, p->l_obj, p->c_obj);
}

void prt_obj_full(struct object* p)
{
  int j, k = 0;
  char c[51];

  printf("object: %s  parent: %s base: %s type: %s\n   ma+mi time: %d %d\n", 
         p->key, p->par_name, p->base_name, p->obj_type, 
         p->ma_time, p->mi_time);
  printf("%s%d%s%d\n", "i-length:     ", p->l_int, 
                       "    i-occupation:     ", p->c_int);
  for (j = 0; j < p->c_int; j++)
    {
     if (j-10*(j/10)==0) printf("\n");
     printf("%s%d","  ", p->a_int[j]);
    }
  printf("\n");
  printf("%s%d%s%d\n", "dp-length:    ", p->l_dble,
                       "    dp-occupation:    ", p->c_dble);
  for (j = 0; j < p->c_dble; j++)
    {
     if (j-5*(j/5)==0) printf("\n");
     printf("%s%e","  ", p->a_dble[j]);
    }
  printf("\n");
  printf("%s%d%s%d\n", "char-length:  ", p->l_char,
                       "    char-occupation:  ", p->c_char);
  for (j = 0; j < p->c_char; j++)
    {
     c[k++] = p->a_char[j] == '\0' ? '|' : p->a_char[j];
     if (k == 50)
       {
        c[k] = '\0'; puts(c); k = 0;
       }
    }
  if (k > 0)
    {
    c[k] = '\0'; puts(c);
    }
  printf("%s%d%s%d\n", "object-length:  ", p->l_obj,
                       "    object-occupation:  ", p->c_obj);
}

void prt_i_list(struct object* p)
{
  int j;

  printf("%s%d\n", "integer list occupation: ", p->c_int);
  for (j = 0; j < p->c_int; j++) printf("%d\n", p->a_int[j]);
}

void prt_o_list(struct object* l)
{
  int j;
  if (debug_flag > 1)
    {
     printf("%s%d\n", "object list occupation: ", l->c_obj);
    }
  for (j = 0; j < l->c_obj; j++) puts(l->names[j]);
}

void put_obj(       /* puts an object into db, ignoring object list */
             struct object* p)  /* pointer to object */
{
  int counts[4];

  counts[0] = CONT_OBJ;
  counts[1] = p->c_int;
  counts[2] = p->c_dble;
  counts[3] = p->c_char;

#ifdef _TREE_
       puts("++++++ enter put_obj");
#endif
  if (db_is_open() == 0)
     puts("+++DOOM - put_obj: db not open");
  else
    {
     if (db_put(p->key, counts, &p->ma_time, p->a_int, p->a_dble, 
                p->a_char) == 0)
         printf("%s%s\n", "DOOM - failure to store: ", p->key);
    }
}

void rec_get  /* gets one record from the buffer */
            (int* counts,     /* input: max. length of control, integer,
                                 d.p., char arrays;
                                 output: actual occupations */
             int* a_cont,     /* output: control array */
             int* a_int,      /* output: integer array */
             double* a_dble,  /* output: d.p. array */
             char* a_char)    /* output: string */
{
  char *tmp;
  int *pint, *tint;
  int j, lcounts[4];
  int aux=0;

#ifdef _TREE_
       puts("++++++ enter rec_get");
#endif
  tint = (int*) buffer;
  for (j = 0; j < REC_CONT; j++)
    {
     lcounts[j] = *tint++;
    }
  if (swap_flag)
    {
     swap_int(REC_CONT, lcounts);
     aux = counts[0] == CONT_OBJ ? CONT_INT : counts[0];
    }
  for (j = 0; j < REC_CONT; j++)
    {
     if (counts[j] > lcounts[j])  counts[j] = lcounts[j];
    }
  for (j = 0; j < counts[0]; j++)
    {
     a_cont[j] = *tint++;
    }
  if (swap_flag)  swap_int(aux, a_cont);
  tint += lcounts[0] - counts[0];
   for (j = 0; j < counts[1]; j++)
    {
     a_int[j] = *tint++;
    }
  if (swap_flag) swap_int(counts[1], a_int);
  tint += lcounts[1] - counts[1];
  pint = (int*) a_dble;
  for (j = 0; j < 2 * counts[2]; j++)
    {
     *pint++ = *tint++;
    }
  if (swap_flag) swap_double(counts[2], (int*) a_dble);
  tint += 2 * (lcounts[2] - counts[2]);
  pc = a_char;
  tmp = (char*) tint;
  for (j = 0; j < counts[3]; j++)
    {
     *pc++ = *tmp++;
    }
}

void rec_put /* copies one record -> output buffer */
            (int* counts,     /* input: length of control, integer,
                                 d.p., char arrays */
             int* a_cont,     /* input: control array */
             int* a_int,      /* input: integer array */
             double* a_dble,  /* input: d.p. array */
             char* a_char)    /* input: string */
{
  char *tmp;
  int *tint, *pint, *pkeep;
  int j;
  int aux=0;

#ifdef _TREE_
       puts("++++++ enter rec_put");
#endif
  tint = (int*) buffer;
  for (j = 0; j < REC_CONT; j++)
  	 {
       *tint++ = counts[j];
  	 }
  if (swap_flag)
    {
     swap_int(REC_CONT, (int*) buffer);
     aux = counts[0] == CONT_OBJ ? CONT_INT : counts[0];
    }
  pkeep = tint;
  for (j = 0; j < counts[0]; j++)
  	 {
       *tint++ = a_cont[j];
  	 }
  if (swap_flag) swap_int(aux, pkeep);
  pkeep = tint;
  for (j = 0; j < counts[1]; j++)
  	 {
       *tint++ = a_int[j];
  	 }
  if (swap_flag) swap_int(counts[1], pkeep);
  pint = (int*) a_dble;
  pkeep = tint;
  for (j = 0; j < 2 * counts[2]; j++)
  	 {
       *tint++ = *pint++;
  	 }
  if (swap_flag) swap_double(counts[2], pkeep);
  tmp = (char*) tint;
  pc = a_char;
  for (j = 0; j < counts[3]; j++)
  	 {
       *tmp++ = *pc++;
  	 }
}

void shift_char(  /* shifts character array in itself */
     char* start,
     int diff,    /* offset, > 0: shift up, < 0 : shift down */
     int count)
{
  char *source, *target;

  if (diff > 0)
    {
     source = start + count; target = source + diff;
     while (count > 0)
       {
        *--target = *--source; count--;
       }
    }
  else if(diff < 0)
    {
     source = start; target = source + diff;
     while (count > 0)
       {
        *target++ = *source++; count--;
       }
    }
}

int start_get /* gets one partial record from the database, keeps it */
            (char* rec_key,     /* input: item rec_key */
             int free_buff,    /* input: if not zero, buffer is freed */
             int* counts,     /* input: max. length of control, integer,
                                 d.p., char arrays;
                                 output: actual occupations, 
                                         or -input if record not found */
             int* a_cont,    /* output: control array */
             int* a_int,     /* output: integer array */
             double* a_dble,   /* output: d.p. array */
             char* a_char)   /* output: string */
{
#ifdef _TREE_
       puts("++++++ enter start_get");
#endif
  if (db_is_open() == 0)
    {
      puts("+++DOOM - start_get: db not open");
      return 0;
    }
  else
    {
     if (db_exists(rec_key))
       {
        db_get(rec_key, counts, a_cont, a_int, a_dble, a_char);
        return 1;
       }
     else return 0;
    }
}

void upper_conv(char* str)  /* convert str to upper case */
{
  char* cp = str;
  char c;
  /* leave this as is ! byte swapping screws it up otherwise */
  while (*cp != '\0') {c = *cp; *cp++ = (char) toupper((int) c);}
}

int var_list_pos( /* gives -pos-1 in ordered list, or place to go if absent */
               char* str, /* lookup string */
               struct object* list)
{
  int pos, low, mid, high, last;
  low = 0; last = 0; high = list->c_obj - 1;

#ifdef _TREE_
       puts("++++++ enter var_list_pos");
#endif
  while (low <= high)
    {
     mid = (low + high) / 2;
     if ((pos = strcmp(str, list->names[mid])) < 0)
                       {high = mid - 1; last = mid;}
     else if (pos > 0) {low  = mid + 1; last = low;}
     else              return (-mid-1);
    }
    return last;
}

void error_dump(  /*  dumps all objects starting with a letter */
              int* flag) /* 0: names, 1: partial, 2: full */
{
  int j, index;
  struct object* p=0;

  for (j = 0; j < o_list_sys[0][0]->c_obj; j++)
    {
     index = o_list_sys[0][0]->a_int[j];
     if (o_list_sys[0][0]->p_obj[index] == NULL)
        p = doom_fetch(o_list_sys[0][0]->names[index]);
     if (p != NULL && strcmp(p->obj_type, "FIELD_ERROR") == 0
         && p->c_dble > 0)
       {
              puts("+++ ");
        if (*flag == 0) 
           printf("%d  %s\n", j, p->key);
        else if(*flag == 1)
           prt_obj(p);
        else if(*flag > 1)
           prt_obj_full(p);
       }
    }
}

void swap_int(int n_int, int* t_int)
{
  int temp;
  char *p_in, *p_out;
  int j, k;

  for (j = 0; j < n_int; j++)
    {
     temp = t_int[j];
     p_in = (char*)&temp + 3; p_out = (char*)&t_int[j];
     for (k = 0; k < 4; k++)  *p_out++ = *p_in--;
     p_in = (char*)&temp; p_out = (char*)&t_int[j];
    }
}

void swap_double(int n_double, int* t_double) 
                /* swap int to avoid alignment problems */
{
  int j, temp;

  for (j = 0; j < 2 * n_double; j +=2)
    {
     temp = t_double[j]; t_double[j] = t_double[j+1]; t_double[j+1] = temp; 
     swap_int(1, &t_double[j]); swap_int(1, &t_double[j+1]);
    }
}
