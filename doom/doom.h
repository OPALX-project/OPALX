/* (some) constants and structure definitions for DOOM */

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

/*
  element definition: (F: Fortran, C: C/C++)
  d.p. array
  word / e_type = 1            = 2                               = 3
 C    F
 0    1    l      [m]          l [m]                             l [m]
 1    2    rhoinv [1/m]        volt [MV]                          deltap/p
 2    3    e1                  ex [MV/m]                          kick
 3    4    e2                  ey [MV/m]                           .
 4    5    h1                  freq [MHz]                          .
 5    6    h2                  lag [2 Pi]                          .
 6    7    tilt                tilt                                .
 7    8    ks                  betrf                              kick
 8    9    hgap [m]            pg {MW]                 C:8-43 F:9-44: rm
 9   10    fint [Tm]           shunt [MOhm/m]          C:44-259 F:45-260: tm
10   11    angle = K_0*l       tfill [micro sec]
11   12    lrad                harmon
12   13    k0 or k0*l (l=0)    xsize (coll.) or xma (beam-beam) or x (mon.)
13   14    k0s or k0s*l        ysize (coll.) or yma (beam_beam) or y (mon.)
14   15    k1 or k1*l (l=0)    sigx
15   16    k1s or k1s*l        sigy
16   17    k2 or k2*l          fractional charge
17   18    k2s  etc.           npart (# particles in opposite beam)

int array: as d.p. array, containing expression flag:
ex_flag = 1   value
ex_flag > 1   expression

name array:as d.p. array, pointers to parameter names if ex_flag > 0
*/

/*
  parameter definition:
  int  array:
    1    exflag            1 if value, > 1 if expression
  d.p. array:
    1    value             (always)
  char array:
         string            expression as read if exflag > 1
*/

struct object
{
  char key[48];          /* d.b. key */
/* The order of the first 11 variables below is FIXED */
  int ma_time,            /* start of control part = 
                             major time at creation or last modification */
      mi_time,            /* minor time at creation or last modification */
      l_int,              /* length of integer array */
      l_dble,             /* length of double array */
      l_char,             /* length of string */
      l_obj,              /* length of object pointer array */
      c_int,              /* occupation of integer array */
      c_dble,             /* occupation of double array */
      c_char,             /* occupation of string */
      c_obj;              /* occupation of object and names pointer array */
  char par_name[24],      /* parent name */
       base_name[24],     /* basic type name (e.g. QUADRUPOLE, DRIFT,..) */
       obj_type[24];      /* object type such as ELEMENT, TWISS_SUMMARY etc. */

  int* a_int;             /* integer array */
  double* a_dble;         /* d.p. array */
  char* a_char;           /* string */
  struct object* parent;   /* pointer to parent object */
  struct object** p_obj;  /* object pointer array */
  char** names;           /* name pointers into a_char */
};

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif
