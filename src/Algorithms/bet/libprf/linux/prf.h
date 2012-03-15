/*
 * Copyright (c) 1988 SPECS GmbH, Berlin
 * Not derived from licensed software.
 *
 * Permission is granted to freely use, copy, modify, and distribute
 * this software, provided that neither the author nor SPECS GmbH is
 * construed to be liable for any results of using this software,
 * alterations are marked as such, and the above copyright and this notice
 * are not modified.
 *
 * Author:  Heiner Marxen, SPECS GmbH.
 *
 * prf.h
 *
 * Interface definitions for ``prf'' library.
 *
 * $Log: prf.h,v $
 * Revision 1.2  2005/07/18 12:32:59  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:58  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.2  89/04/15  17:58:04  heiner
 * New manifest parameters PRF_NO{STDIO,LIBC}.
 *
 * Revision 1.1  88/12/27  09:13:58  heiner
 * Initial revision
 *
 */

#ifndef PRF_h_INCLUDED          /* guard against multiple inclusion */
#define PRF_h_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#ifndef va_arg
#  if ! PRF_VA_UNUSED || defined(c_plusplus) /* except explicitly not wanted */
#    include <stdarg.h>     /* useful to have it always */
#  endif
#endif

#if defined(PRF_NOLIBC) && !defined(PRF_NOSTDIO)
#  define PRF_NOSTDIO           /* PRF_NOLIBC implies PRF_NOSTDIO */
#endif

#if !defined(PRF_NOSTDIO) && !defined(EOF)
#  include <stdio.h>
#endif

#ifndef NULL
#  define NULL  0
#endif

#ifndef c_plusplus

    typedef struct PrfDest      PrfDest;
    typedef struct PrfSpec      PrfSpec;
    typedef struct PrfLimBuf    PrfLimBuf;
    typedef struct PrfIndir     PrfIndir;
    typedef struct PrfCtab      PrfCtab;

    typedef void (*PrfCvtFPtr)();       /* pointer to conversion function */
    typedef int (*PrfFlsFPtr)();    /* pointer to buffer flush function */
    typedef int (*PrfPutFPtr)();    /* pointer to putchar function */
    typedef void (*PrfIcputFPtr)();     /* ptr to fcn like prf__cput */
    typedef void (*PrfIncputFPtr)();    /* ptr to fcn like prf__ncput */
    typedef void (*PrfImputFPtr)();     /* ptr to fcn like prf__mput */

#else   /* def c_plusplus */

    struct      PrfDest;
    struct      PrfSpec;
    struct      PrfLimBuf;
    struct      PrfCtab;

    typedef void (*PrfCvtFPtr)(PrfDest *, PrfSpec *, va_list *);
    /* pointer to conversion function */
    typedef int (*PrfFlsFPtr)(PrfLimBuf *);
    /* pointer to buffer flush function */
    typedef int (*PrfPutFPtr)(int); /* pointer to putchar function */
    typedef void (*PrfIcputFPtr)(PrfDest *, char);          /* prf__cput */
    typedef void (*PrfIncputFPtr)(PrfDest *, char, int);        /* prf__mcput */
    typedef void (*PrfImputFPtr)(PrfDest *, const char *, int);     /* prf__mput */

#endif  /* def c_plusplus */


    /*
     * Types of "prf" destinations found in PDtype:
     */
#define PRFT_NULL   0       /* to waste paper basket, >/dev/null */
#ifndef PRF_NOSTDIO
#  define PRFT_STREAM   1       /* to a stream */
#endif
#define PRFT_BUFF   2       /* to a buffer */
#define PRFT_LIMBUFF    3       /* to a limited buffer */
#define PRFT_USER   4       /* via user defined putchar function */
#define PRFT_INDIR  5       /* via user defined cput/ncput/mput */

    /*
     * Predefined bits found in PDflags:
     */
#define PRFF_ERROR  0x01        /* some error did occur */


    struct PrfLimBuf {      /* describes limited "prf" memory buffer */
        int      PLBsize;       /* current size of free space */
        char    *PLBptr;        /* current free space */
        char    *PLBbase;       /* begin of total buffer */
        PrfFlsFPtr   PLBflush;      /* user routine on overflow */
    };


    struct PrfIndir {
        PrfIcputFPtr     PIcput;    /* effective prf__cput */
        PrfIncputFPtr    PIncput;   /* effective prf__ncput */
        PrfImputFPtr     PImput;    /* effective prf__mput */
        void        *PIsomep;   /* for convenience */
    };


    struct PrfDest {        /* destination description */
        short        PDtype;    /* type of destination */
        unsigned short   PDflags;   /* for convenience */
        PrfCtab     *PDctab;    /* null, or preferred conversions */
        union {
#ifndef PRF_NOSTDIO
            FILE        *PDDstream; /* target stream */
#endif
            char        *PDDbuff;   /* target buffer pointer */
            PrfLimBuf   *PDDlimb;   /* limited target buffer ptr */
            PrfPutFPtr   PDDufuncp; /* putchar user function */
            PrfIndir    *PDDindir;  /* -> effective cput/ncput/mput */
        }       PDdest;
    };

# ifndef PRF_NOSTDIO
#define PDstream    PDdest.PDDstream
# endif
#define PDbuff      PDdest.PDDbuff
#define PDufuncp    PDdest.PDDufuncp
#define PDindir     PDdest.PDDindir
#define PDlimb      PDdest.PDDlimb

#define PDlimbase   PDlimb->PLBbase
#define PDlimbuff   PDlimb->PLBptr
#define PDlimsize   PDlimb->PLBsize
#define PDlimflush  PDlimb->PLBflush


    struct PrfSpec {        /* describes conversion specification */
        char    *PSpercp;   /* points to initiating '%' */
        char    *PSconvp;   /* points to conversion char */
        int      PSwidth;   /* specified width or 0 */
        int      PSprec;    /* specified precision or 0 */
        int      PSflags;   /* flag bits defined below */
    };

    /*
     * Bits found in PSflags and eventually set by "prf__parse":
     */
#define PRF_ALT     0x001   /* whether '#' seen */
#define PRF_MINUS   0x002   /* whether '-' seen */
#define PRF_PLUS    0x004   /* whether '+' seen */
#define PRF_BLANK   0x008   /* whether ' ' seen */
#define PRF_ZPAD    0x010   /* whether '0' started the width */
#define PRF_WIDTH   0x020   /* whether width was specified */
#define PRF_PREC    0x040   /* whether precision was specified */
#define PRF_LONG    0x080   /* whether 'l' preceded conversion char */
#define PRF_HALF    0x100   /* whether 'h' preceded conversion char */

#define PRF_ALL ( PRF_ALT | PRF_MINUS | PRF_PLUS | PRF_BLANK | PRF_ZPAD \
        | PRF_WIDTH | PRF_PREC | PRF_LONG | PRF_HALF )
    /*
     * Further bits in PSflags may be used by conversion functions:
     */
#define PRF_CAPS    0x200   /* user defined, for convenience */

    /*
     * In order define sufficiently large buffers for conversions we need to
     * know the [maximum] number of bits in a (long), or in an (unsigned long):
     *  PRF_LONGBITS    is the maximum of both
     * NOTE: machine dependant, i.e. check `mklocal.c' and "prflocal.h" when porting
     */
#include "Algorithms/bet/libprf/linux/prflocal.h"

#define PRF_MAXODIGS    ((PRF_LONGBITS + 2) / 3)    /* octal limit */
#define PRF_MAXDDIGS    PRF_MAXODIGS            /* decimal limit */
#define PRF_MAXXDIGS    ((PRF_LONGBITS + 3) / 4)    /* hex limit */

    /*
     * Range of legal conversion characters (cf. mapping to conversion functions).
     * Currently zero and negative values disallowed.
     * We recommend to use letters only (upper case and lower case).
     */
#define PRF_CC_MIN  (1)     /* smallest conversion character */
#define PRF_CC_MAX  PRF_MAXCHAR /* largest conversion character */

#define PRF_CC_LEGAL(c) (((unsigned)(c)-PRF_CC_MIN) <= (PRF_CC_MAX-PRF_CC_MIN))
#define PRF_CC_RANGE    (PRF_CC_MAX - PRF_CC_MIN + 1)
#define PRF_CC_INDEX(c) ((unsigned)(c) - PRF_CC_MIN + 1)    /* into PCtab */

    struct PrfCtab {    /* defines a table of conversion functions */
        PrfCvtFPtr  PCtab[1+PRF_CC_RANGE];  /* [0] is default function */
    };


    /*
     * Global variables of "prf" package:
     */

    extern PrfCtab      prf_thectab;    /* standard conversion function table */


    /*
     * Function definitions of "prf" package:
     */

    /*
     * Note: if a user conversion function modifies the format string,
     *  some the the `const' below aren't true.
     */

# ifndef PRF_NOSTDIO
    extern void prf(const char *, ...);         /* to stdout */
    extern void prfv(const char *, va_list *);      /* to stdout */
    extern void fprf(FILE *, const char *, ...);    /* to stream */
    extern void fprfv(FILE *, const char *, va_list *); /* to stream */
# endif
    extern void sprf(char *, const char *, ...);    /* to memory */
    extern void sprfv(char *, const char *, va_list *); /* to memory */
    extern void lprf(PrfLimBuf *, const char *, ...);
    extern void lprfv(PrfLimBuf *, const char *, va_list *);
    extern void uprf(PrfPutFPtr, const char *, ...);
    extern void uprfv(PrfPutFPtr, const char *, va_list *);
    extern void iprf(PrfIndir *, const char *, ...); /* general extension */
    extern void iprfv(PrfIndir *, const char *, va_list *);
    extern void dprf(PrfDest *, const char *, ...); /* general */
    extern void dprfv(PrfDest *, const char *, va_list *);

    /* set a function, returns old one: */
    extern PrfCvtFPtr   setprf(char, PrfCvtFPtr);
    /* multiple "setprf": */
    extern void defprf(const char *, PrfCvtFPtr, ...);
    extern void clrprf();   /* set all functions undefined */
    extern void stdprf();   /* set standard functions */
    extern void extprf();   /* set extended standard functions */
    extern void fltprf();   /* set floating standard functions */
    /* ... with PrfCtab argument: */
    /* set a function, returns old one: */
    extern PrfCvtFPtr   settprf(PrfCtab *, char, PrfCvtFPtr);
    /* multiple "setprf": */
    extern void deftprf(PrfCtab *, const char *, PrfCvtFPtr, ...);
    extern void clrtprf(PrfCtab *); /* set all functions undefined */
    extern void stdtprf(PrfCtab *); /* set standard functions */
    extern void exttprf(PrfCtab *); /* set extended standard functions */
    extern void flttprf(PrfCtab *); /* set floating standard functions */

    extern void prf__cput(const PrfDest *, char);   /* put a char */
    /* put some memory: */
    extern void prf__mput(const PrfDest *, const char *, int);
    /* put a char multiply */
    extern void prf__ncput(const PrfDest *, char, int);

    /* convert unsigned into digits: */
    extern int  prf__udigs(unsigned long val, int base, const char *digs,
                           char *buf, int len);

    /* regular insertions to buffer: */
    extern int  prf__stretch(char *buf, int buflen, int usedlen, int interval,
                             const char *str, int slen);

    /* puts a numeric conversion string: */
    extern void prf__num(const PrfDest *, PrfSpec *, const char *str, int len,
                         const char *pref, int preflen);

    /* kernel of %d conversion: */
    extern void prf__d(const PrfDest *, PrfSpec *, long);

    /* parse a conversion specification: */
    extern const char  *prf__parse(const char *, PrfSpec *, va_list *);
    /* put a parsed conversion: */
    extern void prf__conv(const PrfDest *, PrfSpec *, va_list *);

    /* standard for unknown conversions: */
    extern void prf_unknown(const PrfDest *, PrfSpec *, va_list *);

    extern void prf_d(const PrfDest *, PrfSpec *, va_list *); /* %d */
    extern void prf_u(const PrfDest *, PrfSpec *, va_list *); /* %u */
    extern void prf_o(const PrfDest *, PrfSpec *, va_list *); /* %o */
    extern void prf_x(const PrfDest *, PrfSpec *, va_list *); /* %x */
    extern void prf_X(const PrfDest *, PrfSpec *, va_list *); /* %X */
    extern void prf_c(const PrfDest *, PrfSpec *, va_list *); /* %c */
    extern void prf_s(const PrfDest *, PrfSpec *, va_list *); /* %s */

    /*
     * Extended conversions:
     */
    extern void prf_b(const PrfDest *, PrfSpec *, va_list *); /* %b */
    extern void prf_C(const PrfDest *, PrfSpec *, va_list *); /* %C */
    extern void prf_p(const PrfDest *, PrfSpec *, va_list *); /* %p */
    extern void prf_R(const PrfDest *, PrfSpec *, va_list *); /* %R */
    extern void prf_S(const PrfDest *, PrfSpec *, va_list *); /* %S */
    extern void prf_m(const PrfDest *, PrfSpec *, va_list *); /* %m */
    extern void prf_repeat(const PrfDest *, PrfSpec *, va_list *); /* %@ */
    extern void prf_plural(const PrfDest *, PrfSpec *, va_list *); /* %P */

    /*
     * Floating point conversions:
     */
    /* NaN or Inf: */
    extern void prf__naninf(const PrfDest *, PrfSpec *, double);
    /* similar to [ef]cvt(3) */
    extern void prf__efcvt(double, char *, int, int, int, int *, int *);
    /* kernel function for %e: */
    extern void prf__e(const PrfDest *, PrfSpec *, double val, const char *mant,
                       int expo, int isneg);
    /* kernel function for %f: */
    extern void prf__f(const PrfDest *, PrfSpec *, double val, const char *mant,
                       int expo, int isneg);

    extern void prf_e(const PrfDest *, PrfSpec *, va_list *); /* %e */
    extern void prf_f(const PrfDest *, PrfSpec *, va_list *); /* %f */
    extern void prf_g(const PrfDest *, PrfSpec *, va_list *); /* %g */

    /* flt standard for %E via prf_e: */
    extern void     prf_eE(const PrfDest *, PrfSpec *, va_list *);
    /* flt standard for %G via prf_g: */
    extern void     prf_gG(const PrfDest *, PrfSpec *, va_list *);


#ifdef __cplusplus
}
#endif


#endif  /* ndef PRF_h_INCLUDED */
