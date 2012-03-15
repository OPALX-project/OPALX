/*
 * Copyright (c) 1989 SPECS GmbH, Berlin
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
 * prf_f.c  conversion: put floating point simple style
 *
 * The major part of the work is done in the digit string converter,
 * and in the support function "prf__f()".
 *
 * $Log: prf_f.c,v $
 * Revision 1.2  2005/07/18 12:33:03  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  89/07/02  21:24:24  heiner
 * Initial revision
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/prf_f.c,v 1.2 2005/07/18 12:33:03 birke Exp $";
#endif

#include "prfflt.h"

#ifndef min
#  define min(a,b)  (((a)<(b)) ? (a) : (b))
#endif

#define arrelems(a) (sizeof(a)/sizeof*(a))
#define eobuf(b)    ((b) + arrelems(b))

/*
 * Note, that we should not have any static data here,
 * in order to allow a bit more recursion via signal handlers.
 * Unfortunately, ecvt() and fcvt() return pointers to static data.
 */

#define WANT_LIBC_FCVT  0       /* whether we like "fcvt()" in libc */

#if defined(PRF_NOLIBC) || !WANT_LIBC_FCVT
#  define USE_LIBC_FCVT 0
#else
#  define USE_LIBC_FCVT 1
extern char  *fcvt();
#endif

void
prf_f(const PrfDest *dp, PrfSpec *psp, va_list *argp) {
    double   val;
    int      decpt;
    int      isneg;
#if USE_LIBC_FCVT
    char    *cvtbuf;        /* fcvt() tells us where it is */
#else   /* ! USE_LIBC_FCVT */
    char     cvtbuf[PRF_MAXFCVT+1]; /* space for mantissa */
#endif  /* ! USE_LIBC_FCVT */

    /*
     * Fetch value and check for non-finiteness:
     */
    val = va_arg(*argp, double);
#ifdef is_nonfinite
    if(is_nonfinite(val)) {
        prf__naninf(dp, psp, val);
        return;
    }
#endif  /* def is_nonfinite */
    /*
     * Adjust precision.  Negative is senseless and changed.
     */
    if(!(psp->PSflags & PRF_PREC) || (psp->PSprec < 0)) {
        psp->PSprec = 6;        /* default precision */
        psp->PSflags |= PRF_PREC;   /* made valid */
    }
    /*
     * Let a library function do the conversion to decimal digits:
     */
#if USE_LIBC_FCVT
    cvtbuf = fcvt(val, min(PRF_MAXFCVT, psp->PSprec), &decpt, &isneg);
#else
    prf__efcvt(val, cvtbuf, arrelems(cvtbuf), psp->PSprec, 'f', &decpt, &isneg);
#endif
    /*
     * For the rest we use a support function:
     */
    prf__f(dp, psp, val, cvtbuf, decpt, isneg);
}
