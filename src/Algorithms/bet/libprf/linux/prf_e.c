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
 * prf_e.c  conversion: put floating point scientific style
 *
 * The major part of the work is done in the digit string converter,
 * and in the support function "prf__e()".
 *
 * $Log: prf_e.c,v $
 * Revision 1.2  2005/07/18 12:33:02  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  89/07/02  21:24:23  heiner
 * Initial revision
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/prf_e.c,v 1.2 2005/07/18 12:33:02 birke Exp $";
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

#define WANT_LIBC_ECVT  0   /* whether we like "ecvt(3)" from libc */

#if defined(PRF_NOLIBC) || !WANT_LIBC_ECVT
#  define USE_LIBC_ECVT 0
#else
#  define USE_LIBC_ECVT 1
extern char  *ecvt();
#endif

void
prf_e(const PrfDest *dp, PrfSpec *psp, va_list *argp) {
    double   val;
    int      decpt;
    int      isneg;
#if USE_LIBC_ECVT
    char    *cvtbuf;        /* ecvt() tells us where it is */
#else
    char     cvtbuf[PRF_MAXECVT+1]; /* mantissa space */
#endif

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
     * Let a library function do the conversion to decimal digits.
     * We want 1 digit before the decimal point and `prec' many behind,
     * but at most PRF_MAXECVT significant digits at all.
     */
#if USE_LIBC_ECVT
    cvtbuf = ecvt(val, min(PRF_MAXECVT, 1 + psp->PSprec), &decpt, &isneg);
#else
    prf__efcvt(val, cvtbuf, arrelems(cvtbuf),
               min(PRF_MAXECVT, 1 + psp->PSprec), 'e', &decpt, &isneg);
#endif
    /*
     * For the rest we use a support function:
     */
    prf__e(dp, psp, val, cvtbuf, decpt, isneg);
}
