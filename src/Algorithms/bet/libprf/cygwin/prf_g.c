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
 * prf_g.c  conversion: put floating point appropriate style
 *
 * The alternate flag suppresses stripping of trailing zeroes.
 *
 * $Log: prf_g.c,v $
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
    "$Header: /opt/csr/CVS/GUI/libprf/prf_g.c,v 1.2 2005/07/18 12:33:03 birke Exp $";
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
prf_g(const PrfDest *dp, PrfSpec *psp, va_list *argp) {
    register double  val;
    register int     sig;
    register int     n;
    register char   *p;
    int          expo;
    int          isneg;
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
    {
        double      nrval;      /* not a register */
        nrval = val;
        if(is_nonfinite(nrval)) {
            prf__naninf(dp, psp, val);
            return;
        }
    }
#endif  /* def is_nonfinite */
    /*
     * Check for %d style conversion.
     * This is wanted if the FP value really can be represented in a (long).
     * When doing %d style conversion the specified precision
     * is considered meaningless and cleared.
     * The alternate flag is also cleared as prf__d might interpret it.
     * Avoid casting of too large values, as that may raise exceptions.
     */
    if((((double)PRF_MINLONG) <= val) && (val <= (double)PRF_MAXLONG)) {
        long    ival;
        ival = (long)val;
        if(val == (double)ival) {    /* well, do %d like style */
            psp->PSflags &= ~(PRF_PREC | PRF_ALT);
            psp->PSprec = 1;
            prf__d(dp, psp, ival);
            return;
        }
    }
    /*
     * Adjust precision.
     * It specifies the wanted total number of significant digits.
     */
    sig = psp->PSprec;          /* wanted significant digits */
    if(!(psp->PSflags & PRF_PREC) || (sig < 0)) {
        sig = 6;            /* default precision */
        psp->PSflags |= PRF_PREC;
    } else if(psp->PSprec == 0) {     /* zero is forbidden */
        sig = 1;            /* one digit needed */
    }
    /*
     * We now have to decide for either %e or %f.
     * Details of this decision seem to be different on every system.
     * The BSD manual page for %g reads
     *  `... whichever gives full precision in minimum space.'
     * The manual of System V.3 is more specific:
     *  `... style e will be used only is the exponent
     *   resulting from the conversion is less than -4
     *   or greater than the precision.'
     * For the sake of compatability we use the latter.
     */
#if USE_LIBC_ECVT
    cvtbuf = ecvt(val, min(PRF_MAXECVT, sig), &expo, &isneg);
#else
    prf__efcvt(val, cvtbuf, arrelems(cvtbuf), sig, 'e', &expo, &isneg);
#endif
    /*
     * If not hindered by the alternate flag, reduce `sig'
     * to reflect the present significant digits:
     */
    if(!(psp->PSflags & PRF_ALT)) {
        n = 0;
        p = cvtbuf;
        while((n < sig) && *p) {
            ++n;
            ++p;
        }
        while((n > 0) && (*--p == '0')) {
            --n;
        }
        sig = n;
    }
    /*
     * %e style would print exponent (expo-1) with a 1 smaller precision.
     */
    if((expo < -3) || (expo > sig)) {    /* %e style */
        n = 0;
        if(--sig < 0) {          /* precision bias -1 for 1st digit */
            sig = 0;
        }
    } else {             /* %f style */
        n = 1;
        if(expo < sig) {
            sig -= expo;        /* subtract those before dec. point */
        } else {
            sig = 0;
        }
    }
    psp->PSprec = sig;
    (n ? prf__f : prf__e)(dp, psp, val, cvtbuf, expo, isneg);
}
