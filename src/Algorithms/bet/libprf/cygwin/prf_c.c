/*
 * Copyright (c) 1988 SPECS, GmbH, Berlin
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
 * prf_c.c  conversion: put character
 *
 * The precision is ignored.
 * Zero padding is done (this is funny). In alternate form with dots.
 * Left adjustment is known.
 * Alternate form is to replace nonprintable chars (not spacing 1 column)
 * by a dot ('.').
 *
 * $Log: prf_c.c,v $
 * Revision 1.2  2005/07/18 12:33:02  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.2  89/04/15  18:04:55  heiner
 * Obey PRF_NOLIBC: avoid <ctype.h> to determine printability.
 *
 * Revision 1.1  88/12/27  09:15:11  heiner
 * Initial revision
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/prf_c.c,v 1.2 2005/07/18 12:33:02 birke Exp $";
#endif

#include "prf.h"

#if PRF_ASCII
#  define is_a_print(c)     (((unsigned)((c)-' ')) <= (unsigned)('~'-' '))
#else   /* ! PRF_ASCII */
#  ifndef PRF_NOLIBC
#    include "ctype.h"
/*
 * When your isprint() macro does not require isascii() to be true,
 * the following define may be shortened.
 */
#    define is_a_print(c)   (isascii(c) && isprint(c))
#  else /* def PRF_NOLIBC */
/*
 * Well, we will IGNORE the alternate flag.
 */
#    undef is_a_print
#  endif /* def PRF_NOLIBC */
#endif  /* ! PRF_ASCII */

void
prf_c(const PrfDest *dp, PrfSpec *psp, va_list *argp) {
    long        lval;
    char        val;
    register int    len;
    register int    flags;
    char        padc;

    flags = psp->PSflags;
#if ! PRF_INTISLONG
    if(flags & PRF_LONG) {
        lval = va_arg(*argp, long int);
    } else
#endif
    {
        lval = va_arg(*argp, int);
    }
    val = lval;
#ifdef is_a_print
    if((flags & PRF_ALT) && ! is_a_print(val)) {
        val = '.';
    }
#endif  /* def is_a_print */
    len = psp->PSwidth;
    if(len <= 1) {
        prf__cput(dp, val);
        return;
    }
    --len;      /* for the char to be put */
    padc = ((flags & PRF_ZPAD) ? ((flags & PRF_ALT) ? '.' : '0') : ' ');
    if(flags & PRF_MINUS) {
        prf__cput(dp, val);
        prf__ncput(dp, padc, len);
    } else {
        prf__ncput(dp, padc, len);
        prf__cput(dp, val);
    }
}
