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
 * prf_plural.c conversion: put if not absolute value 1 (one)
 *
 * The long flag is understood.
 * If alternate form is specified, the string to be put
 * follows the conversion character and is terminated by '#',
 * else the string "s" is used.
 * All other parameters are ignored.
 *
 * $Log: prf_plural.c,v $
 * Revision 1.2  2005/07/18 12:33:04  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  88/12/27  09:15:12  heiner
 * Initial revision
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/prf_plural.c,v 1.2 2005/07/18 12:33:04 birke Exp $";
#endif

#include "prf.h"

void
prf_plural(const PrfDest *dp, PrfSpec *psp, va_list *argp) {
    register char   *p;
    register char   *cp = NULL;
    register long int    val;

    p = psp->PSconvp;
    if(psp->PSflags & PRF_ALT) {
        cp = p;
        while(*++p != '#') {
            if(*p == '\0') {
                prf_unknown(dp, psp, argp);
                return;
            }
        }
    }
#if ! PRF_INTISLONG
    if(psp->PSflags & PRF_LONG) {
        val = va_arg(*argp, long int);
    } else
#endif
    {
        val = va_arg(*argp, int);
    }
    if((val != 1) && (val != -1)) {
        if(psp->PSflags & PRF_ALT) {
            ++cp;
            prf__mput(dp, cp, p - cp);
        } else {
            prf__cput(dp, 's');
        }
    }
    psp->PSconvp = p;
}
