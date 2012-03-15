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
 * prf_b.c  conversion: put binary
 *
 * The alternate form prints characters "-+" instead of "01".
 *
 * $Log: prf_b.c,v $
 * Revision 1.2  2005/07/18 12:33:02  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  88/12/27  09:15:11  heiner
 * Initial revision
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/prf_b.c,v 1.2 2005/07/18 12:33:02 birke Exp $";
#endif

#include "prf.h"

#define arrelems(a) (sizeof(a)/sizeof*(a))
#define eobuf(b)    ((b) + arrelems(b))

void
prf_b(const PrfDest *dp, PrfSpec *psp, va_list *argp) {
    unsigned long int    uval;
    char         buf[ PRF_LONGBITS ];
    int          len;

#if ! PRF_INTISLONG
    if(psp->PSflags & PRF_LONG) {
        uval = va_arg(*argp, unsigned long int);
    } else
#endif
    {
        uval = va_arg(*argp, unsigned int);
    }
    len = prf__udigs(uval, 2, ((psp->PSflags & PRF_ALT) ? "-+" : "01"),
                     buf, arrelems(buf));
    prf__num(dp, psp, eobuf(buf) - len, len, (char *)0, 0);
}
