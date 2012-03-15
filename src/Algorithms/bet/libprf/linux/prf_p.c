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
 * prf_p.c  conversion: put pointer appropriate style
 *
 * NOTE: this conversion is machine/system dependant.
 *
 * Pointers are put hexadecimal on this machine.
 * We expect them to fit into an (unsigned long).
 * When no width is specified
 * - an appropriate one is set,
 * - when alternate is specified, it is increased, and
 *   +  zero padding is done,
 *   +  when no precision is given, it is set
 *  to see the maximum possible digits,
 *
 * $Log: prf_p.c,v $
 * Revision 1.2  2005/07/18 12:33:03  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  88/12/27  09:15:13  heiner
 * Initial revision
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/prf_p.c,v 1.2 2005/07/18 12:33:03 birke Exp $";
#endif

#include "prf.h"

#define arrelems(a) (sizeof(a)/sizeof*(a))
#define eobuf(b)    ((b) + arrelems(b))

void
prf_p(const PrfDest *dp, PrfSpec *psp, va_list *argp) {
    register int     flags;
    unsigned long int    uval;
    char         buf[ PRF_MAXXDIGS ];
    int          len;

    flags = psp->PSflags;
    if(!(psp->PSflags & PRF_WIDTH)) {
        psp->PSwidth = PRF_MAXXDIGS;
        if(flags & PRF_ALT) {
            if(!(flags & PRF_PREC)) {
                psp->PSprec = psp->PSwidth;
                flags |= PRF_PREC;
            }
            flags |= PRF_ZPAD;
            psp->PSwidth += 2;
        }
        flags |= PRF_WIDTH;
        psp->PSflags = flags;
    }
    /*  prf_x(dp, psp, argp);   */
    uval = (unsigned long int) va_arg(*argp, void *);
    len = prf__udigs(uval, 16, "0123456789abcdef", buf, arrelems(buf));
    prf__num(dp, psp, eobuf(buf) - len, len, "0x", ((flags & PRF_ALT) ? 2 : 0));
}
