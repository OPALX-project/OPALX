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
 * sprfv.c  interface: print into memory
 *
 * sprfv( bufp, format, argp )
 *
 * Print the arguments at 'argp' under control of 'format'.
 * Put it into memory at 'bufp'.
 * Add a zero byte at the end (but skip back the pointer).
 *
 * $Log: sprfv.c,v $
 * Revision 1.2  2005/07/18 12:33:06  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  88/12/27  09:15:16  heiner
 * Initial revision
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/sprfv.c,v 1.2 2005/07/18 12:33:06 birke Exp $";
#endif

#include "prf.h"

void
sprfv(char *bufp, const char *format, va_list *argp) {
    PrfDest  dest;

    dest.PDtype = PRFT_BUFF;
    dest.PDflags = 0;
    dest.PDctab = (PrfCtab *)0;
    dest.PDbuff = bufp;
    dprfv(&dest, format, argp);
}
