/*
 * Copyright (c) 1990 SPECS, GmbH, Berlin
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
 * iprf.c   interface: print via indirect functions
 *
 * iprf( indirp, format [, arg ] ... )
 *  PrfIndir    *indirp;
 *  char        *format;
 *
 * Print the arguments 'arg' under control of 'format',
 * using the functions specified at 'indirp'.
 *
 * $Log: iprf.c,v $
 * Revision 1.2  2005/07/18 12:32:57  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:58  birke
 * First rev. checked in at Bessy II
 *
 */
#include "prf.h"

#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/iprf.c,v 1.2 2005/07/18 12:32:57 birke Exp $";
#endif

void
iprf(PrfIndir *indirp, const char *format, ...) {
    va_list  ap;
    PrfDest  dest;

    va_start(ap, format);
    dest.PDtype = PRFT_INDIR;
    dest.PDflags = 0;
    dest.PDctab = (PrfCtab *)0;
    dest.PDindir = indirp;
    dprfv(&dest, format, &ap);
    va_end(ap);
}
