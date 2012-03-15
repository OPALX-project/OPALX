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
 * prf.c    interface: print to buffered stream stdout
 *
 * prf( format [, arg ] ... )
 *  char        *format;
 *
 * Print the arguments 'arg' under control of 'format'
 * into the buffered stream "stdout".
 *
 * $Log: prf.c,v $
 * Revision 1.2  2005/07/18 12:32:58  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:58  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.2  89/04/15  18:00:59  heiner
 * Obey PRF_NOSTDIO: no function defined.
 *
 * Revision 1.1  88/12/27  09:15:07  heiner
 * Initial revision
 *
 */
#include "prf.h"

#ifndef PRF_NOSTDIO

#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/prf.c,v 1.2 2005/07/18 12:32:58 birke Exp $";
#endif

void
prf(const char *format, ...) {
    va_list  ap;
    PrfDest  dest;

    va_start(ap, format);
    dest.PDtype = PRFT_STREAM;
    dest.PDflags = 0;
    dest.PDctab = (PrfCtab *)0;
    dest.PDstream = stdout;
    dprfv(&dest, format, &ap);
    va_end(ap);
}
#endif  /* ndef PRF_NOSTDIO */
