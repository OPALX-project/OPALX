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
 * prfv.c   interface: print to buffered stream stdout
 *
 * prfv( format, argp )
 *
 * Print the arguments at 'argp' under control of 'format'
 * into the buffered stream "stdout".
 *
 * $Log: prfv.c,v $
 * Revision 1.2  2005/07/18 12:33:05  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.2  89/04/15  18:01:18  heiner
 * Obey PRF_NOSTDIO: no function defined.
 *
 * Revision 1.1  88/12/27  09:15:15  heiner
 * Initial revision
 *
 */
#include "prf.h"

#ifndef PRF_NOSTDIO

#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/prfv.c,v 1.2 2005/07/18 12:33:05 birke Exp $";
#endif

void
prfv(const char *format, va_list *argp) {
    PrfDest  dest;

    dest.PDtype = PRFT_STREAM;
    dest.PDflags = 0;
    dest.PDctab = (PrfCtab *)0;
    dest.PDstream = stdout;
    dprfv(&dest, format, argp);
}
#endif  /* ndef PRF_NOSTDIO */
