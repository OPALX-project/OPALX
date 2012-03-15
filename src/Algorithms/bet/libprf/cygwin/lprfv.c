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
 * Author:	Heiner Marxen, SPECS GmbH.
 *
 * lprfv.c	interface: print via limited user buffer
 *
 * lprfv( limbp, format, argp )
 *
 * Print the arguments at 'argp' under control of 'format'.
 *
 * $Log: lprfv.c,v $
 * Revision 1.2  2005/07/18 12:32:58  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:58  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  88/12/27  09:15:07  heiner
 * Initial revision
 * 
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/lprfv.c,v 1.2 2005/07/18 12:32:58 birke Exp $";
#endif

#include "prf.h"

    void
lprfv( PrfLimBuf *limbp, const char *format, va_list *argp )
{
    PrfDest	 dest;

    dest.PDtype = PRFT_LIMBUFF;
    dest.PDflags = 0;
    dest.PDctab = (PrfCtab*)0;
    dest.PDlimb = limbp;
    dprfv(&dest, format, argp);
}
