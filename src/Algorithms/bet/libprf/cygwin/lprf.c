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
 * lprf.c	interface: print to limited user buffer
 *
 * lprf( limbp, format [, arg ] ... )
 *	PrfLimBuf	*limbp;
 *	char		*format;
 *
 * Print the arguments 'arg' under control of 'format'.
 *
 * $Log: lprf.c,v $
 * Revision 1.2  2005/07/18 12:32:58  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:58  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  88/12/27  09:15:06  heiner
 * Initial revision
 * 
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/lprf.c,v 1.2 2005/07/18 12:32:58 birke Exp $";
#endif

#include "prf.h"

    void
lprf( PrfLimBuf* limbufp, const char* format, ... )
{
    va_list	 ap;
    PrfDest	 dest;

    va_start(ap, format);
    dest.PDtype = PRFT_LIMBUFF;
    dest.PDflags = 0;
    dest.PDctab = (PrfCtab*)0;
    dest.PDlimb = limbufp;
    dprfv(&dest, format, &ap);
    va_end(ap);
}
