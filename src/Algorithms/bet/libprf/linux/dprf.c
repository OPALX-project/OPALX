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
 * dprf.c	interface: print to some destination
 *
 * dprf( dp, format [, arg ] ... )
 *	PrfDest		*dp;
 *	char		*format;
 *
 * Print the arguments 'arg' under control of 'format'
 * into the destination described at 'dp'.
 *
 * $Log: dprf.c,v $
 * Revision 1.2  2005/07/18 12:32:56  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:58  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  88/12/27  09:15:04  heiner
 * Initial revision
 * 
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/dprf.c,v 1.2 2005/07/18 12:32:56 birke Exp $";
#endif

#include "prf.h"

    void
dprf( PrfDest* dp, const char* format, ... )
{
    va_list	 ap;

    va_start(ap, format);
    dprfv(dp, format, &ap);
    va_end(ap);
}
