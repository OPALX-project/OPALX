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
 * prf_R.c	conversion: put recursively
 *
 * Recursively calls "dprfv" with a format string fetched from the arguments.
 * The alternate form also fetches a "va_list*" argument
 * for the embedded call of "dprfv".
 * All other flags and sizes are ignored.
 *
 * $Log: prf_R.c,v $
 * Revision 1.2  2005/07/18 12:32:59  birke
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
	"$Header: /opt/csr/CVS/GUI/libprf/prf_R.c,v 1.2 2005/07/18 12:32:59 birke Exp $";
#endif

#include "prf.h"

    void
prf_R( const PrfDest *dp, PrfSpec *psp, va_list *argp )
{
    char	*format;

    format = va_arg(*argp, char *);
    if( psp->PSflags & PRF_ALT ) {
	argp = va_arg(*argp, va_list*);
    }
    dprfv(dp, format, argp);
}
