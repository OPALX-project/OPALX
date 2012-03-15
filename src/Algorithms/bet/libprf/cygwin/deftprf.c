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
 * Author:	Heiner Marxen, SPECS GmbH.
 *
 * deftprf.c	setup: define multiple conversion functions
 *
 * void deftprf( tp, str, func1, ... )
 *	PrfCtab		*tp;
 *	char		*str;
 *	PrfCvtFPtr	 func1;
 * Within the conversion table at 'tp' defines 'func1' to be the
 * conversion function for the first char of 'str', etc.
 * (just multiple call of "settprf")
 *
 * $Log: deftprf.c,v $
 * Revision 1.2  2005/07/18 12:32:56  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:58  birke
 * First rev. checked in at Bessy II
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/deftprf.c,v 1.2 2005/07/18 12:32:56 birke Exp $";
#endif

#include "prf.h"

    void
deftprf( PrfCtab* tp, const char* c, PrfCvtFPtr fn, ... )
{
    register char	*cstr = c;
    va_list		 ap;

    va_start(ap, fn);
    (void) settprf(tp, *cstr, fn);
    for( cstr = va_arg(ap, char*) ; *cstr != '\0' ; ++cstr ) {
    	(void) settprf(tp, *cstr, va_arg(ap, PrfCvtFPtr));
    }
    va_end(ap);
}
