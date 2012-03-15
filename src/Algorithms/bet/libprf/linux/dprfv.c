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
 * dprfv.c	interface: print to some destination
 *
 * dprfv( dp, format, argp )
 *
 * Print the arguments at 'argp' under control of 'format'
 * into the destination described at 'dp'.
 * This is the most general interface.
 *
 * Here we handle the special feature of PRFT_BUFF:
 * the null byte terminating the format string is also put to the buffer,
 * and the pointer is skipped back over this terminator.
 *
 * $Log: dprfv.c,v $
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
	"$Header: /opt/csr/CVS/GUI/libprf/dprfv.c,v 1.2 2005/07/18 12:32:56 birke Exp $";
#endif

#include "prf.h"

    void
dprfv( PrfDest* dp, const char* format, va_list* argp )
{
    register char	*fmt;
    register int	 c;			/* logically type (char) */
    register char	*fmt0;
    PrfSpec		 spec;

    fmt = format;
    for(;;) {
	fmt0 = fmt;				/* start point of scanning */
	while( (c = *fmt++) != '%' ) {
	    if( c == '\0' ) {			/* normal end of format */
		goto putrest;
	    }
	}
	/*
	 * Found a percent designator character ...
	 */
	if( *fmt == '%' ) {			/* is "%%", which is special */
	    prf__mput(dp, fmt0, fmt-fmt0);	/* put literal part */
	    ++fmt;				/* skip second percent */
	    continue;
	}
	prf__mput(dp, fmt0, (fmt-fmt0)-1);	/* put literal part */
	fmt = prf__parse(fmt, &spec, argp);
	/*
	 * Now there must be the conversion specification character:
	 */
	if( *fmt == '\0' ) {			/* unexpected end of format */
	    fmt0 = spec.PSpercp;		/* treated as literal output */
	    ++fmt;
	    break;				/* -> putrest */
	}
	prf__conv(dp, &spec, argp);
	fmt = spec.PSconvp + 1;			/* restart scanning behind it */
    }
putrest:;					/* "fmt" is behind the \0 */
    prf__mput(dp, fmt0, --fmt - fmt0);
    if( dp->PDtype == PRFT_BUFF ) {
	*(dp->PDbuff) = '\0';			/* assert terminator */
    }
}
