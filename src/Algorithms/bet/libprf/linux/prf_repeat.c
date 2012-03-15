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
 * prf_repeat.c	conversion: repeat part of format string
 *
 * This conversion is a bit tricky and may be considered a demonstration
 * rather than serious. Nevertheless someone might find it useful.
 *
 * Let '@' be the programmed conversion character. Then the format string
 *	"%9@[%s]%@"
 * fetches 9 character pointers and prints those 9 strings enclosed in brackets.
 * There must not be any size or flag specified at the terminator.
 *
 * Note: "%0@%*d%@" would eat two int's while scanning to find the
 *	terminating "%@", although no output would be done.
 * Note: as start and stop conversions communicate via a (static) global
 *	variable, we may have trouble when some such conversion
 *	is not terminated properly (e.g. via longjump or by missing terminator).
 *
 * $Log: prf_repeat.c,v $
 * Revision 1.2  2005/07/18 12:33:04  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.2  88/12/28  09:39:39  heiner
 * 	Added note for repetition count 0.
 * 
 * Revision 1.1  88/12/27  09:15:13  heiner
 * Initial revision
 * 
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/prf_repeat.c,v 1.2 2005/07/18 12:33:04 birke Exp $";
#endif

#include "prf.h"

static char	**fmt_continue = (char**)NULL;

    void
prf_repeat( const PrfDest *dp, PrfSpec *psp, va_list *argp )
{
    char	 *fmt;
    char	**sav_continue;
    char	 *newfmt;
    int		  count;
    PrfDest	  dest;

    if( (psp->PSflags & PRF_ALL) == 0 ) {
	if( fmt_continue != (char**)NULL ) {
	    /*
	     * Tell our counting incarnation where we stopped:
	     */
	    *fmt_continue = psp->PSconvp;
	    /*
	     * Tell the caller to stop here:
	     */
	    psp->PSconvp = "\0";
	    return;
	}
    }else if( psp->PSflags & PRF_WIDTH ) {
	/*
	 * Set up a repetition. Fetch the count:
	 */
	count = psp->PSwidth;
	if( count <= 0 ) {
	    /*
	     * Do not repeat by faking a temporary NULL destination:
	     */
	    count = 1;
	    dest.PDtype = PRFT_NULL;
	    dest.PDflags = 0;
	    dest.PDctab = dp->PDctab;
	    dp = &dest;
	}
	/*
	 * Save the old "fmt_continue" and set up a new one:
	 */
	newfmt = (char*)NULL;
	sav_continue = fmt_continue;
	fmt_continue = &newfmt;
	/*
	 * Now call "dprfv" repeatedly:
	 */
	fmt = psp->PSconvp + 1;
	do {
	    dprfv(dp, fmt, argp);
	}while( --count );
	/*
	 * Restore old "fmt_continue"
	 * and tell the caller where to continue:
	 */
	fmt_continue = sav_continue;
	psp->PSconvp = ((newfmt != (char*)NULL) ? newfmt : "\0");
	return;
    }
    prf_unknown(dp, psp, argp);
}
