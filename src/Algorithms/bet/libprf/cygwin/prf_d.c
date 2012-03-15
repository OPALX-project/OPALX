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
 * prf_d.c	conversion: put signed decimal
 *
 * $Log: prf_d.c,v $
 * Revision 1.2  2005/07/18 12:33:02  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  88/12/27  09:15:11  heiner
 * Initial revision
 * 
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/prf_d.c,v 1.2 2005/07/18 12:33:02 birke Exp $";
#endif

#include "prf.h"

#define USE__D		1	/* whether to use prf__d() */

#define arrelems(a)	(sizeof(a)/sizeof*(a))
#define eobuf(b)	((b) + arrelems(b))

    void
prf_d( const PrfDest *dp, PrfSpec *psp, va_list *argp )
{
    long int	 val;

#if ! PRF_INTISLONG
    if( psp->PSflags & PRF_LONG ) {
	val = va_arg(*argp, long int);
    }else
#endif
    {	val = va_arg(*argp, int);
    }
#if USE__D
    prf__d(dp, psp, val);
#else	/* ! USE__D */
    {
	    unsigned long int	 uval;
	    char		 buf[ PRF_MAXDDIGS*2 ];
	    int			 len;
	    char		*pref;
	    int			 preflen;
	if( val < 0 ) {
	    uval = -val;
	    pref = "-";
	    preflen = 1;
	}else {
	    uval = val;
	    if( psp->PSflags & PRF_PLUS ) {
		pref = "+"; preflen = 1;
	    }else if( psp->PSflags & PRF_BLANK ) {
		pref = " "; preflen = 1;
	    }else {
		pref = ""; preflen = 0;
	    }
	}
	len = prf__udigs(uval, 10, "0123456789", buf, arrelems(buf));
# if ! STRICT
	if( psp->PSflags & PRF_ALT ) {		/* insert commas */
	    len = prf__stretch(buf, arrelems(buf), len, 3, ",", 1);
	}
# endif	/* ! STRICT */
	prf__num(dp, psp, eobuf(buf) - len, len, pref, preflen);
    }
#endif	/* ! USE__D */
}
