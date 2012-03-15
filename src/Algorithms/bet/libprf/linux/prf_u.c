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
 * prf_u.c	conversion: put unsigned decimal
 *
 * The alternate flag causes commas to be inserted for each 3-digit block.
 *
 * $Log: prf_u.c,v $
 * Revision 1.2  2005/07/18 12:33:04  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  88/12/27  09:15:14  heiner
 * Initial revision
 * 
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/prf_u.c,v 1.2 2005/07/18 12:33:04 birke Exp $";
#endif

#include "prf.h"

#define arrelems(a)	(sizeof(a)/sizeof*(a))
#define eobuf(b)	((b) + arrelems(b))

    void
prf_u( const PrfDest *dp, PrfSpec *psp, va_list *argp )
{
    unsigned long int	 uval;
    char		 buf[ PRF_MAXDDIGS*2 ];
    int			 len;

#if ! PRF_INTISLONG
    if( psp->PSflags & PRF_LONG ) {
	uval = va_arg(*argp, unsigned long int);
    }else
#endif
    {	uval = va_arg(*argp, unsigned int);
    }
    len = prf__udigs(uval, 10, "0123456789", buf, arrelems(buf));
#if ! STRICT
    if( psp->PSflags & PRF_ALT ) {		/* insert commas */
	len = prf__stretch(buf, arrelems(buf), len, 3, ",", 1);
    }
#endif	/* ! STRICT */
    prf__num(dp, psp, eobuf(buf) - len, len, (char*)0, 0);
}
