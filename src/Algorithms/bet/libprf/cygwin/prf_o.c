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
 * prf_o.c	conversion: put octal
 *
 * $Log: prf_o.c,v $
 * Revision 1.2  2005/07/18 12:33:03  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  88/12/27  09:15:12  heiner
 * Initial revision
 * 
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/prf_o.c,v 1.2 2005/07/18 12:33:03 birke Exp $";
#endif

#include "prf.h"

#define arrelems(a)	(sizeof(a)/sizeof*(a))
#define eobuf(b)	((b) + arrelems(b))

    void
prf_o( const PrfDest *dp, PrfSpec *psp, va_list *argp )
{
    unsigned long int	 uval;
    char		 buf[ 1 + PRF_MAXODIGS ];
    char		*p;
    int			 len;

#if ! PRF_INTISLONG
    if( psp->PSflags & PRF_LONG ) {
	uval = va_arg(*argp, unsigned long int);
    }else
#endif
    {	uval = va_arg(*argp, unsigned int);
    }
    len = prf__udigs(uval, 8, "01234567", buf, arrelems(buf));
    p = eobuf(buf) - len;
    if(   (psp->PSflags & PRF_ALT) 
       && (uval  ||  !(psp->PSflags & PRF_PREC)  ||  psp->PSprec)
      ) {
	*--p = '0'; ++len;
    }
    prf__num(dp, psp, p, len, (char*)0, 0);
}
