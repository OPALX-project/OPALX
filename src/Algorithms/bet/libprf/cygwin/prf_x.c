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
 * prf_x.c	conversion: put (unsigned) hexadecimal
 *
 * Note, that the alternate prefix "0x" is suppressed,
 * when the value is 0 and the precision is specified and 0.
 *
 * $Log: prf_x.c,v $
 * Revision 1.2  2005/07/18 12:33:05  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  88/12/27  09:15:15  heiner
 * Initial revision
 * 
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/prf_x.c,v 1.2 2005/07/18 12:33:05 birke Exp $";
#endif

#include "prf.h"

#define arrelems(a)	(sizeof(a)/sizeof*(a))
#define eobuf(b)	((b) + arrelems(b))

    void
prf_x( const PrfDest *dp, PrfSpec *psp, va_list *argp )
{
    unsigned long int	 uval;
    char		 buf[ PRF_MAXXDIGS ];
    int			 len;
    int			 preflen;

#if ! PRF_INTISLONG
    if( psp->PSflags & PRF_LONG ) {
	uval = va_arg(*argp, unsigned long int);
    }else
#endif
    {	uval = va_arg(*argp, unsigned int);
    }
    len = prf__udigs(uval, 16, "0123456789abcdef", buf, arrelems(buf));
    if(   uval
       && (psp->PSflags & PRF_ALT)
       && ( !(psp->PSflags & PRF_PREC) || psp->PSprec)
      ) {
	preflen = 2;
    }else {
	preflen = 0;
    }
    prf__num(dp, psp, eobuf(buf) - len, len, "0x", preflen);
}
