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
 * prf__num.c	support: put converted numeric string
 *
 * prf__num() 
 *
 * Puts an already converted number at 'str' ('len' characters)
 * to the destination described at 'dp'.
 * The conversion specification is at 'psp'.
 * Padding for width and precision is handled.
 *
 * Left padding with zero digits occurs to reach the precision.
 * The default precision here is 1.
 *
 * If the number is to be right justified and zero padding is specified,
 * first the header string (the 'preflen' characters at 'pref') is put,
 * then padding '0' characters, then the value.
 * Else the header is prepended to the value,
 * and left or right padding with blanks is done to reach the specified width.
 *
 * The header typically contains signs or e.g. "0x".
 *
 * $Log: prf__num.c,v $
 * Revision 1.2  2005/07/18 12:33:01  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  88/12/27  09:15:09  heiner
 * Initial revision
 * 
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/prf__num.c,v 1.2 2005/07/18 12:33:01 birke Exp $";
#endif

/*#define PRF_VA_UNUSED	1*/
#include "prf.h"

    void
prf__num( const PrfDest *dp, PrfSpec *psp, const char *str, int len, const char *pref, int preflen )
{
    register int	 width;
    register int	 prec;
    register int	 nz;

    if( psp->PSflags & PRF_PREC ) {
	prec = psp->PSprec;
    }else {
	prec = 1;
    }
    if( prec < len ) {
	prec = len;		/* pull up precision to forced length */
    }
    width = psp->PSwidth;
    if( width < prec ) {
	width = prec;		/* pull up width to forced length */
    }
    /*
     * Now we have:	width >= prec >= len
     */
    if( psp->PSflags & PRF_MINUS ) {
	/*
	 * Left justified. No zero padding for width.
	 */
	if( preflen > 0 ) {
	    prf__mput(dp, pref, preflen);
	    width -= preflen;
	}
	if( prec > len ) {
	    prf__ncput(dp, '0', prec-len);
	}
	if( len > 0 ) {
	    prf__mput(dp, str, len);
	}
	if( width > prec ) {
	    prf__ncput(dp, ' ', width-prec);
	}
    }else {
	/*
	 * Right justified. Zero padding for width may occur.
	 */
	nz = 0;			/* in case of no zero padding */
	if( preflen > 0 ) {
	    if( psp->PSflags & PRF_ZPAD ) {
		prf__mput(dp, pref, preflen);
		if( width > (len+preflen) ) {
		    nz = width - (len+preflen);
		}
	    }else {
		if( width > (prec+preflen) ) {
		    prf__ncput(dp, ' ', width - (prec+preflen));
		}
		prf__mput(dp, pref, preflen);
		if( prec > len ) {
		    nz = prec - len;
		}
	    }
	}else {			/* (preflen <= 0):  no prefix present */
	    if( psp->PSflags & PRF_ZPAD ) {
		if( width > len ) {
		    nz = width - len;
		}
	    }else {
		if( width > prec ) {
		    prf__ncput(dp, ' ', width-prec);
		}
		if( prec > len ) {
		    nz = prec - len;
		}
	    }
	}
	if( nz > 0 ) {
	    prf__ncput(dp, '0', nz);
	}
	if( len > 0 ) {
	    prf__mput(dp, str, len);
	}
    }
}
