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
 * prf__mput.c	support: put memory block
 *
 * prf__mput(dp, str, len) 
 *	Puts 'len' characters at 'str' to the destination described at 'dp'.
 *
 * $Log: prf__mput.c,v $
 * Revision 1.2  2005/07/18 12:33:00  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.2  89/04/15  18:03:22  heiner
 * Obey PRF_NOSTDIO: mask stream destination.
 * 
 * Revision 1.1  88/12/27  09:15:09  heiner
 * Initial revision
 * 
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/prf__mput.c,v 1.2 2005/07/18 12:33:00 birke Exp $";
#endif

/*#define PRF_VA_UNUSED	1*/
#include "prf.h"

    void
prf__mput( const PrfDest *cdp, const char *str, int length )
{
   register PrfDest *dp = cdp;
    register char	*p;
    register int	 len;

    if( (len = length) <= 0 ) {		/* nothing to put */
	return;
    }
    p = str;
    switch( dp->PDtype ) {
     default:				/* unknown destination type */
	dp->PDflags |= PRFF_ERROR;
	break;

     case PRFT_NULL:
	break;

#ifndef PRF_NOSTDIO
     case PRFT_STREAM:
	{	register FILE	*fp;
	    fp = dp->PDstream;
	    do {
		if( EOF == putc(*p, fp) ) {
		    dp->PDflags |= PRFF_ERROR;
		}
		++p;
	    }while( --len );
	}
	break;
#endif

     case PRFT_BUFF:
	{	register char	*buf;
	    buf = dp->PDbuff;
	    do {
		*buf++ = *p++;
	    }while( --len );
	    dp->PDbuff = buf;
	}
	break;

     case PRFT_USER:
	{	register PrfPutFPtr	funcp;
	    funcp = dp->PDufuncp;
	    if( funcp != (PrfPutFPtr)NULL ) {
		do {
		    dp->PDflags |= (*funcp)(*p++);
		}while( --len );
	    }
	}
	break;

     case PRFT_LIMBUFF:
	{	register char		*buf;
		register int		 siz;
		register PrfLimBuf	*lbp;
	    lbp = dp->PDlimb;
	    siz = lbp->PLBsize;
	    while( (len > siz) && (len > 0) ) {
					/* Fill rest of buffer: */
		if( siz > 0 ) {
		    len -= siz;
		    buf = lbp->PLBptr;
		    do {
			*buf++ = *p++;
		    }while( --siz );
		    lbp->PLBptr = buf;
		    lbp->PLBsize = 0;
		}
					/* Try flush the buffer: */
		if( lbp->PLBflush == (PrfFlsFPtr)NULL ) {
		    dp->PDflags |= PRFF_ERROR;
		    return;		/* missing flush function */
		}
		dp->PDflags |= (*lbp->PLBflush)(lbp);
					/* Check success of buffer flush: */
		if( (siz = lbp->PLBsize) <= 0 ) {
		    dp->PDflags |= PRFF_ERROR;
		    return;		/* nothing freed by flush function */
		}
	    }
					/* The rest fits into the buffer: */
	    if( len > 0 ) {
		lbp->PLBsize -= len;
		buf = lbp->PLBptr;
		do {
		    *buf++ = *p++;
		}while( --len );
		lbp->PLBptr = buf;
	    }
	}
	break;

     case PRFT_INDIR:
	(*dp->PDindir->PImput)(dp, p, len);
	break;
    }
}
