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
 * prf__conv.c	support: call the correct conversion function
 *
 * void prf__conv( dp, psp, argp )
 *	PrfDest	*dp;
 *	PrfSpec	*psp;
 *	va_list	*argp;
 * Call the conversion routine for the conversion described at 'psp'
 * with the above 3 arguments.
 *
 * $Log: prf__conv.c,v $
 * Revision 1.2  2005/07/18 12:32:59  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:58  birke
 * First rev. checked in at Bessy II
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/prf__conv.c,v 1.2 2005/07/18 12:32:59 birke Exp $";
#endif

#include "prf.h"

#define legalchar(c)	PRF_CC_LEGAL(c)
#define chr_func(tp,c)	((tp)->PCtab[ PRF_CC_INDEX(c) ])
#define def_func(tp)	((tp)->PCtab[ 0 ])

    void
prf__conv( const PrfDest *dp, PrfSpec *psp, va_list *argp )
{
    register PrfCvtFPtr	 funcp;
    register char	 c;
    register PrfCtab	*tp;

    c = *(psp->PSconvp);
    if(   ((tp = dp->PDctab) == (PrfCtab*)NULL)
       || (   (   ! legalchar(c)
               || ((funcp = chr_func(tp, c)) == (PrfCvtFPtr)NULL)
	      )
	   && ((funcp = def_func(tp)) == (PrfCvtFPtr)NULL)
	  )
      ) {
	tp = &prf_thectab;			/* use default table: */
	if(   (   ! legalchar(c)
	       || ((funcp = chr_func(tp, c)) == (PrfCvtFPtr)NULL)
	      )
	   && ((funcp = def_func(tp)) == (PrfCvtFPtr)NULL)
	  ) {
	    funcp = prf_unknown;		/* the last resort */
	}
    }
    (*funcp)(dp, psp, argp);
}
