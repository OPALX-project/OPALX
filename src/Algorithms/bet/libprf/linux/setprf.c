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
 * setprf.c	setup: define a conversion function
 *
 * PrfCvtFPtr setprf( c, funcp )
 *	char		c;
 *	PrfCvtFPtr	funcp;
 * Define the conversion function for character 'c' to be at 'funcp'.
 *
 * $Log: setprf.c,v $
 * Revision 1.2  2005/07/18 12:33:05  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  88/12/27  09:15:16  heiner
 * Initial revision
 * 
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/setprf.c,v 1.2 2005/07/18 12:33:05 birke Exp $";
#endif

#include "prf.h"

#define legalchar(c)	PRF_CC_LEGAL(c)
#define chr_func(tp,c)	((tp)->PCtab[ PRF_CC_INDEX(c) ])
#define def_func(tp)	((tp)->PCtab[ 0 ])


    PrfCvtFPtr			/* former conversion function */
setprf( c, funcp )
    char	c;		/* for this conversion character */
    PrfCvtFPtr	funcp;		/* define this conversion function */
{
    PrfCvtFPtr	 oldfuncp;
    PrfCvtFPtr	*fpp;

    if( ! legalchar(c) ) {
	return( def_func(&prf_thectab) ); /* FFS: should we return NULL ? */
    }
    fpp = &(chr_func(&prf_thectab, c));
    oldfuncp = *fpp;
    *fpp = funcp;
    return( oldfuncp );
}
