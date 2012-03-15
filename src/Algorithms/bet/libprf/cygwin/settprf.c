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
 * settprf.c	setup: define a conversion function
 *
 * PrfCvtFPtr settprf( tp, c, funcp )
 *	PrfCtab		*tp;
 *	char		 c;
 *	PrfCvtFPtr	 funcp;
 * Within the conversion table at 'tp' define the conversion function
 * for character 'c' to be at 'funcp'.
 *
 * $Log: settprf.c,v $
 * Revision 1.2  2005/07/18 12:33:06  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/settprf.c,v 1.2 2005/07/18 12:33:06 birke Exp $";
#endif

#include "prf.h"

#define legalchar(c)	PRF_CC_LEGAL(c)
#define chr_func(tp,c)	((tp)->PCtab[ PRF_CC_INDEX(c) ])
#define def_func(tp)	((tp)->PCtab[ 0 ])


    PrfCvtFPtr			/* former conversion function */
settprf( tp, c, funcp )
    PrfCtab	*tp;		/* within this table */
    char	 c;		/* for this conversion character */
    PrfCvtFPtr	 funcp;		/* define this conversion function */
{
    PrfCvtFPtr	 oldfuncp;
    PrfCvtFPtr	*fpp;

    if( ! legalchar(c) ) {
	return( def_func(tp) );	/* FFS: should we return NULL ? */
    }
    fpp = &(chr_func(tp, c));
    oldfuncp = *fpp;
    *fpp = funcp;
    return( oldfuncp );
}
