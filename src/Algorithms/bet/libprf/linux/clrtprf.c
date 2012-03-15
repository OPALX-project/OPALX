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
 * clrtprf.c	setup: undefine all conversions
 *
 * clrtprf()
 *	Within the specified conversion table sets all "prf" functions
 *	to undefined.  This includes the default function entry.
 *
 * $Log: clrtprf.c,v $
 * Revision 1.2  2005/07/18 12:32:56  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:58  birke
 * First rev. checked in at Bessy II
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/clrtprf.c,v 1.2 2005/07/18 12:32:56 birke Exp $";
#endif

#include "prf.h"

#define arrelems(a)	(sizeof(a)/sizeof*(a))

    void
clrtprf( tp )
    register PrfCtab	*tp;
{
    register unsigned	 i;

    for( i=0 ; i<arrelems(tp->PCtab) ; ++i ) {
	tp->PCtab[i] = (PrfCvtFPtr)NULL;
    }
}
