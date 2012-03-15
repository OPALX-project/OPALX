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
 * flttprf.c	setup: define floating point conversions in a conversion table
 *
 * flttprf()
 *	Within the specified table sets the floating point "prf" functions
 *	to standard values.  Does not define other functions.
 *
 * $Log: flttprf.c,v $
 * Revision 1.2  2005/07/18 12:32:57  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:58  birke
 * First rev. checked in at Bessy II
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/flttprf.c,v 1.2 2005/07/18 12:32:57 birke Exp $";
#endif

#include "prfflt.h"

    void
flttprf( tp )
    register PrfCtab	*tp;		/* in this table */
{
    (void) settprf(tp, 'e', prf_e);		/* scientific style */
    (void) settprf(tp, 'E', prf_eE);		/* scientific style */
    (void) settprf(tp, 'f', prf_f);		/* simple style */
    (void) settprf(tp, 'g', prf_g);		/* appropriate style */
    (void) settprf(tp, 'G', prf_gG);		/* appropriate style */
}
