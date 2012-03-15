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
 * stdtprf.c	setup: define standard conversions
 *
 * stdtprf()
 *	Within the specified conversion table sets the "prf" functions
 *	to standard values.  Does not touch nonstandard functions.
 *	Note: floating-point conversions are not included.
 *
 * $Log: stdtprf.c,v $
 * Revision 1.2  2005/07/18 12:33:06  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:45:00  birke
 * First rev. checked in at Bessy II
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/stdtprf.c,v 1.2 2005/07/18 12:33:06 birke Exp $";
#endif

#include "prf.h"

    void
stdtprf( tp )
    register PrfCtab	*tp;
{
    (void) settprf(tp, 'd', prf_d);
    (void) settprf(tp, 'u', prf_u);
    (void) settprf(tp, 'o', prf_o);
    (void) settprf(tp, 'x', prf_x);
    (void) settprf(tp, 'X', prf_x);
    (void) settprf(tp, 'c', prf_c);
    (void) settprf(tp, 's', prf_s);
}
