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
 * Author:  Heiner Marxen, SPECS GmbH.
 *
 * stdprf.c setup: define standard conversions
 *
 * stdprf()
 *  Sets the "prf" functions to standard values.
 *  Does not touch nonstandard functions.
 *  Note: floating-point conversions are not included.
 *
 * $Log: stdprf.c,v $
 * Revision 1.2  2005/07/18 12:33:06  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:45:00  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  88/12/27  09:15:17  heiner
 * Initial revision
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/stdprf.c,v 1.2 2005/07/18 12:33:06 birke Exp $";
#endif

#include "prf.h"

void
stdprf() {
    stdtprf(&prf_thectab);
}
