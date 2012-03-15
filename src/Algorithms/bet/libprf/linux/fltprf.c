/*
 * Copyright (c) 1989 SPECS, GmbH, Berlin
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
 * fltprf.c setup: define floating point conversions
 *
 * fltprf()
 *  Sets the floating point "prf" functions to standard values.
 *  Does not define other functions.
 *
 * $Log: fltprf.c,v $
 * Revision 1.2  2005/07/18 12:32:57  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:58  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  89/07/02  21:20:32  heiner
 * Initial revision
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/fltprf.c,v 1.2 2005/07/18 12:32:57 birke Exp $";
#endif

#include "prfflt.h"

void
fltprf() {
    flttprf(&prf_thectab);
}
