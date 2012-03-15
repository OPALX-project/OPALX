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
 * Author:  Heiner Marxen, SPECS GmbH.
 *
 * exttprf.c    setup: define extended conversions
 *
 * exttprf()
 *  Within the specified conversion table sets the extended "prf" functions
 *  to standard values.  Does not define other functions.
 *
 * $Log: exttprf.c,v $
 * Revision 1.2  2005/07/18 12:32:57  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:58  birke
 * First rev. checked in at Bessy II
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/exttprf.c,v 1.2 2005/07/18 12:32:57 birke Exp $";
#endif

//#define PRF_VA_UNUSED 1
#include "prf.h"

void
exttprf(tp)
register PrfCtab    *tp;
{
    (void) settprf(tp, 'b', prf_b);     /* binary conversion */
    (void) settprf(tp, 'R', prf_R);     /* recursive */
    (void) settprf(tp, 'C', prf_C);     /* C character */
    (void) settprf(tp, 'S', prf_S);     /* C string */
    (void) settprf(tp, '@', prf_repeat);    /* format repetition */
    (void) settprf(tp, 'P', prf_plural);    /* plural construction */
    (void) settprf(tp, 'p', prf_p);     /* pointer conversion */
    (void) settprf(tp, 'm', prf_m);     /* errno conversion */
}
