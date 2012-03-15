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
 * extprf.c setup: define extended conversions
 *
 * extprf()
 *  Sets the extended "prf" functions to standard values.
 *  Does not define other functions.
 *
 * $Log: extprf.c,v $
 * Revision 1.2  2005/07/18 12:32:56  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:58  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  88/12/27  09:15:05  heiner
 * Initial revision
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/extprf.c,v 1.2 2005/07/18 12:32:56 birke Exp $";
#endif

//#define PRF_VA_UNUSED 1
#include "prf.h"

void
extprf() {
    exttprf(&prf_thectab);
}
