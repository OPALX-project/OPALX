/*
 * Copyright (c) 1989 SPECS GmbH, Berlin
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
 * prf_gG.c conversion: put floating point appropriate
 *
 * This just sets the PRF_CAPS flag, and calls "prf_g".
 *
 * $Log: prf_gG.c,v $
 * Revision 1.2  2005/07/18 12:33:03  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  89/07/02  21:24:24  heiner
 * Initial revision
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/prf_gG.c,v 1.2 2005/07/18 12:33:03 birke Exp $";
#endif

#include "prfflt.h"

void
prf_gG(const PrfDest *dp, PrfSpec *psp, va_list *argp) {
    psp->PSflags |= PRF_CAPS;
    prf_g(dp, psp, argp);
}
