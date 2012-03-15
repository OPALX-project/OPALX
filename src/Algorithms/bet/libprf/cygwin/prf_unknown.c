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
 * prf_unknown.c	conversion: put unknown conversion
 *
 * prf_unknown()
 *	Does not fetch and convert a value, but puts the format string itself.
 *
 * $Log: prf_unknown.c,v $
 * Revision 1.2  2005/07/18 12:33:04  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  88/12/27  09:15:14  heiner
 * Initial revision
 * 
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/prf_unknown.c,v 1.2 2005/07/18 12:33:04 birke Exp $";
#endif

#include "prf.h"

    void
prf_unknown( const PrfDest *dp, PrfSpec *psp, va_list *argp )
{
    prf__mput( dp, psp->PSpercp, psp->PSconvp - psp->PSpercp + 1 );
}
