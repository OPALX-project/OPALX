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
 * prf__thectab.c	support: initialization for "prf_thectab"
 *
 * $Log: prf__thectab.c,v $
 * Revision 1.2  2005/07/18 12:33:01  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/prf__thectab.c,v 1.2 2005/07/18 12:33:01 birke Exp $";
#endif

/*#define PRF_VA_UNUSED	1*/
#include "prf.h"

PrfCtab		prf_thectab		/* THE conversion function table */
			={ {(PrfCvtFPtr)0, /*... filled with same */ } }
			;
