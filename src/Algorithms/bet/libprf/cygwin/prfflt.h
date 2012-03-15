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
 * prfflt.h
 *
 * Extension of "prf.h" for floating point oriented stuff.
 * This file is inherently machine dependant and thus
 * needs careful examination when being ported to another system.
 *
 * $Log: prfflt.h,v $
 * Revision 1.2  2005/07/18 12:33:05  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  89/07/02  21:19:50  heiner
 * Initial revision
 *
 */

#include "prf.h"            /* for the sloppy user */

/*
 * The following manifests define some limits for floating point conversions.
 * They are machine dependant, i.e. the stated limits must hold.
 * They do so for IEEE double precision (64 bit, 11 bit exponent).
 */
#define PRF_MAXEXP  330 /* maximum absolute decimal exponent value */
#define PRF_MAXFSIG 17  /* max significant digits of double */

/*
 * The following manifests define some layout decisions:
 */
#define PRF_MAXFCVT PRF_MAXFSIG /* max significant digits for %f */
#define PRF_MAXECVT PRF_MAXFSIG /* max significant digits for %e & %g */
#define PRF_MINEEDIGS   2   /* min significant digits for %e exponent */

/*
 * Two macros are to be defined here:
 *  is_NaN(x)       TRUE iff x is not a number
 *  is_nonfinite(x) TRUE iff x is not a finite number,
 *          i.e. a NaN or some infinity
 *
 * IEEE 754 defines that NaN values are unequal to themselves.
 * IEEE 754 defines nonfinite values to have the exponent bits all 1.
 */
#define is_NaN(x)       ((x) != (x))

#ifdef clipper
/*
 * IEEE 754, least significant bits at lowest address
 */
#  define is_nonfinite(x)   (((((short*)&(x))[3]) & 0x7ff0) == 0x7ff0)
#endif  /* def clipper */

#ifdef m68k
/*
 * IEEE 754, least significant bits at highest address
 */
#  define is_nonfinite(x)   (((((short*)&(x))[0]) & 0x7ff0) == 0x7ff0)
#endif  /* def m68k */

#ifndef is_nonfinite            /* approximation while we don't know */
#  define is_nonfinite(x)   is_NaN(x)
#endif
