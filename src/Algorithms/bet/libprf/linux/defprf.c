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
 * defprf.c setup: define multiple conversion functions
 *
 * void defprf( str, func1, ... )
 *  char        *str;
 *  PrfCvtFPtr   func1;
 * Defines 'func1' to be the conversion function for the first char of 'str',
 * etc. (just multiple call of "setprf")
 *
 * $Log: defprf.c,v $
 * Revision 1.2  2005/07/18 12:32:56  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:58  birke
 * First rev. checked in at Bessy II
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/defprf.c,v 1.2 2005/07/18 12:32:56 birke Exp $";
#endif

#include "prf.h"

void
defprf(const char *fmt, PrfCvtFPtr fn, ...) {
    register char   *cstr;
    va_list      ap;

    va_start(ap, fn);
    (void) setprf(*fmt, fn);
    for(cstr = va_arg(ap, char *) ; *cstr != '\0' ; ++cstr) {
        (void) setprf(*cstr, va_arg(ap, PrfCvtFPtr));
    }
    va_end(ap);
}
