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
 * Author:	Heiner Marxen, SPECS GmbH.
 *
 * prf_m.c	extended conversion: put errno, readable string
 *
 * The alternate flag causes the errno to be fetched from the argument list.
 * Else the external variable "errno" is fetched.
 *
 * NOTE: This conversion function uses "%d" and "%s".
 * NOTE: This might be BSD-special (i.e. need be changed for other systems).
 *
 * $Log: prf_m.c,v $
 * Revision 1.2  2005/07/18 12:33:03  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  90/04/12  21:10:57  heiner
 * Initial revision
 * 
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/prf_m.c,v 1.2 2005/07/18 12:33:03 birke Exp $";
#endif

#include "prf.h"
#include <errno.h>
#include <string.h>

    void
prf_m( const PrfDest *dp, PrfSpec *psp, va_list *argp )
{
/*    extern char	*sys_errlist[];			readable error strings */
/*    extern int	 sys_nerr;			 length of "sys_errlist[]" */
    long int	 val;
    PrfSpec	 spec;

    char         ebuf[256];

    spec = *psp;				/* conservative: get a copy */
    if( spec.PSflags & PRF_ALT ) {		/* fetch argument */
#if ! PRF_INTISLONG
	if( spec.PSflags & PRF_LONG ) {
	    val = va_arg(*argp, long int);
	}else
#endif
	{
	    val = va_arg(*argp, int);
	}
	spec.PSflags &= ~PRF_LONG;		/* we did interpret it */
    }else {
	val = errno;
    }
    spec.PSflags &= ~PRF_ALT;			/* we did interpret it */
    if( (val >= 0) && (strerror_r(val,ebuf,256) != EINVAL) ) {	/* normal case */
      /* RB: deprecated:  dprf(dp, "%!s", &spec, sys_errlist[val]);  */
      dprf(dp, "%!s", &spec, strerror(val));
    }else {
	/*
	 * In case of left adjustment and large width, the following
	 * is not always what we want.
	 * But then, this is the rare case ...
	 */
	static char	suffix[] = " (errno)";

	if( spec.PSwidth > 0 ) {
	    if( 0 > (spec.PSwidth -= (sizeof suffix) - 1) ) {
		spec.PSwidth = 0;
	    }
	}
	if( spec.PSprec > 0 ) {
	    if( 0 > (spec.PSprec -= (sizeof suffix) - 1) ) {
		spec.PSprec = 0;
	    }
	}
	dprf(dp, "%!ld%s", &spec, val, suffix);
    }
}
