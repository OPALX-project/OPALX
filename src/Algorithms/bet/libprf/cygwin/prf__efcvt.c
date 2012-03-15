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
 * prf__efcvt.c	support: convert FP value, useful for %f and %e
 *
 * This function performs an operation similar to "fcvt(3)"/"ecvt(3)",
 * although we here pass a buffer and its size instead of returning
 * a pointer to static data.
 *
 * This stuff is useful for those who do not have "fcvt()" and "ecvt()"
 * or do not want to use it.
 * Together with "prf_f.c" and "prf_e.c" it also reflects what we
 * expect "fcvt(3)" and "ecvt(3)" do.
 * Note that there is no guarantee for any compatabilty
 * (but testing these against the library functions is encouraging).
 *
 * $Log: prf__efcvt.c,v $
 * Revision 1.2  2005/07/18 12:33:00  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  89/05/02  13:50:17  heiner
 * Initial revision
 * 
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/prf__efcvt.c,v 1.2 2005/07/18 12:33:00 birke Exp $";
#endif

#include "prflocal.h"

#if PRF_EASYDECDIGS
#  define decdig(i)	('0' + (i))		/* the easy way */
#else
   static char	decdigs[10] = { '0','1','2','3','4','5','6','7','8','9' };
#  define decdig(i)	(decdigs[i])		/* the truely portable way */
#endif


    void
prf__efcvt( val, buf, siz, round, style, expop, isnegp )
    double	 val;		/* convert this value */
    char	*buf;		/* fill digits here */
    int		 siz;		/* that much space available */
    int		 round;		/* round for this significance */
    int		 style;		/* in this style, 'f' or 'e' */
    int		*expop;		/* fill how many digits before decimal point */
    int		*isnegp;	/* fill whether negative number */
{
    register double	 v;
    register char	*p;
    register char	*q;
    register int	 expo;
    int			 isneg;

    p = buf;
				/* Check for negative; force non-negative */
    v = val;
    if( v < 0.0 ) {
	isneg = 1;
	v = -v;
	val = v;
    }else {
	isneg = 0;
    }
    /*
     * In order to produce decimal digits we scale `v'
     * into the range 1.0 <= v < 10.0,
     * and count the divisions by ten in `expo'.
     * This process may be slow, and can be done more efficiently
     * by using larger powers of ten.
     * Unfortunately, this process can produce inaccurate results.
     */
    expo = 0;
    if( v >= 10.0 ) {				/* scale down */
	    register double	pospow;
#if 1 && (PRF_MAXEXP >= 49)
	pospow = 1e49;
	while( v >= pospow ) {
	    v /= pospow;
	    expo += 49;
	}
#endif
#if 1 && (PRF_MAXEXP >= 7)
	pospow = 1e7;
	while( v >= pospow ) {
	    v /= pospow;
	    expo += 7;
	}
#endif
	pospow = 10.0;
	while( v >= pospow ) {
	    v /= pospow;
	    expo += 1;
	}
    }else if( v < 1.0 ) {			/* scale up */
	    register double	pospow;
	    register double	negpow;
	if( v == 0.0 ) {
		goto out;
	}
#if 1 && (PRF_MAXEXP >= 49)
	pospow = 1e49;
	negpow = 10.0 / pospow;
	while( v < negpow ) {
	    v *= pospow;
	    expo -= 49;
	}
#endif
#if 1 && (PRF_MAXEXP >= 7)
	pospow = 1e7;
	negpow = 10.0 / pospow;
	while( v < negpow ) {
	    v *= pospow;
	    expo -= 7;
	}
#endif
	pospow = 10.0;
	negpow = 1.0;
	while( v < negpow ) {
	    v *= pospow;
	    expo -= 1;
	}
    }
    expo += 1;			/* significant digits before decimal point */
    if( siz <= 1 ) {
	goto out;		/* no space for digits */
    }
    /*
     * From style and rounding precision determine index
     * into the buffer where rounding has to occur.
     *  style 'e':	'round' == wanted digits at all
     *  style 'f':	'round' == wanted digits after true dec point
     */
    if( round < 0 ) {
	round = -1;			/* do not round at all */
    }else {
	if( style == 'f' ) {
	    round += expo;		/* index in p by which to round */
	}
	if( round > (siz-1) ) {
	    round = siz-1;		/* round last in buffer */
	}else if( round < 1 ) {
	    if( round == 0 ) {
		if( v >= 5.0 ) {
		    expo += 1;
		    *p = '1';
		}else {
		    *p = '0';
		}
		++p; --siz;
	    }
	    goto out;
	}
    }
    /*
     * Collect digits binary i.e. not yet really characters.
     *
     * Note, that C does not define whether or how a double
     * value is rounded/truncated when it is casted into int.
     * The following code is totally portable (at least I hope so),
     * but can be accelerated by knowledge of your C semantics.
     */
#define KNOWN_ROUNDTONEAREST	0
#define KNOWN_TRUNCTOZERO	0

#undef NEXT
#if KNOWN_TRUNCTOZERO
#    define NEXT	\
		*p = (int)v;		/* trunc is what we want */	\
		v -= (double)*p;	/* subtract the digit */
#else
#  if KNOWN_ROUNDTONEAREST
#    define NEXT	\
		*p = (int)(v - 0.5);	/* simulate trunc by round */	\
		v -= (double)*p;	/* subtract the digit */
#  endif	/* KNOWN_ROUNDTONEAREST */
#endif	/* KNOWN_TRUNCTOZERO */
#ifndef NEXT
#  define NEXT		\
		*p = (int)v;		/* trunc or round to int */	\
		v -= (double)*p;	/* subtract the digit */	\
		if( v < 0.0 ) {		/* was rounded, we need trunc */\
			v += 1.0;					\
			*p -= 1;					\
		}
#endif	/* ndef NEXT */

    /*
     * Currently "v" is divided by 10^-(expo-1).
     * If expo is greater than 1, we can correct rounding errors
     * after (expo-1) constructed digits,
     * provided we do not stop before that many digits anyhow.
     */
#define WANT_ACCURACY	1
#if WANT_ACCURACY
    if( (expo > 1) && (expo <= siz) && ((round < 0) || (round >= (expo-1))) ) {
	    register double	vput;
	    register int	n;
	n = expo - 1;
	vput = 0.0;
	do {
	    NEXT
	    vput *= 10.0;
	    vput += (double)*p;
	    --n;
	    v *= 10.0;			/* scale up for next digit */
	    --siz; --round; ++p;
	}while( (n > 0) && (siz > 1) && (round != 0) );
	if( n == 0 ) {
	    v = val - 10.0*vput;	/* more correct fraction */
	    while( v >= 10.0 ) {	/* correct like round up */
		v -= 10.0;
		q = p;
		while( (q > buf) && ((q[-1] += 1) >= 10) ) {
		    *--q = 0;
		}
		if( q == buf ) {	/* overflow */
		    *q = 1;
		    expo += 1;
		    v /= 10;
		    if( style == 'f' ) {
			++round;
		    }
		}
	    }
	    while( v < 0.0 ) {		/* correct like rounddown */
		v += 10.0;
		q = p;
		while( *--q == 0 ) {
		    *q = 9;
		}
		if( ((*q -= 1) == 0) && (q == buf) ) {
		    *q = 9;
		    expo -= 1;
		    --p; ++siz;
		    if( style != 'f' ) {
			++round;
		    }
		}
	    }
	}
    }
#endif	/* WANT_ACCURACY */
    while( (siz > 1) && (round != 0) && (v != 0.0) ) {
	NEXT
	v *= 10.0;			/* scale up for next digit */
	--siz; --round; ++p;
    }
    if( round == 0 ) {			/* use the next `digit' for rounding */
	if( v >= 5.0 ) {		/* round up */
	    q = p;
	    while( q > buf ) {
		if( (q[-1] += 1) < 10 ) {
		    break;
		}
		*--q = 0;
	    }
	    if( q == buf ) {		/* overflow by rounding */
		*q = 1;
		expo += 1;
	    }
	}
    }
#if 0			/* activate when trailing '0's are always wanted */
     else {				/* not stopped for rounding */
	while( siz > 1 ) {		/* append some '0's */
	    *p++ = 0;
	    --siz;
	}
    }
#endif
    /*
     * Convert the binary values into real digits:
     */
    for( q=buf ; q<p ; ++q ) {
	*q = decdig(*q);
    }
out:	;
    if( siz > 0 ) {
	*p = '\0';
    }
    if( expop != (int*)0 ) {
	*expop = expo;
    }
    if( isnegp != (int*)0 ) {
	*isnegp = isneg;
    }
}
