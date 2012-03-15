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
 * prf_s.c	conversion: put string
 *
 * Zero padding is done (this is funny). In alternate form with dots.
 * Left adjustment is known.
 * Alternate form is to replace nonprintable chars (not spacing 1 column)
 * by a dot ('.').
 *
 * Special feature for convenience and compatability with other UNIX systems:
 * A NULL pointer given to prf_s is printed as "(null)".
 *
 * $Log: prf_s.c,v $
 * Revision 1.2  2005/07/18 12:33:04  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.2  89/04/15  18:08:05  heiner
 * Obey PRF_NOLIBC: avoid <ctype.h> to determine printability,
 * and avoid <string[s].h> and "strlen()".
 * 
 * Revision 1.1  88/12/27  09:15:14  heiner
 * Initial revision
 * 
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/prf_s.c,v 1.2 2005/07/18 12:33:04 birke Exp $";
#endif

#include "prf.h"

#ifndef PRF_NOLIBC
/*
 * As we use `strlen()', we have to include the appropriate header file
 * (although it is common practice to be sloppy here).
 */
#  ifdef BSD
#    include "string.h"				/* BSD location */
#  else
#    ifdef USG
#      include "strings.h"			/* SYS V location */
#    else	/* ndef BSD && ndef USG */
       extern int	strlen();		/* ahem */
#    endif
#  endif
#endif	/* def PRF_NOLIBC */

#if PRF_ASCII
#  define is_a_print(c)		(((unsigned)((c)-' ')) <= (unsigned)('~'-' '))
#else	/* ! PRF_ASCII */
#  ifndef PRF_NOLIBC
#    include "ctype.h"
     /*
      * When your isprint() macro does not require isascii() to be true,
      * the following define may be shortened.
      */
#    define is_a_print(c)	(isascii(c) && isprint(c))
#  else	/* def PRF_NOLIBC */
     /*
      * Well, we will IGNORE the alternate flag.
      */
#    undef is_a_print
#  endif /* def PRF_NOLIBC */
#endif	/* ! PRF_ASCII */

    void
prf_s( const PrfDest *dp, PrfSpec *psp, va_list *argp )
{
    register char	*p;
    register int	 len;
    register char	*str;
    register int	 flags;
    char		 padc;
    int			 width;

    str = va_arg(*argp, char*);
    if( str == (char*)NULL ) {
	str = "(null)";			/* instead of nil pointer dereference */
    }
    width = psp->PSwidth;
    flags = psp->PSflags;
    padc = ((flags & PRF_ZPAD) ? ((flags & PRF_ALT) ? '.' : '0') : ' ');
    /*
     * Find out length of source data:
     */
    if( flags & PRF_PREC ) {
	len = psp->PSprec;
	p = str;
	while( (len > 0) && *p ) {
	    ++p; --len;
	}
	len = p - str;
    }else {			/* no precision specified */
#ifndef PRF_NOLIBC
	len = strlen(str);
#else	/* def PRF_NOLIBC */
	p = str;
	while( *p ) {
	    ++p;
	}
	len = p - str;
#endif	/* def PRF_NOLIBC */
    }
    /*
     * Check for left padding:
     */
    if( (width > len) && !(flags & PRF_MINUS) ) {
	prf__ncput(dp, padc, width-len);
    }
    /*
     * Put the source itself:
     */
    if( len > 0 ) {
#ifdef is_a_print
	if( flags & PRF_ALT ) {
		register int	 i;
	    /*
	     * Scan source for nonprintable characters:
	     */
	    p = str; i = len;
	    do {
		while( is_a_print(*p) ) {
		    ++p;
		    if( --i <= 0 ) break;
		}
		if( p > str ) {
		    prf__mput(dp, str, p-str);
		}
		if( i <= 0 ) break;
		str = p;
		do {
		    ++p;
		    if( --i <= 0 ) break;
		}while( ! is_a_print(*p) );
		prf__ncput(dp, '.', p-str);
		str = p;
	    }while( i > 0 );
	}else
#endif	/* def is_a_print */
	{
	    prf__mput(dp, str, len);
	}
    }
    /*
     * Check for right padding:
     */
    if( (width > len) && (flags & PRF_MINUS) ) {
	prf__ncput(dp, padc, width-len);
    }
}
