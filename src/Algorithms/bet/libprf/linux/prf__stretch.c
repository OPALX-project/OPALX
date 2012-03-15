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
 * prf__stretch.c	support: stretch buffer for regular insertions
 *
 * prf__stretch()
 *	A buffer filled at its end is stretched somewhat to the front
 *	in order to insert a specified string before each fixed length
 *	part of the input, except at the very start.
 *	The resulting filled buffer length is returned.
 *	This is especially useful to insert commas before each 3-digit
 *	block in numbers.
 *	When buffer space is not sufficient, some insertions at the front
 *	are silently suppressed.
 *
 * $Log: prf__stretch.c,v $
 * Revision 1.2  2005/07/18 12:33:01  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  90/04/16  18:41:44  heiner
 * Initial revision
 * 
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/prf__stretch.c,v 1.2 2005/07/18 12:33:01 birke Exp $";
#endif

    int
prf__stretch( buf, buflen, usedlen, interval, str, slen )
    char	*buf;		/* buffer to be modified */
    int		 buflen;	/* allocated size */
    int		 usedlen;	/* that many filled at end */
    int		 interval;	/* before that many chars insert */
    char	*str;		/* ... this string */
    int		 slen;		/* ... of this size */
{
    int		 insert;

    if( (slen <= 0) || (usedlen >= buflen) || (interval <= 0) )
	return( usedlen );	/* nothing inserted */
    /*
     * Determine, how often we want and can insert:
     */
    insert = 0;			/* count insertions */
    {
	    register int	have	= buflen - usedlen;
	    register int	ilen	= usedlen;
	while( (ilen > interval) && (have >= slen) ) {
	    ++insert;
	    ilen -= interval;
	    have -= slen;
	}
    }
    if( insert ) {		/* something wanted and possible */
	/*
	 * Scan the source in blocks of "interval" chars,
	 * and copy down until the last block is reached,
	 * which need not be moved.
	 */
	    register char	*dst;
	    register char	*src;
	    register char	*isrc;
	    register char	*p;
	    register int	 n;
	    register char	*ep;
	ep = &(buf[buflen]);		/* behind last char */
	src = ep - usedlen;		/* first input char */
	isrc = ep - insert*interval;	/* next insertion before this */
	insert *= slen;			/* that much total insertion */
	dst = src - insert;		/* first output char */
	do {
	    do {
		*dst++ = *src++;	/* copy source */
	    }while( src < isrc );
	    p = str;
	    n = slen;
	    do {
		*dst++ = *p++;		/* insert "str" */
	    }while( --n );
	}while( (isrc += interval) < ep );
    }
    return( usedlen + insert );		/* that much now in buffer */
}
