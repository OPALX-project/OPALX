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
 * prf__udigs.c	support: convert unsigned into digits
 *
 * prf__udigs()
 *
 * Converts the unsigned 'val' into a string of digits to the 'base'.
 * The 'base' characters are fetched from 'digs', and filled into
 * the 'len' chars at 'buf', right justified and NOT zero-byte terminated.
 * For a zero no digit is put into the buffer.
 * The number of filled characters is returned.
 * When the buffer space is not sufficient high digits are silently suppressed.
 *
 * $Log: prf__udigs.c,v $
 * Revision 1.2  2005/07/18 12:33:02  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  88/12/27  09:15:10  heiner
 * Initial revision
 * 
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/prf__udigs.c,v 1.2 2005/07/18 12:33:02 birke Exp $";
#endif

#define is2pow(X)	(((X) & ((X) - 1)) == 0)

    int
prf__udigs( val, base, digs, buf, len )
    register unsigned long	 val;
    int				 base;
    register char		*digs;
    char			*buf;
    register int		 len;
{
    register char	*p;
    char		*ep;

    if( (val == 0) || (base <= 0) || (len <= 0) ) {
	return( 0 );		/* nothing put */
    }
    ep = p = &(buf[len]);
    if( base == 1 ) {
	/*
	 * Base 1 (unary) is very special.
	 * This code is here just to be complete.
	 */
	do {
	    *--p = digs[0];
	    if( --len <= 0 ) {
		return( ep - p );
	    }
	}while( --val );
    }else if( is2pow(base) ) {
	    register unsigned int	shift;
	    register unsigned long	mask;
	/*
	 * Powers of two can be handled by shifts and masks.
	 */
	switch( base ) {
	 case 2:	shift = 1; mask = 0x01; break;
	 case 4:	shift = 2; mask = 0x03; break;
	 case 8:	shift = 3; mask = 0x07; break;
	 case 16:	shift = 4; mask = 0x0f; break;
	 default:
	    --base;
	    mask = base;
	    shift = 0;
	    do {
		++shift;
	    }while( mask >>= 1 );
	    mask = base;
	    break;
	}
	do {
	    *--p = digs[val & mask];
	    if( --len <= 0 ) {
		return( ep - p );
	    }
	}while( val >>= shift );
    }else {
	/*
	 * Base is not known to be a power of two.
	 * We have to do divisions:
	 */
	    register unsigned	ubase;
	    register unsigned	dig;
	ubase = base;
	while( val >= ubase ) {
	    dig = val % ubase;
	    *--p = digs[ dig ];
	    if( --len <= 0 ) {
		return( ep - p );
	    }
	    val -= dig;		/* are there rounding divisions ? */
	    val /= ubase;
	}
	*--p = digs[val];
    }
    return( ep - p );
}
