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
 * prf__e.c	support: for floating point scientific style
 *
 * The alternate form forces the decimal point.
 * The CAPS flag changes the exponent introducer 'e' into capital 'E'.
 *
 * The value is already fetched, precision checked, the value converted
 * and rounded (e.g. via "prf__efcvt").
 * The conversion also determined exponent and need of a minus sign.
 * Here we do the rest of the operation.
 * This is expected to be used by "%e" and "%g".
 *
 * $Log: prf__e.c,v $
 * Revision 1.2  2005/07/18 12:33:00  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  89/07/02  21:22:50  heiner
 * Initial revision
 * 
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
	"$Header: /opt/csr/CVS/GUI/libprf/prf__e.c,v 1.2 2005/07/18 12:33:00 birke Exp $";
#endif

#include "prfflt.h"

#if PRF_MINEEDIGS > PRF_MAXDDIGS
	>>>> no No NO <<<<
#endif

#define USE_UDIGS	1		/* whether to use "prf__udigs()" */

#define arrelems(a)	(sizeof(a)/sizeof*(a))
#define eobuf(b)	((b) + arrelems(b))

    void
prf__e( const PrfDest *dp, PrfSpec *psp, double val, const char *mant, int expo, int isneg )
{
    register char	*tp;
    register char	*fp;
    register int	 n;
    register int	 k;
    register char	*suff;
    char		 prefix[3];	/* sign, first digit, dec. point */
    char		 suffix[1+1+PRF_MAXDDIGS+1];	/* stuff expo here */

    fp = mant;
    /*
     * Construct the exponent in the suffix buffer (not null terminated).
     * As we put the first significant digit before the decimal point,
     * the exponent to print is one smaller than `expo'.
     */
    if( (fp[0] == '\0') || (fp[0] == '0') ) {	/* val == 0.0 */
	n = expo = 0;			/* expo for effective 0.0 is +00 */
    }else {
	n = --expo;			/* bias -1 for first significant */
	if( n < 0 ) {
	    n = -n;			/* |expo| */
	}
    }
    tp = eobuf(suffix);			/* behind the buffer */
    /*
     * Prepend significant decimal digits of exponent.
     */
#if USE_UDIGS
    tp -= prf__udigs((unsigned long)n, 10, "0123456789",
					suffix, arrelems(suffix));
#else
    /*
     * For efficiency we here may want to do unsigned arithmetic,
     * but on the other hand exponents are not very large.
     */
    while( n ) {
	k = n % 10;			/* 0 <= k <= 9, we hope */
	*--tp = "0123456789"[k];
	n -= k;				/* are there any rounding divisions ? */
	n /= 10;
    }
#endif
    /*
     * Expand to the minimum wanted exponent digits with '0's:
     */
    while( tp > (eobuf(suffix) - PRF_MINEEDIGS) ) {
	*--tp = '0';
    }
    /*
     * Prepend the exponent sign, which is always present and non-blank,
     * and the 'e' or 'E' according to `PRF_CAPS':
     */
    *--tp = ((expo < 0) ? '-' : '+');
    *--tp = ((psp->PSflags & PRF_CAPS) ? 'E' : 'e');
    suff = tp;				/* save start of suffix */
    /*
     * Construct sign, leading digit and decimal point in "prefix":
     */
    tp = prefix;
    if( isneg ) {
	*tp++ = '-';
    }else if( psp->PSflags & PRF_PLUS ) {
	*tp++ = '+';
    }else if( psp->PSflags & PRF_BLANK ) {
	*tp++ = ' ';
    }
    if( *fp ) {
	*tp++ = *fp++;			/* copy first digit */
    }else {
	*tp++ = '0';			/* simulate first digit */
    }
    n = psp->PSprec;			/* wanted precision */
    if( n < 0 ) {			/* should not happen */
	n = 0;				/* force precision nonnegative */
    }
    if( (n > 0) || (psp->PSflags & PRF_ALT) ) {
	*tp++ = '.';			/* needed or forced */
    }
    mant = fp;				/* mantissa behind decimal point */
    /*
     * Scan rest of mantissa, find where we stop to use it:
     */
    while( (n > 0) && (*fp != '\0') ) {
	++fp; --n;			/* use this significant digit */
    }
    /*
     * `n' is the number of trailing '0's still needed for the precision.
     * Determine the amount of blank padding needed for field width:
     */
    k = (tp - prefix) + (fp - mant) + n + (eobuf(suffix) - suff);
    if( psp->PSwidth > k ) {
	k = psp->PSwidth - k;		/* so much blank padding */
    }else {
	k = 0;				/* no more blank padding */
    }
    /*
     * Now we are really ready for output ...
     */
    if( k && ! (psp->PSflags & PRF_MINUS) ) {	/* right adjust */
	prf__ncput(dp, ' ', k);
    }
    prf__mput(dp, prefix, tp-prefix);
    prf__mput(dp, mant, fp-mant);
    prf__ncput(dp, '0', n);
    prf__mput(dp, suff, eobuf(suffix) - suff);
    if( k && (psp->PSflags & PRF_MINUS) ) {	/* left adjust */
	prf__ncput(dp, ' ', k);
    }
}
