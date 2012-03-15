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
 * prf__f.c support: for floating point simple style
 *
 * The alternate flag forces a decimal point.
 *
 * The value is already fetched, precision checked, the value converted
 * and rounded (e.g. via "prf__efcvt").
 * The conversion also determined exponent and need of a minus sign.
 * Here we do the rest of the operation.
 * This is expected to be used by "%f" and "%g".
 *
 * $Log: prf__f.c,v $
 * Revision 1.2  2005/07/18 12:33:00  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  89/07/02  21:22:51  heiner
 * Initial revision
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/prf__f.c,v 1.2 2005/07/18 12:33:00 birke Exp $";
#endif

#include "prfflt.h"

struct out {        /* describes part of the output */
    int      o_len;     /* so many characters take */
    char     o_chr;     /* here, or, if '\0' ... */
    char    *o_str;     /* ... from here */
};

/*
 * In the buffer at `mant' are (\0-terminated) the significant digits.
 * `expo' has the position of the decimal point relative
 * to the start of the buffer (negative means left of it).
 * `isneg' is nonzero iff we want a minus sign printed.
 * Rounding has been done for the above specified precision.
 */
void
prf__f(const PrfDest *dp, PrfSpec *psp, double val, const char *mant, int expo, int isneg) {
    register char   *p;
    register struct out *op;
    register int     n;
    register int     prec;
    register struct out *eop;
    struct out       od[6];

    op = od;            /* put descriptors here */
    p = mant;           /* scan digits from here */
    prec = psp->PSprec;     /* wanted precision */
    if(prec < 0) {       /* negative precision is forbidden */
        prec = 0;       /* and forced to zero (needed below) */
    }
    /*
     * Check the sign and optionally put a character.
     * Avoid printing -0.
     * An effective zero output is recognized at the first character.
     */
    if(isneg && (expo > -prec) && (*p != '\0') && (*p != '0')) {
        op->o_len = 1;
        op->o_chr = '-';
        ++op;
    } else if(psp->PSflags & PRF_PLUS) {
        op->o_len = 1;
        op->o_chr = '+';
        ++op;
    } else if(psp->PSflags & PRF_BLANK) {
        op->o_len = 1;
        op->o_chr = ' ';
        ++op;
    }
    /*
     * Scan the significant digits.
     * First the digits before the decimal point (at least one):
     */
    n = expo;           /* so many before decimal point (or negative) */
    if(n <= 0) {         /* nothing before decimal point: force a '0' */
        op->o_len = 1;
        op->o_chr = '0';
        ++op;
    } else {         /* there are digits before '.' */
        op->o_str = p;
        while((n > 0) && *p) {
            ++p;
            --n;
        }
        if((op->o_len = p - op->o_str) != 0) {
            op->o_chr = '\0';
            ++op;
        }
        if(n > 0) {          /* decimal point not yet reached */
            op->o_len = n;
            op->o_chr = '0';
            ++op;
            n = 0;
        }
    }
    /*
     * Decide the decimal point (alternate form forces it).
     */
    if((prec > 0) || (psp->PSflags & PRF_ALT)) {
        op->o_len = 1;
        op->o_chr = '.';
        ++op;
    }
    /*
     * While we have not yet reached the significant digits in the buffer
     * and the precision demands output, we need some more '0's.
     */
    if((n < 0) && (prec > 0)) {
        op->o_chr = '0';
        op->o_len = ((-n < prec) ? -n : prec);
        n += op->o_len;     /* so much nearer to significant digits */
        prec -= op->o_len;  /* so much less precision to put */
        ++op;
    }
    /*
     * If we now really reached the significant digits, then append
     * digits while the precision demands and the buffer provides them.
     */
    if(n >= 0) {
        op->o_str = p;
        while((prec > 0) && (*p != '\0')) {
            ++p;
            --prec;     /* one less to put */
        }
        if((op->o_len = p - op->o_str) != 0) {
            op->o_chr = '\0';
            ++op;
        }
    }
    /*
     * `prec' specifies not yet put but wanted precision.
     * Those must be implicit zero digits.  Note that padding
     * to reach the field width is done with blanks, not with zero digits.
     * Sum the width we will put anyhow:
     */
    eop = op;           /* end of out descriptors */
    for(n = prec, op = od ; op < eop ; ++op) {
        n += op->o_len;
    }
    if(psp->PSwidth > n) {
        n = psp->PSwidth - n;   /* so much padding needed */
    } else {
        n = 0;          /* no blank padding needed */
    }
    /*
     * Well, now lets really produce output:
     */
    if(n && !(psp->PSflags & PRF_MINUS)) {       /* right adjust */
        prf__ncput(dp, ' ', n);
    }
    for(op = od ; op < eop ; ++op) {     /* put body of output */
        if(op->o_len) {
            if(op->o_chr != '\0') {
                prf__ncput(dp, op->o_chr, op->o_len);
            } else {
                prf__mput(dp, op->o_str, op->o_len);
            }
        }
    }
    prf__ncput(dp, '0', prec);          /* further demanded precision */
    if(n && (psp->PSflags & PRF_MINUS)) {    /* left adjust */
        prf__ncput(dp, ' ', n);
    }
}
