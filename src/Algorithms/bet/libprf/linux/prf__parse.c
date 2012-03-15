/*
 * Copyright (c) 1988 SPECS GmbH, Berlin
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
 * prf__parse.c support: parse conversion specification
 *
 * char * prf__parse( fmt, psp, argp )
 *
 * Fills the PrfSpec at 'psp' from the format specification at 'fmt',
 * which points behind the percent character.
 * Length arguments may be fetched from varargs list at 'argp'.
 * The pointer to the conversion character is returned.
 *
 * The parsed syntax is
 *  spec:   '%' ['!'] flag* ['0'] [width] [ '.' [prec] ] [lenmod]
 *  flag:   '#' | '-' | '+' | ' '
 *  width:  digit+ | '*'
 *  prec:   digit+ | '*'
 *  lenmod: 'l' | 'h'
 *  digit:  '0' | '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9'
 * Note that the zero padding indicator is NOT part of the width.
 *
 * $Log: prf__parse.c,v $
 * Revision 1.2  2005/07/18 12:33:01  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.2  89/04/15  18:03:29  heiner
 * Obey PRF_NOLIBC: avoid <ctype.h> to check out decimal digits.
 *
 * Revision 1.1  88/12/27  09:15:10  heiner
 * Initial revision
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/prf__parse.c,v 1.2 2005/07/18 12:33:01 birke Exp $";
#endif

#include "prf.h"

#if PRF_EASYDECDIGS
#  define digit_val(c)      ((c) - '0')
#else   /* ! PRF_EASYDECDIGS */
#  ifndef PRF_NOLIBC
#    include <ctype.h>
/*
 * When your isdigit() macro does not require isascii() to be true,
 * the following definition may be shortened.
 */
#    define is_a_digit(c)   (isascii(c) && isdigit(c))
#  endif /* ndef PRF_NOLIBC */
#endif  /* ! PRF_EASYDECDIGS */

#ifndef is_a_digit
#  define is_a_digit(c)     (((unsigned)(digit_val(c))) <= (unsigned)9)
#endif  /* ndef is_a_digit */

#ifndef digit_val
static int
digit_val(c)
char    c;
{
    static char      digs[11] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 0};
    register char   *p;

    p = digs;
    while(*p && (*p != c)) {
        ++p;
    }
    return(p - digs);
}
#endif  /* ndef digit_val */


const char *
prf__parse(const char *fmt, PrfSpec *psp, va_list *argp) {
    register int     c;     /* logically type (char) */

    if(*fmt == '!') {
        register PrfSpec    *inipsp = va_arg(*argp, PrfSpec *);
        *psp = *inipsp;
        ++fmt;
    } else {
        psp->PSflags = 0;
        psp->PSwidth = 0;
        psp->PSprec = 0;
    }
    psp->PSpercp = fmt - 1;
getflag:                /* Collect the flag characters: */
    switch(c = *fmt++) {
        case '#':
            psp->PSflags |= PRF_ALT;
            goto getflag;
        case '-':
            psp->PSflags |= PRF_MINUS;
            goto getflag;
        case '+':
            psp->PSflags |= PRF_PLUS;
            goto getflag;
        case ' ':
            psp->PSflags |= PRF_BLANK;
            goto getflag;
        default:
            break;
    }
    if(c == '0') {
        psp->PSflags |= PRF_ZPAD;
        c = *fmt++;
    }
    if(is_a_digit(c)) {          /* explicit width */
        register int    width;
        width = 0;
        do {
            width *= 10;
            width += digit_val(c);
            c = *fmt++;
        } while(is_a_digit(c));
        psp->PSwidth = width;
        psp->PSflags |= PRF_WIDTH;
    } else if(c == '*') {         /* fetch argument for width */
        c = *fmt++;
        psp->PSwidth = va_arg(*argp, int);
        if(psp->PSwidth < 0) {   /* simulate minus flag */
            psp->PSwidth = - psp->PSwidth;
            psp->PSflags |= PRF_MINUS;
        }
        psp->PSflags |= PRF_WIDTH;
    }
    if(c == '.') {           /* now expect precision: */
        c = *fmt++;
        if(is_a_digit(c)) {          /* explicit precision */
            register int    prec;
            prec = 0;
            do {
                prec *= 10;
                prec += digit_val(c);
                c = *fmt++;
            } while(is_a_digit(c));
            psp->PSprec = prec;
            psp->PSflags |= PRF_PREC;
        } else if(c == '*') {         /* fetch argument for precision */
            c = *fmt++;
            psp->PSprec = va_arg(*argp, int);
            psp->PSflags |= PRF_PREC;
        }
    }
    if(c == 'l') {
        psp->PSflags |= PRF_LONG;
        c = *fmt++;
    } else if(c == 'h') {
        psp->PSflags |= PRF_HALF;
        c = *fmt++;
    }
    psp->PSconvp = --fmt;
    return(fmt);
}
