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
 * prf_S.c  conversion: put string, C string style
 *
 * Alternate form additionally embeds the string into double quotes.
 * Left adjustment is known.
 * Zero padding is done (this is funny).
 *
 * Octal escape sequences are always given with maximum length.
 * This makes the result unique even when concatenated with arbitrary
 * stuff of the same kind.
 *
 * Special feature for convenience and compatability with other UNIX systems:
 * A NULL pointer given to prf_S is printed as "(null)".
 *
 * $Log: prf_S.c,v $
 * Revision 1.2  2005/07/18 12:32:59  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:58  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.1  89/07/02  21:21:31  heiner
 * Initial revision
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/prf_S.c,v 1.2 2005/07/18 12:32:59 birke Exp $";
#endif

#include "prf.h"

#define OESC_DIGS   3   /* so many octal digits form an octal escape */

#if PRF_ASCII
#  define is_a_print(c)     (((unsigned)((c)-' ')) <= (unsigned)('~'-' '))
#else   /* ! PRF_ASCII */
#  ifndef PRF_NOLIBC
#    include "ctype.h"
/*
 * When your isprint() macro does not require isascii() to be true,
 * the following define may be shortened.
 */
#    define is_a_print(c)   (isascii(c) && isprint(c))
#  else /* def PRF_NOLIBC */
/*
 * Well, we cannot do the wanted job, so ...
 */
#    undef is_a_print
#  endif /* def PRF_NOLIBC */
#endif  /* ! PRF_ASCII */

#ifdef is_a_print
void
prf_S(const PrfDest *dp, PrfSpec *psp, va_list *argp) {
    register char   *p;
    register int     len;
    register char   *str;
    register int     olen;
    register int     c;     /* logically type (char) */
    register int     flags;
    register int     pad;
    char         buf[OESC_DIGS];

    flags = psp->PSflags;
    str = va_arg(*argp, char *);
    if(str == (char *)NULL) {
        str = "(null)";         /* instead of dereferencing nil */
    }
    /*
     * Determine length of source string (obey precision):
     */
    p = str;
    if(flags & PRF_PREC) {
        len = psp->PSprec;      /* so many allowed */
        while((len > 0) && *p) {
            ++p;
            --len;
        }
    } else {
        while(*p) {
            ++p;
        }
    }
    len = p - str;          /* so many we convert */
    /*
     * Scan the source and sum the effective output width:
     */
    for(olen = len, p = str + len ; p > str ;) {
        switch(c = *--p) {
            case '\'':
                if(flags & PRF_ALT) {    /* inside double quotes it ... */
                    break;          /* ... ok to have a single quote */
                }
                /* else we do not know, and are careful/conservative */
            case '\t':
            case '\b':
            case '\r':
            case '\f':
            case '\n':
            case '"':
            case '\\':
                olen += 1;          /* needs a backslash */
                break;
            default:
                if(! is_a_print(c)) {
                    olen += OESC_DIGS;
                }
                break;
        }
    }
    if(flags & PRF_ALT) {
        olen += 2;          /* two enclosing double quotes */
    }
    /*
     * Determine amount of padding:
     */
    if((pad = psp->PSwidth) > olen) {
        pad -= olen;
    } else {
        pad = 0;
    }
    /*
     * Well, it's time to really do the output:
     */
    if(pad && !(flags & PRF_MINUS)) {
        prf__ncput(dp, ((flags & PRF_ZPAD) ? '0' : ' '), pad);
    }
    if(flags & PRF_ALT) {
        prf__cput(dp, '"');
    }
    for(; len ; --len) {
        switch(c = *str++) {
            case '\'':
                if(flags & PRF_ALT) {    /* inside double quotes it ... */
                    prf__cput(dp, c);   /* ... is ok to have single quotes */
                    break;
                }
                goto backsl;
            case '\t':
                c = 't';
                goto backsl;
            case '\b':
                c = 'b';
                goto backsl;
            case '\r':
                c = 'r';
                goto backsl;
            case '\f':
                c = 'f';
                goto backsl;
            case '\n':
                c = 'n';    /* fall through */
            case '"':          /* fall through */
            case '\\':         /* fall through */
backsl:
                prf__cput(dp, '\\');
                prf__cput(dp, c);
                break;
            default:
                if(is_a_print(c)) {
                    prf__cput(dp, c);
                } else {
                    prf__cput(dp, '\\');
                    p = buf + OESC_DIGS
                        - prf__udigs((unsigned long)(unsigned char)c, 8,
                                     "01234567", buf, OESC_DIGS);
                    while(p > buf) {
                        *--p = '0';
                    }
                    prf__mput(dp, buf, OESC_DIGS);
                }
                break;
        }
    }
    if(flags & PRF_ALT) {
        prf__cput(dp, '"');
    }
    if(pad && (flags & PRF_MINUS)) {
        prf__ncput(dp, ((flags & PRF_ZPAD) ? '0' : ' '), pad);
    }
}
#endif  /* def is_a_print */
