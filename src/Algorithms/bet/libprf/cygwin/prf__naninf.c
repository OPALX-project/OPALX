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
 * prf__naninf.c    support: put already fetched nonfinite FP value
 *
 * The precision and the alternate flag are ignored here.
 * Zero padding is done (this is funny).
 *
 * $Log: prf__naninf.c,v $
 * Revision 1.2  2005/07/18 12:33:01  birke
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
    "$Header: /opt/csr/CVS/GUI/libprf/prf__naninf.c,v 1.2 2005/07/18 12:33:01 birke Exp $";
#endif

#include "prfflt.h"

void
prf__naninf(const PrfDest *dp, PrfSpec *psp, double val) {
    register char   *str;
    register int     len;
    register int     width;
    char         padc;
    char         signc;

    /*
     * We put one of the following strings:
     *      (NaN)
     *      (+Inf)
     *      (-Inf)
     */
    str = "NaN";
    signc = '\0';
#ifdef is_NaN
    if(! is_NaN(val)) {
        str = "Inf";
        signc = ((val < 0.0) ? '-' : '+');
    }
#endif
    /*
     * The rest is a simplified and hacked version of "prf_s":
     */
    width = psp->PSwidth;
    padc = ((psp->PSflags & PRF_ZPAD) ? '0' : ' ');
    /*
     * Find out length of source data.  Here is a loop for robustness,
     * it may be replaced by e.g a constant (if appropriate).
     * Later we add 2 or 3 characters to the string.
     */
    {
        register char   *p;
        p = str;
        while(*p) {
            ++p;
        }
        len = 2 + (signc != '\0') + (p - str);
    }
    /*
     * Check for left padding:
     */
    if((width > len) && !(psp->PSflags & PRF_MINUS)) {
        prf__ncput(dp, padc, width - len);
    }
    /*
     * Put the source itself:
     */
    prf__cput(dp, '(' /*)*/);
    if(signc != '\0') {
        prf__cput(dp, signc);
    }
    prf__mput(dp, str, len);
    prf__cput(dp, /*(*/ ')');
    /*
     * Check for right padding:
     */
    if((width > len) && (psp->PSflags & PRF_MINUS)) {
        prf__ncput(dp, padc, width - len);
    }
}
