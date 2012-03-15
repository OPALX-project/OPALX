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
 * Author:  Heiner Marxen, SPECS GmbH.
 *
 * prf__cput.c  support: put character
 *
 * prf__cput(dp, c)
 *  Puts a character 'c' to the destination described at 'dp'.
 *
 * $Log: prf__cput.c,v $
 * Revision 1.2  2005/07/18 12:33:00  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:58  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.2  89/04/15  18:02:21  heiner
 * Obey PRF_NOSTDIO: mask stream destination.
 *
 * Revision 1.1  88/12/27  09:15:08  heiner
 * Initial revision
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/prf__cput.c,v 1.2 2005/07/18 12:33:00 birke Exp $";
#endif

/*#define PRF_VA_UNUSED 1*/
#include "prf.h"

void
prf__cput(const PrfDest *cdp, char c) {
    PrfDest *dp = cdp;
    switch(dp->PDtype) {
        default:               /* unknown destination type */
            dp->PDflags |= PRFF_ERROR;
            break;

        case PRFT_NULL:
            break;

#ifndef PRF_NOSTDIO
        case PRFT_STREAM:
            if(EOF == putc(c, dp->PDstream)) {
                dp->PDflags |= PRFF_ERROR;
            }
            break;
#endif

        case PRFT_BUFF:
            *(dp->PDbuff)++ = c;
            break;

        case PRFT_LIMBUFF:
            if(dp->PDlimsize <= 0) {
                if(dp->PDlimflush != (PrfFlsFPtr)NULL) {
                    dp->PDflags |= (*dp->PDlimflush)(dp->PDlimb);
                }
            }
            if(dp->PDlimsize > 0) {
                *(dp->PDlimbuff)++ = c;
                --(dp->PDlimsize);
            } else {
                dp->PDflags |= PRFF_ERROR;
            }
            break;

        case PRFT_USER:
            if(dp->PDufuncp != (PrfPutFPtr)NULL) {
                dp->PDflags |= (*dp->PDufuncp)(c);
            }
            break;

        case PRFT_INDIR:
            (*dp->PDindir->PIcput)(dp, c);
            break;
    }
}
