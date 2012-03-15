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
 * prf__ncput.c support: put character multiply
 *
 * prf__ncput(dp, c, count)
 *  Puts 'count' characters 'c' to the destination described at 'dp'.
 *
 * $Log: prf__ncput.c,v $
 * Revision 1.2  2005/07/18 12:33:01  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:59  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.2  89/04/15  18:03:27  heiner
 * Obey PRF_NOSTDIO: mask stream destination.
 *
 * Revision 1.1  88/12/27  09:15:09  heiner
 * Initial revision
 *
 */
#if !defined(lint) && !defined(NO_IDENT)
static char rcsid[] =
    "$Header: /opt/csr/CVS/GUI/libprf/prf__ncput.c,v 1.2 2005/07/18 12:33:01 birke Exp $";
#endif

/*#define PRF_VA_UNUSED 1*/
#include "prf.h"

void
prf__ncput(const PrfDest *cdp, char c, int count) {
    PrfDest *dp = cdp;
    if(count <= 0) {             /* nothing to put */
        return;
    }
    switch(dp->PDtype) {
        default:               /* unknown destination type */
            dp->PDflags |= PRFF_ERROR;
            break;

        case PRFT_NULL:
            break;

#ifndef PRF_NOSTDIO
        case PRFT_STREAM: {
            register FILE   *fp;
            fp = dp->PDstream;
            do {
                if(EOF == putc(c, fp)) {
                    dp->PDflags |= PRFF_ERROR;
                }
            } while(--count);
        }
        break;
#endif

        case PRFT_BUFF: {
            register char   *buf;
            buf = dp->PDbuff;
            do {
                *buf++ = c;
            } while(--count);
            dp->PDbuff = buf;
        }
        break;

        case PRFT_USER: {
            register PrfPutFPtr funcp;
            funcp = dp->PDufuncp;
            if(funcp != (PrfPutFPtr)NULL) {
                do {
                    dp->PDflags |= (*funcp)(c);
                } while(--count);
            }
        }
        break;

        case PRFT_LIMBUFF: {
            register char       *buf;
            register int         siz;
            register PrfLimBuf  *lbp;
            lbp = dp->PDlimb;
            siz = lbp->PLBsize;
            while((count > siz) && (count > 0)) {
                /*
                 * Fill rest of buffer:
                 */
                if(siz > 0) {
                    count -= siz;
                    buf = lbp->PLBptr;
                    do {
                        *buf++ = c;
                    } while(--siz);
                    lbp->PLBptr = buf;
                    lbp->PLBsize = 0;
                }
                /*
                 * Try flush the buffer:
                 */
                if(lbp->PLBflush == (PrfFlsFPtr)NULL) {
                    dp->PDflags |= PRFF_ERROR;
                    return;     /* missing flush function */
                }
                dp->PDflags |= (*lbp->PLBflush)(lbp);
                /*
                 * Check success of buffer flush:
                 */
                if((siz = lbp->PLBsize) <= 0) {
                    dp->PDflags |= PRFF_ERROR;
                    return;     /* nothing freed by flush function */
                }
            }
            /*
             * The rest fits into the buffer:
             */
            if(count > 0) {
                lbp->PLBsize -= count;
                buf = lbp->PLBptr;
                do {
                    *buf++ = c;
                } while(--count);
                lbp->PLBptr = buf;
            }
        }
        break;

        case PRFT_INDIR:
            (*dp->PDindir->PIncput)(dp, c, count);
            break;
    }
}
