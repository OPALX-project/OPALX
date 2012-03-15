/*
 * (Incomplete) example demonstrating the use of "prf".
 * $Header: /opt/csr/CVS/GUI/libprf/Example.c,v 1.2 2005/07/18 12:32:55 birke Exp $
 */

#include <stdlib.h>
#include <prf.h>

static char     *myname;
static int       dbglev;

static void /*ARGSUSED*/
prf_I(PrfDest *dp, PrfSpec *psp, va_list *argp) {
    dprf(dp, "%S", myname);
}

static void
usage(int xcode) {
    prf("Usage: %I [-d?] [number]\n");
    exit(xcode);
}

void    /*VARARGS*/
error(char *fmt, ...) {
    va_list  argp;

    va_start(argp, fmt);
    fprf(stderr, "(%I:) %#R\n", fmt, &argp);
    va_end(argp);
    (void) exit(1);
}

void    /*VARARGS*/
msg(char *fmt, ...) {
    va_list  argp;

    va_start(argp, fmt);
    fprfv(stderr, fmt, &argp);
    va_end(argp);
    (void) fflush(stderr);
}

static void
work(int arg) {
    if(dbglev >= 2) msg("%I: entering work(%d)\n", arg);

    /* ... whatsoever ... */

    if(arg > 0) {
        work(arg >> 1);
    }
    if(dbglev >= 2) msg("%I: leaving work(%d)\n", arg);
}


main(int argc, char **argv) {
    int       ac    = 0;
    char     *cp;
    int       arg   = 1;

    myname = argv[ac++];        /* what we are called today */
    stdprf();
    extprf();
    //    setprf('I', prf_I);
    while((ac < argc) && (argv[ac][0] == '-')) {
        cp = argv[ac++];
        while(*++cp) switch(*cp) {
                case 'd':
                    dbglev++;
                    break;
                case '?':
                    usage(0);
                default:
                    error("unknown option -%C, try -?", *cp);
            }
    }
    if(ac < argc) {
        arg = atoi(argv[ac++]);
        if((arg <= 0) || (arg > 99)) {
            error("bad numeric argument %#S", argv[ac-1]);
        }
    }
    if(ac < argc) usage(1);
    if(dbglev) msg("%I: going to work ...\n");
    work(arg);
    if(dbglev) msg("%I: ready\n");
    exit(0);
}
