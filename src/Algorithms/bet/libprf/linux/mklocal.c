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
 * mklocal.c	determine and print manifests that are machine dependant
 *
 * The result of running this program is expected to be put in "prflocal.h",
 * which is included in "prf.h".
 *
 * Currently the following manifests are determined and printed:
 *  PRF_MININT		most negative value of type (int)
 *  PRF_MAXINT		most positive value of type (int)
 *  PRF_MINCHAR		most negative value of type (char)
 *  PRF_MAXCHAR		most positive value of type (char)
 *  PRF_MINLONG		most negative value of type (long)
 *  PRF_MAXLONG		most positive value of type (long)
 *  PRF_LONGBITS	number of bits which is sufficient to code
 *			all values in the range [PRF_MINLONG..PRF_MAXLONG]
 *			and all values of type (unsigned long).
 *  PRF_INTISLONG	Whether the types (int) and (long) are the same.
 *  PRF_EASYDECDIGS	Whether the codes of '0'..'9' are linear and dense.
 *  PRF_ASCII		Whether our native character set is ASCII,
 *			at least in so far as all characters that
 *			are printable in ASCII have the expected codes.
 *			For suppression of this see below.
 *  PRF_NIZ_CFP		Whether a conversion function pointer filled
 *			will zero bits is a NULL pointer of that type.
 *
 * Those of the above manifests that code a "whether ..." are defined
 * as (1) for TRUE and are #undef'ed (!) for FALSE.
 * This allows to use both `#if' and `#ifdef' to check the condition.
 *
 * There are some more values computed here, for which manifest definitions
 * could be printed.  So if you like you can extend the output ...
 *
 * Some parts of this program reflect paranoia, which is not a bug.
 *
 * $Log: mklocal.c,v $
 * Revision 1.2  2005/07/18 12:32:58  birke
 * before start of Linux-migration
 *
 * Revision 1.1.1.1  1995/06/09  08:44:58  birke
 * First rev. checked in at Bessy II
 *
 * Revision 1.4  90/04/22  21:11:16  heiner
 * edited for indentation style
 * 
 * Revision 1.3  89/07/02  22:57:38  heiner
 * Added min/max for type (char).
 * 
 * Revision 1.2  89/07/02  21:46:11  heiner
 * Extended & polished.
 * 
 * Revision 1.1  89/04/22  00:35:51  heiner
 * Initial revision
 * 
 */

#define USE_PERROR		/* defined when we have perror(3) */
#define TRY_ASCII		/* defined when we try to set PRF_ASCII */

static char	*myname;

#include <stdio.h>

    static void
exi( xcode )
    int	xcode;
{
    (void) fflush(stderr);
    exit(xcode);
}

    static void
err( xcode, str )
    int		 xcode;
    char	*str;
{
    fprintf(stderr, "(%s:) %s", myname, str);
#ifdef USE_PERROR
    if( xcode == 1 ) {
	perror("");
    }else
#endif
    if( xcode ) {
	fprintf(stderr, "\n");
    }
    (void) fflush(stderr);
    if( xcode ) {
	exi(xcode);
    }
}

    static void
do_flush( fp )
    FILE	*fp;
{
    if( ferror(fp) || (EOF == fflush(fp)) || ferror(fp) ) {
	err(1, "error writing standard output");
    }
}

    static void
pr_header()
{
    printf("%s", "#if 0\n");
    printf("%s", " WARNING: This is produced automatically by 'mklocal.c'\n");
    printf("%s", " Changes made in here will eventually go away !\n");
    printf("%s", "#endif\n");
    printf("%s", "\n");
}

    static void
def_whether( name, value )
    char	*name;
    int		 value;
{
    if( value != 0 ) {
	printf("#define %s	1\n", name);
    }else {
	printf("#undef  %s\n", name);
    }
}

static int		max_int;
static int		min_int;
static unsigned int	max_uint;
static int		long_bits;
static long		max_long;
static long		min_long;
static unsigned long	max_u2pot;
static unsigned long	max_ulong;
static char		min_char;
static char		max_char;

static int		int_islong;
static int		easy_decdigs;
static int		native_ascii	= 0;
static int		niz_cvtfp	= 0;

    static int
upd_int( d )
    int		d;	/* some delta */
{
    int		change;
    int		t;
    unsigned	ut;

    change = 0;
    t = max_int + d;
    if( t > max_int ) {
	max_int = t;
	change |= 1;
    }
    t = min_int - d;
    if( t < min_int ) {
	min_int = t;
	change |= 1;
    }
    ut = max_uint + d;
    if( ut > max_uint ) {
	max_uint = ut;
	change |= 1;
    }
    t = d + d;
    if( t > d ) {
	change |= upd_int(t);
    }
    return(change);
}

    static void
mk_int()
{
    max_int = 0;
    min_int = 0;
    max_uint = 0;
    while( upd_int(1) != 0 ) {
	/*nothing*/;
    }
}

    static int
upd_char( d )
    char	d;	/* some delta */
{
    char	t;
    int		change;

    change = 0;
    t = max_char + d;
    if( t > max_char ) {
	max_char = t;
	change |= 1;
    }
    t = min_char - d;
    if( t < min_char ) {
	min_char = t;
	change |= 1;
    }
    t = d + d;
    if( t > d ) {
	change |= upd_char(t);
    }
    return(change);
}

    static void
mk_char()
{
    max_char = 0;
    min_char = 0;
    while( upd_char((char)1) != 0 ) {
	/*nothing*/;
    }
}

    static int
upd_long( d )
    long		d;	/* some delta */
{
    int			change;
    long		t;
    unsigned long	ut;

    change = 0;
    t = max_long + d;
    if( t > max_long ) {
	max_long = t;
	change |= 1;
    }
    t = min_long - d;
    if( t < min_long ) {
	min_long = t;
	change |= 1;
    }
    ut = max_ulong + d;
    if( ut > max_ulong ) {
	max_ulong = ut;
	change |= 1;
    }
    t = d + d;
    if( t > d ) {
	change |= upd_long(t);
    }
    return(change);
}

    static void
no_lfit( lval )
    long	lval;
{
    err(0, "long value does not fit into unsigned long and back:");
    fprintf(stderr, "%ld", lval);
    (void) fflush(stderr);
    exit(2);
}

    static void
chk_lfit( lval )
    long	lval;
{
    long		l;
    unsigned long	ul;

    ul = lval;
    l = ul;
    if( (lval != l) || (lval != (long)(unsigned long)lval) ) {
	no_lfit(lval);
    }
}

    static void
mk_long()
{
    /*
     * Find largest unsigned power of two.
     * Count the number of bits this indicates.
     */
    max_u2pot = 1L;
    long_bits = 1;
    while( (max_u2pot + max_u2pot) > max_u2pot ) {
	max_u2pot += max_u2pot;
	long_bits += 1;
    }
    /*
     * Find exact range of long and unsigned long:
     */
    max_long = 0L;
    min_long = 0L;
    max_ulong = 0L;
    while( upd_long(1L) != 0 ) {
	/*nothing*/;
    }
    /*
     * Check that all longs fit into unsigned long and back:
     */
    chk_lfit(max_long);
    chk_lfit(min_long);
}

    static void
chk_binary()
{
    /*
     * Check that we have a binary representation.
     * This is necessary for some bit fiddling (cf. "prf__udigs.c").
     */
    if( max_ulong < max_u2pot ) {
	err(0, "max_ulong < max_u2pot ?? ");
	fprintf(stderr, " %lx %lx\n", max_ulong, max_u2pot);
	exi(2);
    }
    if( (max_ulong - max_u2pot) != (max_u2pot - 1) ) {
	err(2, "this is not a binary machine, sorry");
    }
}

    static void
pr_int()
{
    printf("#define PRF_MININT	(0x%x)\n", min_int);
    printf("#define PRF_MAXINT	(0x%x)\n", max_int);
}

    static void
pr_char()
{
    printf("#define PRF_MINCHAR	(0x%x)\n", min_char);
    printf("#define PRF_MAXCHAR	(0x%x)\n", max_char);
}

    static void
pr_long()
{
    printf("#define PRF_MINLONG	(0x%lxL)\n", min_long);
    printf("#define PRF_MAXLONG	(0x%lxL)\n", max_long);
    printf("#define PRF_LONGBITS	(%d)\n", long_bits);
}

    static int
intislong()
{
    int	sint;
    int	suint;

    sint = sizeof(int);
    suint = sizeof(unsigned int);
    if( (max_long == (long)max_int)
     && (min_long == (long)min_int)
     && (max_ulong == (unsigned long)max_uint)
     && (sint == sizeof(long))
     && (suint == sizeof(unsigned long))
      ) {
	return(1);
    }
    return(0);
}

    static void
mk_intislong()
{
    int_islong = intislong();
}

    static void
pr_intislong()
{
    def_whether("PRF_INTISLONG", int_islong);
}

    static int
dense_codes( str, len )
    char	*str;
    int	 len;
{
    int	 i;

    for( i=0 ; i < len ; ++i ) {
	if( str[i] != (str[0] + i) ) {
	    return(0);			/* not dense */
	}
    }
    return(1);				/* dense */
}

    static void
mk_easy_decdigs()
{
    easy_decdigs = dense_codes("0123456789", 10);
}

    static void
pr_easy_decdigs()
{
    def_whether("PRF_EASYDECDIGS", easy_decdigs);
}

    static void
mk_ascii()
{
#ifndef TRY_ASCII
    native_ascii = 0;
#else	/* def TRY_ASCII */
    static char	asc[] ={
		' ', '!', '"', '#', '$', '%', '&', '\'',
		'(', ')', '*', '+', ',', '-', '.', '/',
		'0', '1', '2', '3', '4', '5', '6', '7',
		'8', '9', ':', ';', '<', '=', '>', '?',
		'@', 'A', 'B', 'C', 'D', 'E', 'F', 'G',
		'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
		'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W',
		'X', 'Y', 'Z', '[', '\\', ']', '^', '_',
		'`', 'a', 'b', 'c', 'd', 'e', 'f', 'g',
		'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o',
		'p', 'q', 'r', 's', 't', 'u', 'v', 'w',
		'x', 'y', 'z', '{', '|', '}', '~'
    };
    int		 i;

    i = (sizeof asc) / sizeof(char);
    native_ascii = 0;
    if( (asc[0] == 32) && dense_codes(asc, i) && (i == (127-32)) ) {
	native_ascii = 1;
    }
#endif	/* def TRY_ASCII */
}

    static void
pr_ascii()
{
    def_whether("PRF_ASCII", native_ascii);
}

/*
 * The following code tests some pointer types, whether a global declaration
 * of a variable of that type is guaranteed to be initialized with a zero
 * pointer.  May or may not be useful.
 */

#if unix
# define BSS_INITS_WITH_ZERO	1
#endif	/* unix */

#if BSS_INITS_WITH_ZERO
# define CHK_PTYP(accum,ptype)	\
	{	static ptype	just_declared;		\
	    (accum) &= (just_declared == (ptype)NULL);	\
	}
#endif	/* BSS_INITS_WITH_ZERO */

/*
 * Repeat some pointer types from the prf package to test their null pointers:
 */
typedef void	(*PrfCvtFPtr)();	/* pointer to conversion function */
typedef int	(*PrfFlsFPtr)();	/* pointer to buffer flush function */
typedef int	(*PrfPutFPtr)();	/* pointer to putchar function */

    static void
mk_nilptr()
{
#if BSS_INITS_WITH_ZERO
    niz_cvtfp = 1;
    CHK_PTYP(niz_cvtfp, PrfCvtFPtr);
#endif	/* BSS_INITS_WITH_ZERO */
}

    static void
pr_nilptr()
{
    def_whether("PRF_NIZ_CFP", niz_cvtfp);
}

main( argc, argv )
    int		  argc;
    char	**argv;
{
    int		ac;

    ac = 0;
    myname = argv[ac++];

    if( ac != argc ) {
	err(2, "too many arguments");
    }

    mk_int();
    mk_char();
    mk_long();
    mk_intislong();
    mk_easy_decdigs();
    mk_ascii();
    mk_nilptr();

    chk_binary();

    pr_header();
    pr_int();
    pr_char();
    pr_long();
    pr_intislong();
    pr_easy_decdigs();
    pr_ascii();
    pr_nilptr();

    do_flush(stdout);
    exit(0);
}
