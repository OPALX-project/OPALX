THIS	= LibMake
DESTDIR	= 
DEFS	= 
XDEFS	= 
XXDEFS	= # to be set on make command line
CC	= gcc
CPP	= cpp
CCM	= cc -M
CPLUS	= CC
CPLUSM	= CC -M
INCPATH	= 
CFLAGS	=  $(DEFS) $(XDEFS) $(XXDEFS) $(INCPATH)
DFLAGS	=  $(DEFS) $(XDEFS) $(XXDEFS) $(INCPATH)
LD	= ld
LDFLAGS	= 
LEX	= lex
LFLAGS	= 
LINT	= lint
LIFLAGS	= -axbh
YACC	= yacc
YFLAGS	= 
EQUEL	= equel
QFLAGS	= 
F77	= f77
FFLAGS	= 
AR	= /bin/ar
AS	= as
AS_R	= $(AS) -R
STRIP	= strip
SHELL	= /bin/sh
PATH	= /bin:/usr/bin:/etc:/usr/ucb
INSTALL	= install
IFLAGS	= 
MAKEARGS= 'CC=$(CC)' 'CCM=$(CCM)' \
		  'CPLUS=$(CPLUS)' 'CPLUSM=$(CPLUSM)' 'AS=$(AS)' 'LD=$(LD)' \
		  'DEFS=$(DEFS)' 'DESTDIR=$(DESTDIR)' 'INSTALL=$(INSTALL)' \
		  -$(MAKEFLAGS)

all:	# just to make it the first dependency
FRC:	# dummy target

oAn	= objs
oAp	= objs-p
oLn	= objs.nolibc
oLp	= objs.nolibc-p
oSn	= objs.nostdio
oSp	= objs.nostdio-p

oFs     = clrprf.o clrtprf.o defprf.o deftprf.o dprf.o dprfv.o extprf.o exttprf.o \
	iprf.o iprfv.o lprf.o lprfv.o prf_C.o prf_R.o prf_S.o prf_X.o prf__conv.o \
	prf__cput.o prf__d.o prf__mput.o prf__ncput.o prf__num.o prf__parse.o \
	prf__thectab.o prf__udigs.o prf_b.o prf_c.o prf_d.o prf_m.o prf_o.o prf_p.o \
	prf_plural.o prf_repeat.o prf_s.o prf__stretch.o prf_u.o prf_unknown.o \
	prf_x.o setprf.o settprf.o sprf.o sprfv.o stdprf.o stdtprf.o uprf.o uprfv.o \
	prf__naninf.o prf__efcvt.o prf__e.o prf__f.o prf_e.o prf_f.o prf_g.o \
	prf_eE.o prf_gG.o fltprf.o flttprf.o fprf.o fprfv.o prf.o prfv.o


prflocal.h:	./prflocal.h
./prflocal.h:	mklocal.c
	$(CC) mklocal.c
	rm -f "$@"
	./a.out > "$@"
	rm -f a.out

libprf.a: $(oFs)
	$(AR) cr $@ $(oFs)

libprf.sl: $(oFs)
	ld -b -o $@ $(oFs)

$(oDIR)/clrprf.o: clrprf.c
	$(CC) $(CFLAGS) -c -o $@ clrprf.c

$(oDIR)/clrtprf.o: clrtprf.c
	$(CC) $(CFLAGS) -c -o $@ clrtprf.c

$(oDIR)/defprf.o: defprf.c
	$(CC) $(CFLAGS) -c -o $@ defprf.c

$(oDIR)/deftprf.o: deftprf.c
	$(CC) $(CFLAGS) -c -o $@ deftprf.c

$(oDIR)/dprf.o: dprf.c
	$(CC) $(CFLAGS) -c -o $@ dprf.c

$(oDIR)/dprfv.o: dprfv.c
	$(CC) $(CFLAGS) -c -o $@ dprfv.c

$(oDIR)/extprf.o: extprf.c
	$(CC) $(CFLAGS) -c -o $@ extprf.c

$(oDIR)/exttprf.o: exttprf.c
	$(CC) $(CFLAGS) -c -o $@ exttprf.c

$(oDIR)/iprf.o: iprf.c
	$(CC) $(CFLAGS) -c -o $@ iprf.c

$(oDIR)/iprfv.o: iprfv.c
	$(CC) $(CFLAGS) -c -o $@ iprfv.c

$(oDIR)/lprf.o: lprf.c
	$(CC) $(CFLAGS) -c -o $@ lprf.c

$(oDIR)/lprfv.o: lprfv.c
	$(CC) $(CFLAGS) -c -o $@ lprfv.c

$(oDIR)/prf_C.o: prf_C.c
	$(CC) $(CFLAGS) -c -o $@ prf_C.c

$(oDIR)/prf_R.o: prf_R.c
	$(CC) $(CFLAGS) -c -o $@ prf_R.c

$(oDIR)/prf_S.o: prf_S.c
	$(CC) $(CFLAGS) -c -o $@ prf_S.c

$(oDIR)/prf_X.o: prf_X.c
	$(CC) $(CFLAGS) -c -o $@ prf_X.c

$(oDIR)/prf__conv.o: prf__conv.c
	$(CC) $(CFLAGS) -c -o $@ prf__conv.c

$(oDIR)/prf__cput.o: prf__cput.c
	$(CC) $(CFLAGS) -c -o $@ prf__cput.c

$(oDIR)/prf__d.o: prf__d.c
	$(CC) $(CFLAGS) -c -o $@ prf__d.c

$(oDIR)/prf__mput.o: prf__mput.c
	$(CC) $(CFLAGS) -c -o $@ prf__mput.c

$(oDIR)/prf__ncput.o: prf__ncput.c
	$(CC) $(CFLAGS) -c -o $@ prf__ncput.c

$(oDIR)/prf__num.o: prf__num.c
	$(CC) $(CFLAGS) -c -o $@ prf__num.c

$(oDIR)/prf__parse.o: prf__parse.c
	$(CC) $(CFLAGS) -c -o $@ prf__parse.c

$(oDIR)/prf__thectab.o: prf__thectab.c
	$(CC) $(CFLAGS) -c -o $@ prf__thectab.c

$(oDIR)/prf__udigs.o: prf__udigs.c
	$(CC) $(CFLAGS) -c -o $@ prf__udigs.c

$(oDIR)/prf_b.o: prf_b.c
	$(CC) $(CFLAGS) -c -o $@ prf_b.c

$(oDIR)/prf_c.o: prf_c.c
	$(CC) $(CFLAGS) -c -o $@ prf_c.c

$(oDIR)/prf_d.o: prf_d.c
	$(CC) $(CFLAGS) -c -o $@ prf_d.c

$(oDIR)/prf_m.o: prf_m.c
	$(CC) $(CFLAGS) -c -o $@ prf_m.c

$(oDIR)/prf_o.o: prf_o.c
	$(CC) $(CFLAGS) -c -o $@ prf_o.c

$(oDIR)/prf_p.o: prf_p.c
	$(CC) $(CFLAGS) -c -o $@ prf_p.c

$(oDIR)/prf_plural.o: prf_plural.c
	$(CC) $(CFLAGS) -c -o $@ prf_plural.c

$(oDIR)/prf_repeat.o: prf_repeat.c
	$(CC) $(CFLAGS) -c -o $@ prf_repeat.c

$(oDIR)/prf_s.o: prf_s.c
	$(CC) $(CFLAGS) -c -o $@ prf_s.c

$(oDIR)/prf__stretch.o: prf__stretch.c
	$(CC) $(CFLAGS) -c -o $@ prf__stretch.c

$(oDIR)/prf_u.o: prf_u.c
	$(CC) $(CFLAGS) -c -o $@ prf_u.c

$(oDIR)/prf_unknown.o: prf_unknown.c
	$(CC) $(CFLAGS) -c -o $@ prf_unknown.c

$(oDIR)/prf_x.o: prf_x.c
	$(CC) $(CFLAGS) -c -o $@ prf_x.c

$(oDIR)/setprf.o: setprf.c
	$(CC) $(CFLAGS) -c -o $@ setprf.c

$(oDIR)/settprf.o: settprf.c
	$(CC) $(CFLAGS) -c -o $@ settprf.c

$(oDIR)/sprf.o: sprf.c
	$(CC) $(CFLAGS) -c -o $@ sprf.c

$(oDIR)/sprfv.o: sprfv.c
	$(CC) $(CFLAGS) -c -o $@ sprfv.c

$(oDIR)/stdprf.o: stdprf.c
	$(CC) $(CFLAGS) -c -o $@ stdprf.c

$(oDIR)/stdtprf.o: stdtprf.c
	$(CC) $(CFLAGS) -c -o $@ stdtprf.c

$(oDIR)/uprf.o: uprf.c
	$(CC) $(CFLAGS) -c -o $@ uprf.c

$(oDIR)/uprfv.o: uprfv.c
	$(CC) $(CFLAGS) -c -o $@ uprfv.c

$(oDIR)/prf__naninf.o: prf__naninf.c
	$(CC) $(CFLAGS) -c -o $@ prf__naninf.c

$(oDIR)/prf__efcvt.o: prf__efcvt.c
	$(CC) $(CFLAGS) -c -o $@ prf__efcvt.c

$(oDIR)/prf__e.o: prf__e.c
	$(CC) $(CFLAGS) -c -o $@ prf__e.c

$(oDIR)/prf__f.o: prf__f.c
	$(CC) $(CFLAGS) -c -o $@ prf__f.c

$(oDIR)/prf_e.o: prf_e.c
	$(CC) $(CFLAGS) -c -o $@ prf_e.c

$(oDIR)/prf_f.o: prf_f.c
	$(CC) $(CFLAGS) -c -o $@ prf_f.c

$(oDIR)/prf_g.o: prf_g.c
	$(CC) $(CFLAGS) -c -o $@ prf_g.c

$(oDIR)/prf_eE.o: prf_eE.c
	$(CC) $(CFLAGS) -c -o $@ prf_eE.c

$(oDIR)/prf_gG.o: prf_gG.c
	$(CC) $(CFLAGS) -c -o $@ prf_gG.c

$(oDIR)/fltprf.o: fltprf.c
	$(CC) $(CFLAGS) -c -o $@ fltprf.c

$(oDIR)/flttprf.o: flttprf.c
	$(CC) $(CFLAGS) -c -o $@ flttprf.c

$(oDIR)/fprf.o: fprf.c
	$(CC) $(CFLAGS) -c -o $@ fprf.c

$(oDIR)/fprfv.o: fprfv.c
	$(CC) $(CFLAGS) -c -o $@ fprfv.c

$(oDIR)/prf.o: prf.c
	$(CC) $(CFLAGS) -c -o $@ prf.c

$(oDIR)/prfv.o: prfv.c
	$(CC) $(CFLAGS) -c -o $@ prfv.c

clean:
	rm -f *.o *~

empty: clean
	rm -f *.a
