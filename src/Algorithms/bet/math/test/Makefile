# Makefile for test programs
# linux
######################################################################
CC = g++
CFLAGS = -I.. -I../..
LFLAGS = -L../../libprf
LIBES  = -lprf -lm

ERROR  = ../../error
ERRORC = $(ERROR).C
ERRORF = $(ERROR).C $(ERROR).h

BINFILES = xsort xrk4 xodeint xqromb xrtsafe xsplint xsavgol xsvdfit xlinfit

all: $(BINFILES)

xrk4: xrk4.C bessel.C bessel.h \
	../rk.C ../rk.h $(ERRORF)
	$(CC) $(CFLAGS) $(LFLAGS) ../rk.C $(ERRORC) bessel.C xrk4.C \
	$(LIBES) -o $@

xodeint: xodeint.C bessel.C bessel.h \
	../rk.C ../rk.h $(ERRORF)
	$(CC) $(CFLAGS) $(LFLAGS) ../rk.C $(ERRORC) bessel.C xodeint.C \
	$(LIBES) -o $@

xsort: xsort.C \
	../sort.C ../sort.h $(ERRORF)
	$(CC) $(CFLAGS) $(LFLAGS) $(ERRORC) ../sort.C xsort.C \
	$(LIBES) -o $@

xqromb: xqromb.C \
	../integrate.C ../integrate.h $(ERRORF)
	$(CC) $(CFLAGS) $(LFLAGS) ../rk.C $(ERRORC) ../integrate.C xqromb.C \
	$(LIBES) -o $@

xsplint: xsplint.C \
	../interpol.C ../interpol.h $(ERRORF)
	$(CC) $(CFLAGS) $(LFLAGS) ../rk.C $(ERRORC) ../interpol.C xsplint.C \
	$(LIBES) -o $@

xrtsafe: xrtsafe.C bessel.C bessel.h \
	../root.C ../root.h $(ERRORF) 
	$(CC) $(CFLAGS) $(LFLAGS) ../root.C $(ERRORC) bessel.C xrtsafe.C \
	$(LIBES) -o $@

xsavgol: xsavgol.C \
	../savgol.C ../savgol.h $(ERRORF) 
	$(CC) $(CFLAGS) $(LFLAGS) $(ERRORC) ../savgol.C xsavgol.C \
	$(LIBES) -o $@

xsvdfit: xsvdfit.C \
	../svdfit.C ../svdfit.h $(ERRORF) 
	$(CC) $(CFLAGS) $(LFLAGS) $(ERRORC) ../svdfit.C xsvdfit.C \
	$(LIBES) -o $@

xlinfit: xlinfit.C \
	../linfit.C ../linfit.h $(ERRORF) 
	$(CC) $(CFLAGS) $(LFLAGS) $(ERRORC) ../linfit.C xlinfit.C \
	$(LIBES) -o $@

clean:
	rm -f *.o *~ *.stackdump

empty: clean
	rm -f $(BINFILES)
