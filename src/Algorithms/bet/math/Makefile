# Makefile for the math libary
# linux
######################################################################
# 
# Last Revision:
# $Id: Makefile 116 2007-05-18 21:36:28Z bakker $

# HP-UX
######################################################################
CC = g++
LD = ld
AR = ar
CFLAGS = -I.. -O3
LFLAGS = 

OBJS   = \
	integrate.o interpol.o root.o functions.o linfit.o rk.o \
	savgol.o sort.o svdfit.o

ERRORH = ../error.h

LIBFILE = libGmath.a

$(LIBFILE): $(OBJS)
	$(AR) -qc $@ $(OBJS) 

functions.o : functions.C functions.h
	$(CC) -c $(CFLAGS) $< -o $@

linfit.o : linfit.C linfit.h $(ERRORH)
	$(CC) -c $(CFLAGS) $< -o $@

interpol.o : interpol.C interpol.h $(ERRORH)
	$(CC) -c $(CFLAGS) $< -o $@

integrate.o : integrate.C integrate.h $(ERRORH)
	$(CC) -c $(CFLAGS) $< -o $@

rk.o : rk.C rk.h $(ERRORH)
	$(CC) -c $(CFLAGS) $< -o $@

root.o : root.C root.h $(ERRORH)
	$(CC) -c $(CFLAGS) $< -o $@

savgol.o : savgol.C savgol.h $(ERRORH)
	$(CC) -c $(CFLAGS) $< -o $@

sort.o : sort.C sort.h $(ERRORH)
	$(CC) -c $(CFLAGS) $< -o $@

svdfit.o : svdfit.C svdfit.h $(ERRORH)
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o *.stackdump *~
	cd test; make clean

empty: clean
	rm -f *.a
	cd test; make empty
