PROGNAME1=orbit
PROGNAME2=cyclorbit

all:  	ranlib.o	

ranlib.o:	ranlib.cc
	g++ -O3 -c -Wno-deprecated ranlib.cc 

clean:
	rm -rf $(PROGNAME1) $(PROGNAME1).o $(PROGNAME2) $(PROGNAME2).o *~ *.o

