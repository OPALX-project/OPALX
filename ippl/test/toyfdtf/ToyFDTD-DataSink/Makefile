SHELL = /bin/sh

TESTS = ipplToyFDTD2 

OBJECTS = $(ToyFDTD2_OBJECTS)

ToyFDTD2_OBJECTS = ipplToyFDTD2.o\
	DataSink.o

SELF_CHECK = $(TESTS)

include $(IPPL_ROOT)/lib/$(IPPL_ARCH)/Makefile.def

all: ipplToyFDTD2

ipplToyFDTD2: $(OBJECTS)
	$(LINKER) -o $(TESTS) $(LDFLAGS) $(OBJECTS) $(LFLAGS) $(LIBS)

self_check: $(SELF_CHECK)
	for i in $(SELF_CHECK) ; do echo $$i ; $$i ; done

clean:
	rm -f $(TESTS) $(TESTS:=.o) core *~ DataSink.o

clobber: clean
	rm -rf ii_files ti_files .d*

