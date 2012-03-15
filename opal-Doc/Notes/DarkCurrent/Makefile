SHELL = /bin/sh

TESTS = SecondaryEmission

SELF_CHECK = $(TESTS)

include $(IPPL_ROOT)/lib/$(IPPL_ARCH)/Makefile.def

USER_INCLUDES = -DH5_USE_16_API  -I/opt/hdf5/hdf5-1.6.7/include 
USER_LDFLAGS = -L/opt/hdf5/hdf5-1.6.7/lib 
USER_LIBS = -lhdf5 -lz  
USER_CXXOPT = -O2 
USER_OPT = -O2
all: $(TESTS) 

self_check: $(SELF_CHECK)
	for i in $(SELF_CHECK) ; do echo $$i ; $$i ; done

clean:
	rm -f $(TESTS) $(TESTS:=.o) core *~

clobber: clean
	rm -rf ii_files ti_files .d*

