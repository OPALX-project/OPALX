SHELL = /bin/sh

TESTS = hdf5write-1

SELF_CHECK = $(TESTS)

include $(IPPL_ROOT)/lib/Makefile.def
USER_INCLUDES = 
all: $(TESTS)

self_check: $(SELF_CHECK)
	for i in $(SELF_CHECK) ; do echo $$i ; $$i ; done

clean:
	rm -f $(TESTS) $(TESTS:=.o) core *~

clobber: clean
	rm -rf ii_files ti_files .d*

