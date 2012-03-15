SHELL = /bin/sh

TESTS = testPwrSpec

TESTSFAIL = spatial

SELF_CHECK = $(TESTS)

include ../Makefile.def

all: $(TESTS)

self_check: $(SELF_CHECK)
	for i in $(SELF_CHECK) ; do echo $$i ; $$i ; done

clean:
	rm -f $(TESTS) $(TESTS:=.o) core *~

clobber: clean
	rm -rf ii_files ti_files .d*

