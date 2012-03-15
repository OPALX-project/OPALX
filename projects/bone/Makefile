SHELL = /bin/sh

TESTS = typtest PIC2d PIC3d interpolate loadbalance pbconds pbconds3D parexp \
	pdbtest uniform PIC3dBCo-p

TESTSFAIL = spatial

SELF_CHECK = $(TESTS)

include $(IPPL_ROOT)/lib/$(IPPL_ARCH)/Makefile.def

all: $(TESTS)

self_check: $(SELF_CHECK)
	for i in $(SELF_CHECK) ; do echo $$i ; $$i ; done

clean:
	rm -f $(TESTS) $(TESTS:=.o) core *~

clobber: clean
	rm -rf ii_files ti_files .d*

