###########################################################################
#
# The IPPL Framework
# 
# This program was prepared by the Regents of the University of
# Visit http://people.web.psi.ch/adelmann/ for more details
#
###########################################################################

SHELL = /bin/sh
include $(IPPL_ROOT)/lib/Makefile.def

###########################################################################

###
### self-check programs: these test codes can be run as part of a self-check
### suite.  Each prints out text that contains the word "passed" or "failed".
### If any part of the test code output contains "failed" (case-insensitive)
### the program should be considered as having failed due to some logical
### consistency check problem.  If the word failed does not appear anywhere
### in the output, the program should be considered as not showing a problem.
###

SELF_CHECK = \
	TestPeriodicBC \
	TestAverageVC_CV \
	TestBC \
	TestBC2 \
	TestCartesian \
	TestCartesianCentering \
	TestCrossProduct \
	TestEureka \
	TestFieldDebug \
	TestFieldNDIndexing \
	TestParticleDebug \
	TestScalarIndexing \
	TestTensorUnaryOps \
	TestUniformCartesian \
	balance \
	#binarywhere \
	component \
	compressed \
	drive1 \
	drive2 \
	drive3 \
	patches \
	reduceloc \
	single \
	slices \
	soffset \
	sumslice \
	transpose \
	tz

###
### Almost all test codes, even those that don't follow the passed/failed rule
### and that cannot be used as part of the self-check suite.
###

TESTS = \
	bc \
	bc3 \
	bool \
	doof2d_file \
	float \
	funcface \
	minmax \
	p1 \
	p3 \
	t \
	t2 \
	t3 \
	t4 \
	t5 \
	t6 \
	t7 \
	$(SELF_CHECK) 

###
### Yet more test programs.  Many of these requrie input or graphics.
###

MORETESTS = \
	FVecOps \
	Slice \
	TestLayout \
	bc2 \
	#binary \
	ctoc \
	cv \
	doof2d \
	doof2d_a \
	random \
	subdivide \
	tindex \
	vc

###
### Tests for PETE usage in r1.
###

PETETESTS = \
	PETE_test0d \
	PETE_test1d \
	PETE_test2d \
	PETE_test3d 

###########################################################################

###
### Rules to make all these tests
###

all: $(TESTS) $(MORETESTS) $(PETETESTS)

test: tests moretests

tests: $(TESTS)
	for i in $(TESTS) ; do echo $$i ; $$i ; done

moretests: $(MORETESTS)
	for i in $(MORETESTS) ; do echo $$i ; $$i ; done

petetests: $(PETETESTS)
	for i in $(PETETESTS) ; do echo $$i ; $$i ; done


self_check_build: $(SELF_CHECK)

self_check_run:
	for i in $(SELF_CHECK) ; do \
		echo "=========="; \
		echo $$i ; \
		echo "=========="; \
		./$$i ; \
	done

self_check_run_aclmpl:
	for i in $(SELF_CHECK) ; do \
		echo "=========="; \
		echo $$i ; \
		echo "=========="; \
		./$$i --commlib aclmpl --processes 2; \
	done

self_check: self_check_build
	for i in $(SELF_CHECK) ; do \
		echo "=========="; \
		echo $$i ; \
		echo "=========="; \
		./$$i ; \
	done


clean: 
	rm -f $(TESTS) $(TESTS:=.o) $(MORETESTS) $(MORETESTS:=.o) $(PETETESTS) $(PETETESTS:=.o) core *~

clobber: clean
	rm -rf ii_files ti_files .d* text.*


###########################################################################
# $RCSfile: Makefile,v $   $Author: adelmann $
# $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:38 $
###########################################################################
