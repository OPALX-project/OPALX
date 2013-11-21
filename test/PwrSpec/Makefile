###########################################################################
#
# The IPPL Framework
# 
# This program was prepared by the Regents of the University of
# Visit http://people.web.psi.ch/adelmann/ for more details
#
###########################################################################

SHELL = /bin/sh
include $(IPPL_ROOT)/test/Makefile.def

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
	testPwrSpec	

###
### Almost all test codes, even those that don't follow the passed/failed rule
### and that cannot be used as part of the self-check suite.
###

TESTS = \
	testPwrSpec		

###########################################################################

###
### Rules to make all these tests
###

all: $(TESTS) 

test: tests 

tests: $(TESTS)
	for i in $(TESTS) ; do echo $$i ; $$i ; done

clean: 
	rm -f $(TESTS) $(TESTS:=.o) core *~

clobber: clean
	rm -rf ii_files ti_files .d* text.*


###########################################################################
# $RCSfile: Makefile,v $   $Author: adelmann $
# $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:38 $
###########################################################################
