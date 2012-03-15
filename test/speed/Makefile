###########################################################################
#
# The IPPL Framework
# 
# This program was prepared by the Regents of the University of
# Visit http://people.web.psi.ch/adelmann/ for more details
#
###########################################################################

include $(IPPL_ROOT)/lib/Makefile.app

SHELL = /bin/sh
TESTS = speed1 stencil1

all: $(TESTS)

time: stencil1
	./stencil1 1000000 1 10 > time

clean: 
	rm -f $(TESTS) $(TESTS:=.o) core *~

clobber: clean
	rm -rf ii_files ti_files .d*

###########################################################################
# $RCSfile: Makefile,v $   $Author: adelmann $
# $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:40 $
# IPPL_VERSION_ID: $Id: Makefile,v 1.1.1.1 2003/01/23 07:40:40 adelmann Exp $ 
###########################################################################
