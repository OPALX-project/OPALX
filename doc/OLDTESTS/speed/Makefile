// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by PSI. 
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit http://www.acl.lanl.gov/POOMS for more details
 *
 ***************************************************************************/

###########################################################################
#
# The IPPL Framework
# 
#
# Visit http://people.web.psi.ch/adelmann/ for more details
#
###########################################################################

include $(IPPL_ROOT)/lib/$(IPPL_ARCH)/Makefile.app

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
# $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:58 $
# IPPL_VERSION_ID: $Id: Makefile,v 1.1.1.1 2000/11/30 21:06:58 adelmann Exp $ 
###########################################################################
