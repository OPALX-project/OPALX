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

SHELL = /bin/sh

TESTS = volley timesend reduce
SELF_CHECK = $(TESTS)

include $(IPPL_ROOT)/lib/Makefile.app

all: $(TESTS)

self_check: $(SELF_CHECK)
	for i in $(SELF_CHECK) ; do echo $$i ; $$i ; done

clean: 
	rm -f $(TESTS) $(TESTS:=.o) core *~

clobber: clean
	rm -Rf ii_files ti_files .d*

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

