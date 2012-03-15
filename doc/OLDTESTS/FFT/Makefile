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

##########################################################################
#
# Makefile for fft benchmark
#
##########################################################################


USE_PTUSED = 

##########################################################################
### standard Makefile rules and targets for Linac application
##########################################################################

### initialize information about this Makefile environment
SHELL	= /bin/sh
BACK	= `date +%m%d`

### include the Makefile.app file from IPPL, which includes settings on
### how to compile the program and the appropriate suffix rules
include $(IPPL_ROOT)/lib/$(IPPL_ARCH)/Makefile.app

### final collection of all options (do not comment this line out)
USER_CXXFLAGSAPP     = $(USE_PTUSED) -DOLD_FFT

### files used in the Linac application
MAIN		= fftspeed
ALLTARGETS	= $(MAIN) TestRC


##########################################################################
###
### The default rule for Linac (and Linac.o) is sufficient to compile
### the program; it is obtained automatically from the IPPL Makefile.app
###
##########################################################################

all default: $(ALLTARGETS)


clean:
	/bin/rm -f $(ALLTARGETS) *.o core *~

clobber: clean
	/bin/rm -rf ti_files ii_files .d*


##########################################################################
# end of file
##########################################################################


/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

