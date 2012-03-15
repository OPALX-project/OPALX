
##########################################################################
#
# Makefile for fft benchmark
#
##########################################################################


USE_PTUSED = 

### initialize information about this Makefile environment
SHELL	= /bin/sh
BACK	= `date +%m%d`

### include the Makefile.app file from IPPL, which includes settings on
### how to compile the program and the appropriate suffix rules
include $(IPPL_ROOT)/lib/Makefile.def

### final collection of all options (do not comment this line out)
USER_CXXFLAGSAPP     = $(USE_PTUSED)

### files used in the Linac application
MAIN		= TestFFT 
ALLTARGETS	= $(MAIN)


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
