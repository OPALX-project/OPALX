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

################################################################
# Makefile for dcfdump utility program
################################################################

#
# include IPPL make information
#

include $(IPPL_ROOT)/lib/$(IPPL_ARCH)/Makefile.app


#
# the name of the executable
#

MAIN     = dcfdump
MAIN_CPP = $(MAIN).cpp
MAIN_OBJ = $(MAIN).o


#
# the source files
#

HDR = \
	dcfdump.h \
	dcfdump_template.h

SRC = \
	dcfdump_double_1D.cpp \
	dcfdump_double_2D.cpp \
	dcfdump_double_3D.cpp \
	dcfdump_double_4D.cpp \
	dcfdump_float_1D.cpp \
	dcfdump_float_2D.cpp \
	dcfdump_float_3D.cpp \
	dcfdump_float_4D.cpp \
	dcfdump_int_1D.cpp \
	dcfdump_int_2D.cpp \
	dcfdump_int_3D.cpp \
	dcfdump_int_4D.cpp

OBJ = $(SRC:.cpp=.o)


################################################################
# main targets
################################################################

$(MAIN) all: $(MAIN_OBJ) $(OBJ)
	$(LINKER) $(LDFLAGS) -o ${MAIN} $(MAIN_OBJ) ${OBJ} $(LFLAGS) $(LIBS)

clean:
	rm -f *.o *.a *~ core $(MAIN)
	rm -rf ii_files


################################################################
# extra dependencies
################################################################

$(MAIN_OBJ): $(MAIN_CPP) $(HDR) $(SRC)


/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

