###########################################################################
#
# The IPPL Framework
# 
# 
###########################################################################

include $(IPPL_ROOT)/lib/Makefile.def

###########################################################################

USER_INCLUDES   	=	 
USER_DEFINES		=
USER_CXXDEFINES		=
USER_CDEFINES		=
CXXOPT			=
COPT			=
USER_CXXFLAGSLIB	=
USER_CXXFLAGSAPP	=
USER_LDFLAGS		=
USER_LIBS 		=

PROGS = t1 t2 t3 t4 t5 t6 t7 t8 

###########################################################################

###
### Rules to make all these tests
###

all: $(PROGS) 

clean: 
	rm -f $(PROGS) $(PROGS:=.o) core *~

depend:
	$(CXX) $(INCLUDES) $(CXXFLAGS)  -M *.cpp > .depend

include .depend

