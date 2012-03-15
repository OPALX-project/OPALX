###########################################################################
#
# The IPL Framework
# 
# 
###########################################################################

include $(IPL_ROOT)/lib/$(IPL_ARCH)/Makefile.def

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

PROGS = param2proc

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

