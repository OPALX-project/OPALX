SHELL = /bin/sh

TESTS = mc4

SELF_CHECK = $(TESTS)

include $(IPPL_ROOT)/lib/$(IPPL_ARCH)/Makefile.debug.def

USER_INCLUDES   	= -w \
			-I./\
			-I./Parser\
			-I./ChargedParticles\
			-I./FieldSolvers

USER_OBJS               = ./Parser/*.o Utilities/*.o
USER_LIBS               = $(USER_OBJS) 

all: 	$(TESTS)

clean:
	rm -f $(TESTS) $(TESTS:=.o) core *~ 

