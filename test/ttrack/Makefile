SHELL = /bin/sh

TESTS = dsmc 

SELF_CHECK = $(TESTS)

# include $(IPPL_ROOT)/lib/Makefile.def

include ./Makefile.def

USER_INCLUDES   	= -DTTProf -DOLDLAYOUT -DMAPS -DTTRACK -D__cplusplus -I. -I${GTCLASSIC}\
			-I${HOME}/svnwork/gsl/include\
			-I${GTHOME}/GTConst\
			-I${GTHOME}/GTUtilities\
			-I./GTElements\
			-I./GTDistribution\
			-I./GTConfigure\
			-I${GTHOME}/GTIntegrators\
			-I${GTHOME}/GTUtilities\
			-I./GTChargedParticles\
			-I${GTHOME}/GTDataSink\
			-I${GTHOME}/GTQRSolver

USER_OBJS               = ${GTHOME}/GTDistribution/ranlib.o  ${GTHOME}/GTUtilities/*.o \
			./GTDistribution/hammersley.o
USER_LIBS               = $(USER_OBJS) -L${HOME}/svnwork/gsl/lib -lgsl -lgslcblas -L${GTCLASSIC} -lCLASSIC

all: 	$(TESTS)

clean:
	rm -f $(TESTS) $(TESTS:=.o) core *~ 

