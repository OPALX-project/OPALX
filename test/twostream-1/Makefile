SHELL = /bin/sh

TESTS = twostream-2

SELF_CHECK = $(TESTS)

include $(IPPL_ROOT)/lib/Makefile.def

USER_INCLUDES   	=-D__cplusplus -I. -IUtility\
			-I${GTHOME}/GTConst\
			-I${GTHOME}/GTElCl\
			-I${GTHOME}/GTUtilities\
			-I${GTHOME}/GTElements\
			-I${GTHOME}/GTDistribution\
			-I${GTHOME}/GTConfigure\
			-I${GTHOME}/GTIntegrators\
			-I${GTHOME}/GTUtilities\
			-I${GTHOME}/GTChargedParticles\
			-I${CMEEHOME}/src

USER_OBJS               = ${GTHOME}/GTDistribution/ranlib.o Utility/*.o
USER_LIBS               = -L${CMEEHOME}/src/.libs -lcmee $(USER_OBJS) 

all: $(TESTS)

self_check: $(SELF_CHECK)
	for i in $(SELF_CHECK) ; do echo $$i ; $$i ; done

objs:	
	g++ -I./Utility -c  ./Utility/Algorithms/PartData.cpp  -o ./Utility/PartData.o
	g++ -I./Utility -c  ./Utility/Parser/AbsFileStream.cpp -o ./Utility/AbsFileStream.o
	g++ -I./Utility -c  ./Utility/Parser/FileStream.cpp     -o ./Utility/FileStream.o
	g++ -I./Utility -c  ./Utility/Parser/StringStream.cpp  -o ./Utility/StringStream.o
	g++ -I./Utility -c  ./Utility/Parser/TerminalStream.cpp  -o ./Utility/TerminalStream.o
	g++ -I./Utility -c  ./Utility/Parser/TokenStream.cpp  -o ./Utility/TokenStream.o
	g++ -I./Utility -c  ./Utility/Parser/Token.cpp  -o ./Utility/Token.o
	g++ -I./Utility -c  ./Utility/Physics/Physics.cpp      -o ./Utility/Physics.o
	g++ -I./Utility -c  ./Utility/Utilities/ClassicException.cpp -o ./Utility/ClassicException.o		
	g++ -I./Utility -c  ./Utility/Utilities/FormatError.cpp -o ./Utility/FormatErroe.o		
	g++ -I./Utility -c  ./Utility/Utilities/LogicalError.cpp -o ./Utility/LogicalError.o		
	g++ -I./Utility -c  ./Utility/Utilities/ParseError.cpp -o ./Utility/ParseError.o		

clean:
	rm -f $(TESTS) $(TESTS:=.o) core *~ ./Utility/*.o


#ClassicException.o
#FormatError.o
#LogicalError.o
#ParseError.o

#SimpleStatement.o
#Statement.o

