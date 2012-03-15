SHELL = /bin/sh

PROG = cavitytrack
SELF_CHECK = $(PROG)

TRILINOSROOT = /home/geus/install/Trilinos_20050218/MPICH
FEMAXXROOT = /home/geus/svnwork/femaXX_roman/
PBEROOT = /home/geus/svnwork/pbe

# **** COMPILER FLAGS ****
include $(IPPL_ROOT)/lib/Makefile.def
#
#
USER_INCLUDES = -DMPICH_SKIP_MPICXX -DHAVE_CONFIG_H -DnoBINREPART -DNO_USE_PBE  -w \
                -DRLOG_COMPONENT=femaxx \
		-I./BARNES  -I./FieldSolvers -I./ChargedParticles  \
		-I${FEMAXXROOT}/src \
		-I${FEMAXXROOT}/src/fem \
		-I${FEMAXXROOT}/src/tetmesh \
		-I${FEMAXXROOT}/src/colarray \
		-I${TRILINOSROOT}/include \
                -I/usr/local/rlog-1.3.6/include \
		-I/usr/local/superlu_dist-2.0/include \
		-I./FieldSolvers -I./ChargedParticles -I./Collimator -I./Probe -I./Fields -I./Rksuite -I./RK4\
		-I${PBEROOT}/include
#                
#
OBJECTS = 	./CavitySolver.o ./Distribution/ranlib.o

USER_LIBS = 	$(OBJECTS) -L${FEMAXXROOT}/build -L${TRILINOSROOT}/lib \
		-L/usr/local/superlu-2.0/lib64 -L/usr/local/superlu_dist-2.0/lib64 \
		-L/usr/local/zoltan/lib64 -L/opt/acml/gnu64/lib \
		-L/usr/local/rlog-1.3.6/lib \
	  	-lfft -lfemaxx -lml -laztecoo -lepetraext -lzoltan -ltriutils -lteuchos -lamesos \
		-lepetra -lparmetis -lmetis -lsuperlu_dist -lsuperlu -lacml -lrlog -lg2c
# 
$(PROG):	$(OBJECTS)

#
USER_OBJ = 
####################################################################


all:	$(PROG) 
	cd RK4; make; cd ..
	cd Fields; make Fields.o; cd ..
	cd Distribution; make; cd ..
	cd Rksuite; make; cd ..

cavitytrack.o: cavitytrack.cpp cavitytrack.h CavitySolver.h

CavitySolver.o: CavitySolver.cpp
	$(CXX)  $(USER_INCLUDES) -g -c CavitySolver.cpp

clean:
	rm -fr $(PROG) $(PROG:=.o) core *~ ti*
