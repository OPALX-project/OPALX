CC	= g++
#CC	= mpicxx

CFLAGS	= -DDEBUG
MPI_INCLUDE = -I/opt/OpenMPI/openmpi-1.2.4-gcc/ib/include
MPI_LIBS = -L/opt/OpenMPI/openmpi-1.2.4-gcc/ib/lib64 -lmpi -lmpi_cxx

HDR= \
	Definition.h \
	Partition.h \
	CosmoHalo.h \
	CosmoHaloFinder.h \
	CosmoHaloFinderP.h \
	ParticleDistribute.h \
	ParticleExchange.h

OBJ= \
	Partition.o \
	CosmoHaloFinder.o \
	CosmoHaloFinderP.o \
	ParticleDistribute.o \
	ParticleExchange.o \
	HaloTestP.o

HaloTest: HaloTest.cxx CosmoHaloFinder.cxx CosmoHaloFinder.h
	${CC} -O2 HaloTest.cxx CosmoHaloFinder.cxx -o HaloTest

HaloTestP: $(OBJ)
	${CC} ${CFLAGS} $(MPI_INCLUDE) -o HaloTestP $(OBJ) $(MPI_LIBS)

test: HaloTest
	./HaloTest
	diff test.out test.out.binary.REF | wc -l

OneToN: OneToN.cxx
	${CC} -o OneToN OneToN.cxx

BinHalos: BinHalos.cxx
	${CC} -o BinHalos BinHalos.cxx -lm

Partition.o: Partition.cxx $(HDR)
	${CC} ${CFLAGS} -c Partition.cxx

ParticleDistribute.o: ParticleDistribute.cxx $(HDR)
	${CC} ${CFLAGS} -c ParticleDistribute.cxx

ParticleExchange.o: ParticleExchange.cxx $(HDR)
	${CC} ${CFLAGS} -c ParticleExchange.cxx

CosmoHaloFinderP.o: CosmoHaloFinderP.cxx $(HDR)
	${CC} ${CFLAGS} -c CosmoHaloFinderP.cxx

CosmoHaloFinder.o: CosmoHaloFinder.cxx $(HDR)
	${CC} -c CosmoHaloFinder.cxx

HaloTestP.o: HaloTestP.cxx $(HDR)
	${CC} -c HaloTestP.cxx

clean:
	rm *.o HaloTest HaloTestP OneToN BinHalos
