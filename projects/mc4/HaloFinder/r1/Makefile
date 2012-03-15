#CC	= g++
CC	= mpicxx

CFLAGS	= -DDEBUG -DID_32 -DPOSVEL_32 -DGRID_32
#CFLAGS	= -DDEBUG -DID_64 -DPOSVEL_64 -DGRID_64
#MPI_INCLUDE = -I$(MPIHOME)/include
#MPI_LIBS = -L$(MPIHOME)/lib64 -lmpi -lmpi_cxx

HDR= \
	Definition.h \
	Partition.h \
	ParticleDistribute.h \
	ParticleExchange.h \
	GridExchange.h \
	InitialExchange.h \
	CosmoHalo.h \
	CosmoHaloFinder.h \
	CosmoHaloFinderP.h \
	HaloProperties.h

OBJ= \
	Partition.o \
	ParticleDistribute.o \
	ParticleExchange.o \
	CosmoHaloFinder.o \
	CosmoHaloFinderP.o \
	HaloProperties.o \
	HaloTestP.o

GOBJ= \
	Partition.o \
	GridExchange.o \
	GridTestP.o

IOBJ= \
	Partition.o \
	ParticleDistribute.o \
	TestExchange.o \
	InitialExchange.o \
	InitTestP.o

GridTestP: $(GOBJ)
	${CC} ${CFLAGS} ${MPI_INCLUDE}  -O2 -o GridTestP $(GOBJ) ${MPI_LIBS}

HaloTest: HaloTest.cxx CosmoHaloFinder.cxx CosmoHaloFinder.h
	${CC} ${CFLAGS} HaloTest.cxx CosmoHaloFinder.cxx -o HaloTest

HaloTestP: $(OBJ)
	${CC} ${CFLAGS} $(MPI_INCLUDE) -o HaloTestP $(OBJ) $(MPI_LIBS)

InitTestP: $(IOBJ)
	${CC} ${CFLAGS} $(MPI_INCLUDE) -o InitTestP $(IOBJ) $(MPI_LIBS)

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
	${CC} ${CFLAGS} -c CosmoHaloFinder.cxx

HaloProperties.o: HaloProperties.cxx $(HDR)
	${CC} ${CFLAGS} -c HaloProperties.cxx

HaloTestP.o: HaloTestP.cxx $(HDR)
	${CC} ${CFLAGS} -c HaloTestP.cxx

InitialExchange.o: InitialExchange.cxx $(HDR)
	${CC} ${CFLAGS} -c InitialExchange.cxx

TestExchange.o: TestExchange.cxx $(HDR)
	${CC} ${CFLAGS} -c TestExchange.cxx

InitTestP.o: InitTestP.cxx $(HDR)
	${CC} ${CFLAGS} -c InitTestP.cxx

GridExchange.o: GridExchange.cxx $(HDR)
	${CC} ${CFLAGS} -c GridExchange.cxx

GridTestP.o: GridTestP.cxx $(HDR)
	${CC} ${CFLAGS} -c GridTestP.cxx

clean:
	rm *.o HaloTest HaloTestP OneToN BinHalos GridTestP InitTestP
