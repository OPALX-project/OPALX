###### Makefile for the initializer code
##
## FFTW libraries need to be explicitely linked.
##	     Code uses FFTW 2 or FFTW 3version, and 
##      loading modules sets up FFTW_HOME environment
##      variable. If modules are not used, you have to 
##      sets FFTW_HOME to wherever your FFTW resides.
##
## FLAGS:
##      -DLONG_INTEGER
##         If you want to use long integers -- this 
##         may be needed for particle IDs
##
##      -DDOUBLE_REAL
##         If you want positions and velocities in 
##         double precision
##
##      -DNOFFTW
##         For tests that do not require FFT, 
##         above flag can be passed. The code can 
##         than be run on a machine without FFTW.
##
##      -DFFTW2
##         If FFTW version 2 is being used, pass this 
##         flag to the compiler
##
##      -DFFTW3
##         If FFTW version 3 is being used, pass this 
##         flag to the compiler
##
##      -DTESTING 
##         Non-production flag -- makes code perform 
##         several tests, most importantly
##         it checks for reality of the density field.
##
##      -DH5PART
##         For output in HDF5 format
##
##      -DPARALLEL_IO
##         For using parallel IO in HDF5


PROGS = init
CPP = mpicxx
CC = mpicc
CPPFLAGS = -O3 -std=c++0x -funroll-loops -DNOFFTW -DUSENAMESPACE
CFLAGS = -O3 -funroll-loops -std=c99


########## SHOULD NOT CHANGE ANYTHING BELOW ##############################


## Flags:
ifeq (FFTW2,$(findstring FFTW2,$(CPPFLAGS)))  # FFTW2
  FFTWLIBS= -lfftw_mpi  -lrfftw_mpi -lrfftw -lfftw
  LDFLAGS = -L$(FFTW_HOME)/lib
  FFTW_FLAGS = -I$(FFTW_HOME)/include
else
ifeq (FFTW3,$(findstring FFTW3,$(CPPFLAGS)))  # FFTW3
  FFTWLIBS = -lfftw3_mpi -lfftw3
  LDFLAGS = -L$(FFTW_HOME)/lib
  FFTW_FLAGS = -I$(FFTW_HOME)/include
else                                          # No FFTW
  FFTWLIBS = 
endif
endif

ifeq (H5PART,$(findstring H5PART,$(CPPFLAGS))) # HDF5
  LDFLAGS += -L$(H5Part)/src -L$(HDF5HOME)/lib
  FFTW_FLAGS += -I$(H5Part)/src -I$(HDF5HOME)/include
  FFTWLIBS += -lhdf5 -lH5Part -lz
endif
## Compiling:

CPP_SOURCES = Cosmology.cpp MT_Random.cpp DataBase.cpp InputParser.cpp\
              Parallelization.cpp Output.cpp Initializer.cpp main.cpp

C_SOURCES = PerfMon.c distribution.c

INCLUDE = Cosmology.h DataBase.h MT_Random.h Parallelization.h InputParser.h\
          PerfMon.h TypesAndDefs.h Initializer.h distribution.h \
          Output.h 

CPP_OBJ = $(CPP_SOURCES:.cpp=.o)
C_OBJ = $(C_SOURCES:.c=.o)


default: $(C_SOURCES) $(CPP_SOURCES) $(PROGS)

tarball:
	tar -zcf source.tar.gz $(C_SOURCES) $(CPP_SOURCES) \
                               $(INCLUDE) README Makefile input.par cmb.tf
lib:	
	ar -r libMC4Initializer.a PerfMon.o distribution.o Cosmology.o MT_Random.o DataBase.o InputParser.o Parallelization.o Output.o Initializer.o

clean:
	rm -rf *.o $(PROGS)

init: $(C_OBJ) $(CPP_OBJ)
	$(CPP) $(LDFLAGS) $(C_OBJ) $(CPP_OBJ) -o $@ $(FFTWLIBS)

.cpp.o: $(INCLUDE)
	$(CPP) -c $(CPPFLAGS) $(FFTW_FLAGS) $< -o $@

.c.o: $(INCLUDE)
	$(CC) -c $(CFLAGS) $(FFTW_FLAGS) $< -o $@
