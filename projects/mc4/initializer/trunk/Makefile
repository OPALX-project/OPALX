###### Makefile for the initializer code
##
## FFTW libraries need to be explicitely linked.
##	Code uses FFTW 2 or FFTW 3version, and 
##      loading modules 
##      sets up FFTW_HOME environment variable.
##
## FLAGS:
##      -DNOFFTW
##      For tests that do not require FFT, 
##      above flag can be passed. The code can 
##      than be run on a machine without FFTW.
##
##      -DFFTW2
##      If FFTW version 2 is being used, pass this 
##      flag to the compiler
##
##      -DFFTW3
##      If FFTW version 3 is being used, pass this 
##      flag to the compiler
##
##      -DTESTING 
##      non-production flag -- makes code perform 
##      several tests, most importantly
##      it checks for reality of the density field.

H5Part = /home2/ineichen/svnwork/H5Part
#HDF5HOME = /home2/ineichen/svnwork/hdf5-1.8.2/hdf5
FFTW_HOME = /home2/ineichen/extlib/fftw3
EXEC = init
CPP = mpicxx
CC = mpicc
CPPFLAGS = -O3 -funroll-loops -DH5PART -DFFTW3 -DTESTING -DPARALLEL_IO
CFLAGS = -O3 -funroll-loops -std=c99

LDFLAGS = -L$(FFTW_HOME)/lib -L$(H5Part)/src -L$(HDF5HOME)/lib
FFTW_FLAGS = -I$(FFTW_HOME)/include -I$(H5Part)/src -I$(HDF5HOME)/include

## For FFTW2 uncomment this line:
##FFTWLIBS= -lfftw_mpi  -lrfftw_mpi -lrfftw -lfftw
## For FFTW3 uncomment this line:
FFTWLIBS = -lfftw3_mpi -lfftw3 -lH5Part -lhdf5 -lz


########## SHOULD NOT CHANGE ANYTHING BELOW ##############################

CPP_SOURCES = Cosmology.cpp MT_Random.cpp DataBase.cpp \
              Parallelization.cpp Initializer.cpp main.cpp

C_SOURCES = PerfMon.c distribution.c

INCLUDE = Cosmology.h DataBase.h MT_Random.h Parallelization.h \
          PerfMon.h TypesAndDefs.h Initializer.h distribution.h \
          Basedata.h 

CPP_OBJ = $(CPP_SOURCES:.cpp=.o)
C_OBJ = $(C_SOURCES:.c=.o)


default: $(C_SOURCES) $(CPP_SOURCES) $(EXEC)
	rm -rf *.o

tarball:
	tar -zcf source.tar.gz $(C_SOURCES) $(CPP_SOURCES) \
                               $(INCLUDE) README Makefile indat cmb.tf

clean:
	rm -rf *.o $(EXEC) *~

$(EXEC): $(C_OBJ) $(CPP_OBJ)
	$(CPP) $(LDFLAGS) $(C_OBJ) $(CPP_OBJ) -o $@ $(FFTWLIBS)

.cpp.o: $(INCLUDE)
	$(CPP) -c $(CPPFLAGS) $(FFTW_FLAGS) $< -o $@
.c.o: $(INCLUDE)
	$(CC) -c $(CFLAGS) $(FFTW_FLAGS) $< -o $@
