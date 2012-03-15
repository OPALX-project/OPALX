#!/bin/bash
#
mkdir -p config 
aclocal --force
libtoolize --force --copy
automake --force --add-missing --copy
autoheader --force
autoconf --force
autoreconf
CXX=mpicxx ./configure --with-ippl-includedir=$HOME/svnwork/ippl/src --with-ippl-libdir=$HOME/svnwork/ippl/lib/LINUX --with-h5part-includedir=$H5Part/src --with-h5part-libdir=$H5Part/src \
--with-hdf5-includedir=$HDF5HOME/include --with-hdf5-libdir=$HDF5HOME/lib --with-libs="-lz -lm /home2/ineichen/extlib/fftw3/lib/libfftw3_mpi.a /home2/ineichen/extlib/fftw3/lib/libfftw3.a" --with-fftw3-libdir=/home2/ineichen/extlib/fftw3/lib --with-fftw3-includedir=/home2/ineichen/extlib/fftw3/include
#
make -j 10


