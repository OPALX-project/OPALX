#!/bin/bash
#
mkdir -p config
aclocal
autoheader
autoconf 
automake -a -c
autoreconf
make clean   
CXX=CC ./configure --host=x86_64-unknown-linux-gnu CPP="CC" \
 --with-fftw3-includedir=/apps/fftw/fftw-3.1.2_gnu3.3_PE1.5.47/include \
 --with-fftw3-libdir=/apps/fftw/fftw-3.1.2_gnu3.3_PE1.5.47/lib \
 --with-ippl-includedir=$IPPL_ROOT/src \
 --with-hdf5-includedir=$H5HOME/include --with-hdf5-libdir=$H5HOME/lib \
 --with-h5part-includedir=$H5PartHOME/src --with-h5part-libdir=$H5PartHOME/src
make -j 10 
