#!/bin/bash
#
mkdir -p config
aclocal --force  
libtoolize --force --copy
aclocal
autoconf --force
automake --force --add-missing --copy      
make clean 
CXX=CC ./configure  --host=x86_64-unknown-linux-gnu \
  --with-hdf5-includedir=/apps/hdf5-1.6.5/include \
  --with-fftw3-includedir=/apps/fftw/fftw-3.1.2_gnu3.3_PE1.4.48/include \
  --with-h5part-includedir=$H5PartHOME/src \
  --with-h5part-libdir=$H5PartHOME/src \
  --with-classic-includedir=$CLASSIC_ROOT/src \
  --with-classic-libdir=$CLASSIC_ROOT/src \
  --with-doom-includedir=$DOOM_ROOT \
  --with-doom-libdir=$DOOM_ROOT \
  --with-ippl-includedir=$IPPL_ROOT/src \
  --with-ippl-libdir=$IPPL_ROOT/lib/$IPPL_ARCH \
  --with-hdf5-libdir=/apps/hdf5-1.6.5/lib \
  --with-fftw3-libdir=/apps/fftw/fftw-3.1.2_gnu3.3_PE1.4.48/lib 
make -j 10 

