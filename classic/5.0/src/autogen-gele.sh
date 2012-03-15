#!/bin/bash
#
mkdir -p config
aclocal
autoconf 
automake -a -c
make clean   
CXX=CC ./configure --host=x86_64-unknown-linux-gnu \
 --with-fftw3-includedir=/apps/fftw/fftw-3.1.2_gnu3.3_PE1.4.48/include \
 --with-fftw3-libdir=/apps/fftw/fftw-3.1.2_gnu3.3_PE1.4.48/lib \
 --with-ippl-includedir=$IPPL_ROOT/src
make -j 10 

