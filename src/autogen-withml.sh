#!/bin/bash
#
mkdir -p config 
aclocal --force
libtoolize --force --copy
automake --force --add-missing --copy
autoheader --force
autoconf --force
autoreconf
CXX=mpicxx ./configure \
    --with-classic-includedir=$CLASSIC_ROOT/src --with-classic-libdir=$CLASSIC_ROOT/src \
    --with-doom-includedir=$DOOM_ROOT --with-doom-libdir=$DOOM_ROOT \
    --with-ippl-includedir=$IPPL_ROOT/src --with-ippl-libdir=$IPPL_ROOT/lib/$IPPL_ARCH \
    --with-h5part-includedir=$H5Part/src --with-h5part-libdir=$H5Part/src \
    --with-hdf5-includedir=$HDF5_INCLUDE_PATH --with-hdf5-libdir=$HDF5_LIBRARY_PATH \
    --with-gsl-includedir=$GSL_PREFIX/include --with-gsl-libdir=$GSL_PREFIX/lib \
    --with-libdir="-L/opt/parmetis/parmetis-3.1 -L/opt/intel-mkl/mkl-10.0/lib/em64t -L/opt/intel/intel-10.0/fce-10.0/lib" \
    --with-libs="-lz -lm -lifcore -lzoltan -lparmetis -lmetis" \
    --with-blas=mkl --with-lapack=mkl \
    --with-trilinos-includedir=$TRILINOS_INCLUDE_PATH --with-trilinos-libdir=$TRILINOS_LIBRARY_PATH \
    --enable-ml-solver 

make -j 4

