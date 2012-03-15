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
    --with-ippl-includedir=$IPPL_ROOT/src --with-ippl-libdir=$IPPL_ROOT/lib/$IPPL_ARCH \
    --with-h5part-includedir=$H5Part/src --with-h5part-libdir=$H5Part/src \
    --with-hdf5-includedir=$HDF5_INCLUDE_PATH --with-hdf5-libdir=$HDF5_LIBRARY_PATH \
    --with-gsl-includedir=$GSL_INCLUDE --with-gsl-libdir=$GSL_PREFIX/lib \
    --with-libs="-lz -lm"
make -j 4

