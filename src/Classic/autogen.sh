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
    --with-ippl-includedir=$IPPL_ROOT/src \
    --with-h5part-includedir=$H5Part/src --with-h5part-libdir=$H5Part/src \
    --with-gsl-includedir=$GSL_PREFIX/include --with-gsl-libdir=$GSL_PREFIX/lib \
    --with-hdf5-includedir=$HDF5_INCLUDE_PATH --with-hdf5-libdir=$HDF5_LIBRARY_PATH
make -j 10
