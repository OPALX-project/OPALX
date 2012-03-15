#!/bin/bash
#
mkdir -p config
aclocal --force
libtoolize --force --copy
automake --force --add-missing --copy
autoheader --force
autoconf --force
autoreconf   
CXX=mpicxx ./configure  --without-wake-test --without-envelope-solver \
--with-classic-includedir=$OPAL_ROOT/classic/5.0/src --with-classic-libdir=$OPAL_ROOT/classic/5.0/src \
--with-doom-includedir=$OPAL_ROOT/doom/ --with-doom-libdir=$OPAL_ROOT/doom/ \
--with-ippl-includedir=$IPPL_ROOT/src --with-ippl-libdir=$IPPL_ROOT/lib/$IPPL_ARCH \
--with-h5part-includedir=$H5Part/src --with-h5part-libdir=$H5Part/src \
--with-hdf5-includedir=$HDF5HOME/include --with-hdf5-libdir=$HDF5HOME/lib \
--with-gsl-includedir=$GSL_PREFIX/include  --with-gsl-libdir=$GSL_PREFIX/lib \
--with-libs="-lz -lm"
make -j 10 
